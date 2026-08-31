"""
Radiation band groups and frequency-integrated rate coefficients.

This module defines two classes:

- :class:`RadiationGroup` -- a single frequency band (energy interval
  ``[lower, upper]`` in eV) that holds per-reaction rate coefficients and
  cross-section data.
- :class:`Radiation` -- the full collection of bands; responsible for
  computing rate coefficients by integrating tabulated photo cross
  sections (and user-supplied ``dRad`` functions) over each band using a
  power-law photon-number spectrum.

Photon spectrum assumption
--------------------------
The photon number spectrum is assumed to follow a power law in photon energy::

    n(E) ∝ E^(α - 2)

where ``α`` is ``powerlaw_idx``.  The energy-integrated version (energy
density per unit energy interval) is therefore::

    u(E) = E * n(E) ∝ E^(α - 1)

This form is used when computing band-average cross sections and average
photon energies.

Rate coefficient derivation
----------------------------
For a reaction with tabulated cross section σ(E) the *photon-flux-weighted*
average cross section in band *i* is::

    <σ>_i = ∫_{E_lo}^{E_hi} σ(E) n(E) dE  /  ∫_{E_lo}^{E_hi} n(E) dE

The symbolic rate coefficient stored in the radiation density variable
``den[i]`` (either energy density *u_i* in erg/cm³ or photon density *n_i*
in cm⁻³) is::

    k_i = c * den[i] * <σ>_i          (photon density mode)
    k_i = c * den[i] * <σ>_i / <E>_i  (energy density mode)

where *c* is the speed of light and ``<E>_i`` is the band-average photon
energy (``eavg``).

Energy units
------------
Band edges, the photon-energy symbol ``E`` and the cross-section tables are
all in **eV**, and the band-average cross sections are energy-unit-free ratios.
The one quantity carrying a net energy dimension, the band-average photon
energy ``eavg`` (``<E>_i``), is converted to **erg** (see :data:`EV_TO_ERG`)
so the rate coefficients and radiation-moment ODEs are consistent with the CGS
solver; the ``radeden`` field is therefore an energy density in erg/cm³.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, cast

import numpy as np
import sympy as sp
from astropy import units as u

from ...common._integrators import arr_integrate, smart_integrate
from .._typing import RadiationGroupReactionProps
from ._photochemistry import Photochemistry
from ._typing._photochemistry import XsecsProps
from .background_field import BackgroundField

if TYPE_CHECKING:
    from ...core.network import Network
    from ...core.reaction import Reaction


class RadiationGroup:
    """
    A single frequency band in the radiation field discretisation.

    Each group spans the photon-energy interval ``[lower, upper]`` (in eV)
    and accumulates per-reaction rate-coefficient data populated by
    :meth:`Radiation.set_reaction_rate_coefficient` and
    :meth:`Radiation.set_custom_rate`.

    Parameters
    ----------
    lower : float or int
        Lower bound of the energy band in eV.
    upper : float, int, or sympy.Basic
        Upper bound of the energy band in eV.  May be ``sympy.oo`` for the
        uppermost open band.
    index : int
        Zero-based position of this group in the parent :class:`Radiation`
        group list.

    Attributes
    ----------
    index : int
        Band index (same as the constructor argument).
    lower : float or int
        Lower energy bound in eV.
    upper : float, int, or sympy.Basic
        Upper energy bound in eV.
    band : tuple
        ``(lower, upper)`` convenience pair.
    dE : float or sympy.Basic
        Band width ``upper - lower`` in eV.
    props : dict
        Mapping from :class:`~jaff.core.reaction.Reaction` objects to a
        :class:`~jaff.physics._typing.RadiationGroupReactionProps` dict with
        keys:

        - ``"k"``         : symbolic rate coefficient (SymPy ``Expr``) for
          this band.
        - ``"xsec"``      : photon-number-weighted band-average cross section
          (cm²), or ``None`` for custom-rate reactions.
        - ``"xsec_frac"`` : fraction of the total cross section (or total
          ``dRad``) attributed to this band (dimensionless).
        - ``"delta_rad"`` : integrated ``dRad`` over the band -- the radiation
          energy added to the band per reaction event (erg).

    eavg : float or None
        Photon-number-weighted average energy of this band, in **erg** (the
        band-edge integral is in eV and converted via :data:`EV_TO_ERG`), so
        that dividing rate/ODE terms by it stays CGS-consistent.  Computed in
        :class:`Radiation.__init__` and shared across all reactions in the band.
    """

    def __init__(
        self, lower: float | int, upper: float | int | sp.Basic, index: int, sym: sp.Basic
    ):
        """Initialise a single radiation band.

        Parameters
        ----------
        lower : float or int
            Lower energy bound of the band in eV.
        upper : float, int, or sympy.Basic
            Upper energy bound in eV.  May be ``sympy.oo`` for an open band.
        index : int
            Zero-based position of this group in the parent :class:`Radiation`
            group list.
        """
        self.index: int = index
        self.sym: sp.Basic = sym
        self.lower: float | int = lower
        self.upper: float | int | sp.Basic = upper
        self.band: tuple = (lower, upper)

        if not isinstance(self.lower, (int, float)):
            raise ValueError(
                f"Radiation group lower bound must be a float or int. Found {self.lower}"
            )
        # Band width; may be symbolic when upper is sp.oo.
        self.dE: float | None = (
            self.upper - self.lower  # type: ignore
            if all(isinstance(val, (int, float)) for val in [self.upper, self.lower])
            else None
        )
        self.props: dict[Reaction, RadiationGroupReactionProps] = {}
        # Populated on the first call to set_reaction_rate_coefficient for this band.
        self.eavg: float | None = None

    def __repr__(self):
        """Return detailed string representation of this radiation group.

        Returns
        -------
        str
            String including group index and energy band.
        """
        return f"Rad_group({self.index}, band={self.band})"

    def __str__(self):
        """Return human-readable description of this radiation group.

        Returns
        -------
        str
            String of form ``"Radiation group <index>"``.
        """
        return f"Radiation group {self.index}"


class Radiation:
    """
    Collection of frequency bands with integrated photoionisation rate coefficients.

    On construction the band-edge list is parsed (replacing the string
    ``"inf"`` with ``sympy.oo``) and a :class:`RadiationGroup` object is
    created for each consecutive pair of edges.

    Parameters
    ----------
    bands : list of (int, float, str, or sympy.Basic)
        Ordered list of *N+1* band-edge photon energies in eV, defining *N*
        frequency bands.  The string ``"inf"`` is accepted as the last entry
        to represent an open upper boundary.
    powerlaw_idx : int or float
        Power-law index *α* for the assumed photon-number spectrum
        ``n(E) ∝ E^(α-2)``.  Typical values: 1 (flat energy spectrum),
        0 (flat photon spectrum).
    energy_density : bool
        If ``True`` the radiation field is tracked as energy density
        (erg cm⁻³); if ``False`` as photon number density (cm⁻³).  This
        controls the name of the symbolic density variable (``"radeden"`` vs.
        ``"photden"``) and the normalisation of rate coefficients.
    c : float
        Speed of light in cm/s (CGS).  Used in rate-coefficient expressions
        as ``k = c * den * <σ>``.

    Attributes
    ----------
    bands : list of (int, float, or sympy.Basic)
        Parsed band-edge list (``"inf"`` replaced by ``sympy.oo``).
    powerlaw_idx : int or float
        Power-law index passed to the constructor.
    energy_density : bool
        Whether the field is tracked as energy density or photon density.
    c : float
        Speed of light in cm/s.
    nbands : int
        Number of bands (``len(bands) - 1``).
    groups : list of RadiationGroup
        One :class:`RadiationGroup` per band, in ascending energy order.
    """

    def __init__(
        self,
        network: Network,
        bands: list[int | float | str | sp.Basic],
        powerlaw_idx: int | float,
        energy_density: bool,
        c: float | str,
        background_field: str = "draine",
    ):
        """Parse band edges and construct one :class:`RadiationGroup` per band.

        Parameters
        ----------
        bands : list of int, float, str, or sympy.Basic
            Ordered list of *N+1* photon-energy band edges in eV.  The string
            ``"inf"`` is accepted as the last entry to represent an open upper
            boundary (converted to ``sympy.oo``).
        powerlaw_idx : int or float
            Power-law spectral index *α* for ``n(E) ∝ E^(α-2)``.
        energy_density : bool
            If ``True``, radiation is tracked as energy density (erg cm⁻³);
            if ``False``, as photon number density (cm⁻³).
        c : float | str
            Speed of light in cm/s (CGS) or a string to be converted to a symbol.
        """
        self.network: Network = network
        self.bands: list[int | float | sp.Basic] = []
        self.powerlaw_idx: int | float = powerlaw_idx
        self.energy_density: bool = energy_density
        # Speed of light (cm/s) for k = c * σ * n(E) expressions
        self.c: float | sp.Symbol = sp.symbols(c) if isinstance(c, str) else c
        self.background_field = BackgroundField(background_field)

        self.__parse_bands(bands)
        self.nbands: int = len(self.bands) - 1
        # Symbolic radiation density variable: energy density (erg/cm³) or
        # photon number density (cm⁻³), depending on the mode.
        self.den = sp.MatrixSymbol(
            "radeden" if self.energy_density else "photden", self.nbands, 1
        )
        self.groups: list[RadiationGroup] = [
            RadiationGroup(lower, self.bands[i + 1], i, self.den[sp.Idx(i)])  # type: ignore
            for i, lower in enumerate(self.bands[:-1])
        ]
        self.E_sym: sp.Symbol = sp.Symbol("E")
        self.ph_profile_sym: sp.Expr = self.E_sym ** (self.powerlaw_idx - 2)
        self.energy_profile_sym: sp.Expr = self.E_sym * self.ph_profile_sym

        self.photden_tot = smart_integrate(
            self.ph_profile_sym, self.E_sym, (self.bands[0], self.bands[-1])
        )

        for grp in self.groups:
            # Compute the band-average photon energy once per band (shared
            # across all reactions): <E>_i = ∫ E n(E) dE / ∫ n(E) dE
            grp.eavg = (
                smart_integrate(
                    self.energy_profile_sym, self.E_sym, (grp.lower, grp.upper)
                )
                / self.photden_tot
            ) * u.eV.to(u.erg)

    def set_reaction_rate_coefficient(self, reaction: Reaction) -> None:
        """
        Compute and store symbolic band-averaged rate coefficients for a reaction.

        Reads the tabulated photo cross section σ(E) from
        ``reaction.xsecs_dict``, then for each frequency band:

        1. Computes the photon-number-weighted band-average cross section
           ``<σ>_i = ∫ σ n dE / ∫ n dE``.
        2. Computes the band-average photon energy
           ``<E>_i = ∫ E n dE / ∫ n dE`` (stored once per band).
        3. Assembles the symbolic rate coefficient
           ``k_i = c * den[i] * <σ>_i`` (photon-density mode) or
           ``k_i = c * den[i] * <σ>_i / <E>_i`` (energy-density mode).
        4. Stores ``k_i``, the integrated cross section, the cross-section
           fraction, and the integrated ``dRad`` in
           ``self.groups[i].props[reaction]``.

        After iterating over all bands the total symbolic rate coefficient
        (sum over bands, in units of s⁻¹ or cm³ s⁻¹ depending on reaction
        type) is written to ``reaction.rate``.

        If the reaction has no tabulated cross section (no ``xsecs_dict``, or
        no ``photodecay`` array) the method returns silently (no-op).

        Parameters
        ----------
        reaction : Reaction
            The photochemical reaction to process.  ``reaction.dRad`` must
            be a SymPy expression in the symbol ``E`` (photon energy in eV)
            describing the radiation energy the reaction adds to the field per
            unit photon energy (erg/eV); integrated over each band it gives the
            energy added per reaction event (erg).

        Returns
        -------
        None
            Results are stored in ``reaction.rate``, ``reaction.rad_xsecs``,
            ``self.groups[i].props[reaction]``, and ``reaction.rad_groups``
            (back-references to the bands this reaction contributes to)
            in-place.

        Notes
        -----
        The photon-number spectrum used for averaging is
        ``n(E) ∝ E^(powerlaw_idx - 2)``.  For ``powerlaw_idx = 1`` this
        gives a flat energy spectrum; for ``powerlaw_idx = 0`` a flat photon
        spectrum.

        Cross-section integrals (``∫ σ n dE``) are evaluated numerically by
        :func:`~jaff.common._integrators.arr_integrate` over the tabulated
        ``(E, σ)`` arrays.  The remaining analytic integrals over ``n(E)``,
        ``E n(E)`` and ``reaction.dRad`` use
        :func:`~jaff.common._integrators.smart_integrate`, which falls back
        to numerical quadrature when SymPy cannot find a closed form.
        """
        xsec: XsecsProps | None = reaction.xsecs_dict
        if xsec is None:
            return

        # Each reaction carries a single decay-channel cross section.
        pr_xsec = xsec["photodecay"]
        if pr_xsec is None:
            return

        assert isinstance(xsec["photon_energy"], np.ndarray)

        # Photon-number spectrum: n(E) ∝ E^(α-2) used for weighing the cross-section
        # where α = powerlaw_idx.  The factor E^(α-2) arises from
        # n(E) = u(E)/E and u(E) ∝ E^(α-1).
        E = xsec["photon_energy"]  # photon energy array in eV
        energy_profile = E ** (self.powerlaw_idx - 2)
        k_tot = sp.Float(0.0)  # Accumulates total rate coefficient over all bands

        # Total cross section integrated over the full spectrum (cm²),
        # stored on the reaction for later reference
        xsec_tot = (
            arr_integrate(pr_xsec * energy_profile, E, (self.bands[0], self.bands[-1]))
            / self.photden_tot
        )
        reaction.rad_xsecs = xsec_tot

        # Reset any back-references from a previous run before repopulating.
        reaction.rad_groups = []

        for grp in self.groups:
            # ∫ n(E) dE over the band — used as normalisation for averages.
            photden_band = smart_integrate(
                self.ph_profile_sym, self.E_sym, (grp.lower, grp.upper)
            )

            # Photon-number-weighted average cross section in the band:
            # <σ>_i = ∫ σ(E) n(E) dE / ∫ n(E) dE
            pr_xsec_avg = (
                arr_integrate(pr_xsec * energy_profile, E, (grp.lower, grp.upper))
                / photden_band
            )
            rad_xsec_avg = (
                (
                    arr_integrate(
                        xsec["photo_absorption"] * energy_profile,
                        E,
                        (grp.lower, grp.upper),
                    )
                    / photden_band
                )
                if xsec["_equations"]["pa"]
                else pr_xsec_avg
            )

            # Integral of the user-supplied dRad (radiation energy per photon
            # energy, erg/eV) over the band -> energy per reaction event (erg).
            delta_rad_band = smart_integrate(
                reaction.dRad, self.E_sym, (grp.lower, grp.upper)
            )

            # Symbolic rate coefficient: k_i = c · den[i] · <σ>_i
            # (units: s⁻¹ for photon-density mode, cm³ s⁻¹ for two-body)
            k = self.c * self.den[sp.Idx(grp.index)] * rad_xsec_avg
            if "shielding" in reaction._metadata:
                if "value" in reaction._metadata["shielding"]:
                    k *= reaction._metadata["shielding"]["value"]
                else:
                    k *= Photochemistry.shielding(reaction, self.network)

            grp.props[reaction] = {
                "k": k,
                "xsec": rad_xsec_avg,
                "xsec_frac": rad_xsec_avg / xsec_tot,  # fraction of total cross section
                "delta_rad": delta_rad_band,
            }
            reaction.rad_groups.append(grp)

            # In energy-density mode, convert from "per eV" to "per photon"
            # by dividing by the band-average energy <E>_i.
            k_tot += (
                k
                * (1.0 if not xsec["_equations"]["pa"] else (pr_xsec_avg / rad_xsec_avg))
                / (grp.eavg if self.energy_density else 1)
            )

        reaction.rate = k_tot

    def set_custom_rate(self, reaction: Reaction) -> None:
        """
        Partition a user-supplied reaction rate across frequency bands.

        For reactions whose rate is provided analytically (rather than derived
        from a Verner cross section), this method distributes the total rate
        proportionally to the fraction of the ``dRad`` integral that falls
        in each band::

            k_i = reaction.rate * (∫_{E_lo}^{E_hi} dRaddE) / (∫_all dRad dE)

        Parameters
        ----------
        reaction : Reaction
            Reaction with a pre-assigned ``reaction.rate`` (total rate
            coefficient) and ``reaction.dRad`` as a SymPy expression in
            the symbol ``E`` (photon energy in eV).

        Returns
        -------
        None
            Results are stored in ``self.groups[i].props[reaction]`` and
            ``reaction.rad_groups`` in-place.

        Notes
        -----
        ``dRad`` must be expressed per unit eV and the symbol ``E`` must be
        used as the integration variable.

        If the total ``dRad`` integral evaluates to zero (e.g. the
        reaction has no radiation coupling), all band fractions are set to
        zero to avoid division by zero.
        """
        # E is the photon energy symbol used in dRad expressions.
        E = sp.Symbol("E")
        # Integrate dRad over the full spectrum to use as denominator.
        delta_rad_total = smart_integrate(
            reaction.dRad, E, (self.bands[0], self.bands[-1])
        )
        # Guard against zero-denominator case (no radiation coupling).
        delta_rad_total_is_zero = delta_rad_total == 0.0

        # Reset any back-references from a previous run before repopulating.
        reaction.rad_groups = []

        for grp in self.groups:
            # Band-integrated dRad (numerator of the fraction).
            delta_rad_band = smart_integrate(reaction.dRad, E, (grp.lower, grp.upper))
            # Fraction of the total radiation coupling attributed to this band.
            xsec_frac = (
                0.0 if delta_rad_total_is_zero else delta_rad_band / delta_rad_total
            )
            # Scale the total user-supplied rate by the band fraction.
            k = reaction.rate * xsec_frac

            grp.props[reaction] = {
                "k": k,
                "xsec": None,  # No tabulated cross section for custom reactions
                "xsec_frac": xsec_frac,
                "delta_rad": delta_rad_band,
            }
            reaction.rad_groups.append(grp)

    def ordered_index(self, idx: int, order: int) -> tuple[int, int]:
        """
        Map a band index to positions in the flat radiation ODE output array.

        The radiation field is represented by two quantities per band —
        density (energy or photon) and flux — laid out in a flat array.
        Four layout conventions are supported, selected by *order*.

        Parameters
        ----------
        idx : int
            Zero-based band index (0 ≤ ``idx`` < ``self.nbands``).
        order : {0, 1, 2, 3}
            Output-array layout convention:

            - ``0``: ``[den_0, flux_0, den_1, flux_1, ...]``
              -- density and flux interleaved, density first.
            - ``1``: ``[flux_0, den_0, flux_1, den_1, ...]``
              -- density and flux interleaved, flux first.
            - ``2``: ``[den_0, den_1, ..., flux_0, flux_1, ...]``
              -- all densities in the first half, all fluxes in the second.
            - ``3``: ``[flux_0, flux_1, ..., den_0, den_1, ...]``
              -- all fluxes in the first half, all densities in the second.

        Returns
        -------
        ei : int
            Index in the output array for the *density* (energy/photon)
            variable of band *idx*.
        fi : int
            Index in the output array for the *flux* variable of band *idx*.

        Notes
        -----
        This method is used by :func:`~jaff.physics._equations.get_sradodes`
        to place each band's pair of moment equations at the correct position
        in the generated ODE array, matching the memory layout expected by
        the downstream numerical integrator.
        """
        # Default (order 0): interleaved, density first.
        ei = 2 * idx  # density slot
        fi = 2 * idx + 1  # flux slot

        if order == 1:
            # Interleaved, flux first.
            ei = 2 * idx + 1
            fi = 2 * idx
        elif order == 2:
            # Block layout: all densities, then all fluxes.
            ei = idx
            fi = self.nbands + idx
        elif order == 3:
            # Block layout: all fluxes, then all densities.
            ei = self.nbands + idx
            fi = idx

        return ei, fi

    def __parse_bands(self, bands: list[float | int | str | sp.Basic]):
        """
        Validate and store the band-edge list, replacing ``"inf"`` with ``sympy.oo``.

        Also checks that the power-law photon-number spectrum is integrable
        over the supplied band range when energy-density mode is active.
        The average energy ``<E>_i = ∫ E·n(E) dE / ∫ n(E) dE`` must
        converge; this requires:

        - The lower edge to be non-zero when the spectral index is steep
          enough to cause a divergence at ``E → 0``.
        - The upper edge to be finite when the spectral index is shallow
          enough to cause a divergence at ``E → ∞``.

        Parameters
        ----------
        bands : list of (float, int, str, or sympy.Basic)
            Mutable band-edge list; modified in-place to replace any
            ``"inf"`` string with ``sympy.oo``.

        Raises
        ------
        RuntimeError
            If the average-energy integral would diverge given the supplied
            band edges and power-law index.
        """
        # Replace the sentinel string "inf" with SymPy's infinity symbol.
        if "inf" in bands:
            inf_index = bands.index("inf")
            bands[inf_index] = sp.oo

        self.bands = cast(list[int | float | sp.Basic], bands)

        if self.energy_density:
            # The average-energy integral uses the *energy-density* spectrum
            # u(E) ∝ E^(α-1), so the integral ∫ E · u(E) dE ∝ ∫ E^α dE.
            # The effective power-law index for the ∫ E·n(E) dE integral is
            # pl_index = (α-2) + 1 = α - 1.
            pl_index: float = float(self.powerlaw_idx) - 1.0

            if pl_index == -1.0:
                # Integrand ~ E^(-1): log-divergence at both E=0 and E=∞.
                if (
                    isinstance(self.bands[0], (float, int))
                    and float(self.bands[0]) == 0.0
                ):
                    raise RuntimeError(
                        f"The integral for average energy will diverge since the radiation band starts from bands[0]: {self.bands[0]}\n"
                        "Please try a non-zero value"
                    )
                if self.bands[-1] == sp.oo:
                    raise RuntimeError(
                        f'The integral for average energy will diverge since the radiation band ends at bands[{len(self.bands) - 1}]: "inf"\n'
                        "Please try a non-infinite value or change the power_law_index"
                    )
            elif pl_index + 1.0 > 0.0:
                # Integrand ~ E^p with p > -1: diverges at E → ∞.
                if self.bands[-1] == sp.oo:
                    raise RuntimeError(
                        f'The integral for average energy will diverge since the radiation band ends at bands[{len(self.bands) - 1}]: "inf"\n'
                        "Please try a non-infinite value or change the power_law_index"
                    )
            elif pl_index + 1.0 < 0.0:
                # Integrand ~ E^p with p < -1: diverges at E → 0.
                if (
                    isinstance(self.bands[0], (float, int))
                    and float(self.bands[0]) == 0.0
                ):
                    raise RuntimeError(
                        f"The integral for average energy will diverge since the radiation band starts from bands[0]: {self.bands[0]}\n"
                        "Please try a non-zero value"
                    )
