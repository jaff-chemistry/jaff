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

where ``α`` is the spectral index ``profile_idx`` (from
``RadiationProps.profile_index``).  ``α`` may be a single global value or one
value per band, in which case the spectrum is a piecewise power law (a
histogram, intentionally discontinuous at band edges): each band ``i`` uses its
own ``α_i`` and is normalised independently.  The energy-integrated version
(energy density per unit energy interval) is therefore::

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
energy ``eavg`` (``<E>_i``), is converted from eV to **erg** (via
``astropy.units``, ``u.eV.to(u.erg)``) so the rate coefficients and
radiation-moment ODEs are consistent with the CGS solver; the ``radeden``
field is therefore an energy density in erg/cm³.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import sympy as sp
from astropy import units as u

from ...common._integrators import smart_integrate
from .._typing import RadiationGroupReactionProps
from ._photochemistry import Photochemistry
from ._typing._photochemistry import XsecsProps
from .background_field import BackgroundField

if TYPE_CHECKING:
    from ...core.network import Network
    from ...core.reaction import Reaction
    from ...physics import RadiationProps


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
    sym : sympy.Basic
        Symbolic radiation-density entry for this band (the ``self.den``
        matrix element ``den[index]`` supplied by :class:`Radiation`).
    profile_idx : float
        Spectral index *α* of this band's photon-number spectrum
        ``n(E) ∝ E^(α-2)``.

    Attributes
    ----------
    index : int
        Band index (same as the constructor argument).
    sym : sympy.Basic
        Symbolic radiation-density entry for this band (energy density or
        photon number density depending on the parent mode).
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

    profile_idx : float
        Spectral index *α* of this band.
    nph_profile : sympy.Expr
        Photon-number spectral profile ``E^(profile_idx - 2)`` of this band.
    energy_profile : sympy.Expr
        Energy-density spectral profile ``E * nph_profile`` =
        ``E^(profile_idx - 1)`` of this band.
    photden : sympy.Expr or float
        ``∫ nph_profile dE`` over the band, the normalisation for band
        averages.
    eavg : float
        Photon-number-weighted average energy of this band, in **erg** (the
        band-edge integral is in eV and converted via
        ``astropy.units``, ``u.eV.to(u.erg)``), so that dividing rate/ODE
        terms by it stays CGS-consistent.  Computed on construction and
        shared across all reactions in the band.
    """

    E_sym: sp.Symbol = sp.Symbol("E")

    def __init__(
        self,
        lower: float | int,
        upper: float | int | sp.Basic,
        index: int,
        sym: sp.Basic,
        profile_idx: float,
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
        sym : sympy.Basic
            Symbolic radiation-density entry for this band (the ``den[index]``
            matrix element supplied by :class:`Radiation`).
        profile_idx : float
            Spectral index *α* of this band's photon-number spectrum
            ``n(E) ∝ E^(α-2)``.
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
        self.profile_idx: float = profile_idx
        self.nph_profile: sp.Expr = self.E_sym ** (self.profile_idx - 2)
        self.energy_profile: sp.Expr = self.E_sym * self.nph_profile
        # ∫ n(E) dE over the band — used as normalisation for averages.
        self.photden = smart_integrate(
            self.nph_profile, self.E_sym, (self.lower, self.upper)
        )
        # Compute the band-average photon energy once per band (shared
        # across all reactions): <E>_i = ∫ E n(E) dE / ∫ n(E) dE
        self.eavg: float = (
            smart_integrate(self.energy_profile, self.E_sym, (self.lower, self.upper))
            / self.photden
        ) * u.eV.to(u.erg)

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
    network : Network
        The network the radiation field belongs to; retained on
        ``self.network`` and forwarded to shielding look-ups.
    props : RadiationProps
        Radiation configuration object supplying the band edges, spectral
        index, mode, speed of light, and background field (see
        :class:`~jaff.physics.RadiationProps`).

    Attributes
    ----------
    network : Network
        The network passed to the constructor.
    bands : list of (int, float, or sympy.Basic)
        Band-edge list from ``props.bands`` (``"inf"`` already replaced by
        ``sympy.oo`` by :class:`RadiationProps`).
    mode : str
        ``"nph"`` -- field tracked as photon number density (cm⁻³) -- or
        ``"u"`` -- field tracked as energy density (erg cm⁻³).  Controls the
        name of the symbolic density variable (``"radeden"`` vs. ``"photden"``)
        and the normalisation of rate coefficients.
    c : float or sympy.Symbol
        Speed of light in cm/s (CGS), used in rate-coefficient expressions as
        ``k = c * den * <σ>``.  A string value in ``props.c`` (e.g. ``"c_hat"``
        for a reduced speed of light) becomes a SymPy symbol.
    background_field : BackgroundField
        Background radiation field built from ``props.background_field``.
    nbands : int
        Number of bands (``len(bands) - 1``).
    den : sympy.IndexedBase
        Symbolic radiation-density variable, shape ``(nbands,)``, named
        ``"radeden"`` in energy-density mode or ``"photden"`` otherwise.
    groups : list of RadiationGroup
        One :class:`RadiationGroup` per band, in ascending energy order.  Each
        carries its own spectral index *α_i* (``grp.profile_idx``, from
        ``props.profile_index``; a scalar applies to every band).  Typical
        values: 1 (flat energy spectrum), 2 (flat photon spectrum).
    E_sym : sympy.Symbol
        The photon-energy symbol ``E`` (eV) used in the symbolic profiles.
    nph_profile : sympy.Piecewise
        Piecewise photon-number spectral profile: ``E^(α_i - 2)`` in band *i*
        (see :meth:`_piecewise_profile`).
    energy_profile : sympy.Piecewise
        Piecewise energy-density spectral profile: ``E^(α_i - 1)`` in band *i*.
    photden_tot : sympy.Expr or float
        Integral of ``nph_profile`` over the full band range (the sum of the
        per-band ``grp.photden``), used to normalise the full-spectrum
        cross section ``reaction.rad_xsecs``.
    """

    def __init__(
        self,
        network: Network,
        props: RadiationProps,
    ):
        """Parse band edges and construct one :class:`RadiationGroup` per band.

        Parameters
        ----------
        network : Network
            The network the radiation field belongs to; stored on
            ``self.network`` and used for shielding look-ups.
        props : RadiationProps
            Radiation configuration supplying ``bands`` (photon-energy band
            edges in eV, ``"inf"`` already replaced by ``sympy.oo``),
            ``profile_index`` (spectral index *α* for ``n(E) ∝ E^(α-2)``; a
            scalar or one value per band),
            ``mode`` (``"nph"`` photon number density or ``"u"`` energy
            density), ``c`` (speed of light in cm/s, or a string converted to a
            symbol), and ``background_field``.
        """
        self.network: Network = network
        self.bands: list[int | float | sp.Basic] = props.bands
        self.mode: str = props.mode
        # Speed of light (cm/s) for k = c * σ * n(E) expressions
        self.c: float | sp.Symbol = (
            sp.symbols(props.c) if isinstance(props.c, str) else props.c
        )
        self.background_field = BackgroundField(props.background_field)
        # Global photoionization database; reactions may override it.
        self.pi_database: str = props.pi_database

        self.nbands: int = len(self.bands) - 1
        # Symbolic radiation density variable: energy density (erg/cm³) or
        # photon number density (cm⁻³), depending on the mode.
        self.den = sp.IndexedBase(
            "radeden" if self.mode == "u" else "photden", shape=(self.nbands,)
        )
        self.groups: list[RadiationGroup] = [
            RadiationGroup(
                lower=lower,  # type: ignore
                upper=self.bands[i + 1],
                index=i,
                sym=self.den[i],
                profile_idx=props.profile_index
                if not isinstance(props.profile_index, list)
                else props.profile_index[i],
            )
            for i, lower in enumerate(self.bands[:-1])
        ]
        self.E_sym: sp.Symbol = sp.Symbol("E")
        self.nph_profile: sp.Expr = self._piecewise_profile("nph_profile")
        self.energy_profile: sp.Expr = self._piecewise_profile("energy_profile")

        # ∫ n(E) dE over the full range = sum of the per-band integrals.
        self.photden_tot = sum(grp.photden for grp in self.groups)

    def _piecewise_profile(self, attr: str) -> sp.Expr:
        """Join the per-band profiles ``grp.<attr>`` into one ``sympy.Piecewise``.

        Parameters
        ----------
        attr : str
            Name of the :class:`RadiationGroup` profile attribute
            (``"nph_profile"`` or ``"energy_profile"``).

        Returns
        -------
        sympy.Expr
            Piecewise expression in ``E`` selecting the band containing ``E``
            (``lower <= E < upper``); the last band is the catch-all.
        """
        pieces = [
            (getattr(grp, attr), self.E_sym < grp.upper) for grp in self.groups[:-1]
        ]
        pieces.append((getattr(self.groups[-1], attr), True))

        return sp.Piecewise(*pieces)

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

        If the reaction has no cross section (no ``xsecs_dict``, or no
        ``photodecay`` array / ``photodecay_expr``) the method returns silently (no-op).

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
        The photon-number spectrum used for averaging is the piecewise power
        law ``n(E) ∝ E^(α_i - 2)`` in band *i* (``α_i = grp.profile_idx``).
        For ``α = 1`` this gives a flat energy spectrum; for ``α = 2`` a flat
        photon spectrum.

        Cross-section integrals (``∫ σ n dE``) use
        :func:`~jaff.common._integrators.smart_integrate`: tabulated ``(E, σ)``
        arrays go through ``arr_integrate`` and Verner expressions through
        ``sym_integrate``.  The remaining analytic integrals over ``n(E)``,
        ``E n(E)`` and ``reaction.dRad`` use the same function, which falls
        back to numerical quadrature when SymPy cannot find a closed form.
        """
        xsec: XsecsProps | None = reaction.xsecs_dict
        if xsec is None:
            return

        # Each reaction carries a single decay-channel cross section, either as
        # tabulated arrays (NORAD / Leiden) or as a SymPy expression (Verner).
        # Photon-number spectrum: n(E) ∝ E^(α_i-2) used for weighing the cross-section
        # where α_i = profile_idx of the band containing E.  The factor E^(α-2) arises from
        # n(E) = u(E)/E and u(E) ∝ E^(α-1).
        pr_expr = xsec.get("photodecay_expr")
        pa_weighted: np.ndarray | sp.Expr | None = None

        if pr_expr is not None:
            E: np.ndarray | sp.Symbol = self.E_sym
            pr_weighted: np.ndarray | sp.Expr = pr_expr * self.nph_profile
        else:
            if xsec["photodecay"] is None:
                return

            assert isinstance(xsec["photon_energy"], np.ndarray)
            E = xsec["photon_energy"]  # photon energy array in eV
            ph_profile = self.get_photden_profile(E)
            pr_weighted = xsec["photodecay"] * ph_profile

            if xsec["_equations"]["pa"]:
                pa_weighted = xsec["photo_absorption"] * ph_profile

        k_tot = sp.Float(0.0)  # Accumulates total rate coefficient over all bands

        # Total cross section integrated over the full spectrum (cm²),
        # stored on the reaction for later reference
        xsec_tot = (
            smart_integrate(pr_weighted, E, (self.bands[0], self.bands[-1]))
            / self.photden_tot
        )
        reaction.rad_xsecs = xsec_tot

        # Reset any back-references from a previous run before repopulating.
        reaction.rad_groups = []

        for grp in self.groups:
            # Photon-number-weighted average cross section in the band:
            # <σ>_i = ∫ σ(E) n(E) dE / ∫ n(E) dE
            pr_xsec_avg = (
                smart_integrate(pr_weighted, E, (grp.lower, grp.upper)) / grp.photden
            )
            rad_xsec_avg = (
                smart_integrate(pa_weighted, E, (grp.lower, grp.upper)) / grp.photden
                if pa_weighted is not None
                else pr_xsec_avg
            )

            # Integral of the user-supplied dRad (radiation energy per photon
            # energy, erg/eV) over the band -> energy per reaction event (erg).
            delta_rad_band = smart_integrate(
                reaction.dRad, self.E_sym, (grp.lower, grp.upper)
            )

            # Symbolic rate coefficient: k_i = c · den[i] · <σ>_i
            # (units: s⁻¹ for photon-density mode, cm³ s⁻¹ for two-body)
            k = self.c * self.den[grp.index] * rad_xsec_avg
            if "shielding" in reaction._metadata:
                if "value" in reaction._metadata["shielding"]:
                    k *= reaction._metadata["shielding"]["value"]
                else:
                    k *= Photochemistry.shielding(reaction, self.network)

            grp.props[reaction] = {
                "k": k,
                "xsec": rad_xsec_avg,
                "xsec_frac": rad_xsec_avg / xsec_tot
                if xsec_tot != 0.0
                else 0.0,  # fraction of total cross section
                "delta_rad": delta_rad_band,
            }
            reaction.rad_groups.append(grp)

            # In energy-density mode, convert from "per eV" to "per photon"
            # by dividing by the band-average energy <E>_i.
            k_tot += (
                k
                * (
                    1.0
                    if not xsec["_equations"]["pa"]
                    else (pr_xsec_avg / rad_xsec_avg if rad_xsec_avg != 0.0 else 0.0)
                )
                / (grp.eavg if self.mode == "u" else 1)
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

    def get_photden_profile(self, ph_energy: np.ndarray) -> np.ndarray:
        """Evaluate the piecewise photon-number spectral profile on an energy grid.

        Each energy is assigned to the band containing it (``lower <= E < upper``)
        and evaluated with that band's spectral index.  Energies below
        ``bands[0]`` use the first band's index and energies above
        ``bands[-1]`` use the last band's, so endpoint interpolation in
        :func:`~jaff.common._integrators.arr_integrate` is not biased by zeros
        outside the band range.

        Parameters
        ----------
        ph_energy : numpy.ndarray
            Photon energies in eV.

        Returns
        -------
        numpy.ndarray
            The photon-number profile ``E^(profile_idx_i - 2)`` evaluated at
            each energy, same shape as ``ph_energy``.
        """
        lowers = np.array([float(grp.lower) for grp in self.groups])
        alphas = np.array([float(grp.profile_idx) for grp in self.groups])
        band_idx = np.clip(
            np.searchsorted(lowers, ph_energy, side="right") - 1, 0, self.nbands - 1
        )

        return ph_energy ** (alphas[band_idx] - 2)

    def get_eden_profile(self, ph_energy: np.ndarray) -> np.ndarray:
        """Evaluate the energy-density spectral profile on an energy grid.

        Parameters
        ----------
        ph_energy : numpy.ndarray
            Photon energies in eV.

        Returns
        -------
        numpy.ndarray
            The piecewise energy-density profile ``E^(profile_idx_i - 1)``
            evaluated at each energy, using each band's own spectral index
            (see :meth:`get_photden_profile`).
        """
        return ph_energy * self.get_photden_profile(ph_energy)
