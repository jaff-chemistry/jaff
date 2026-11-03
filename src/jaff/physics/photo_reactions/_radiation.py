"""
Radiation band groups and frequency-integrated rate coefficients.

This module defines two classes:

- :class:`RadiationGroup` -- a single frequency band (energy interval
  ``[lower, upper]`` in eV) that holds per-reaction rate coefficients and
  cross-section data.
- :class:`Radiation` -- the full collection of bands; responsible for
  computing rate coefficients by integrating tabulated photo cross
  sections over each band using a power-law photon-number spectrum.

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
energy (``eavg``).  The rate counts photons, so ``<σ>_i`` is photon-number
weighted.

The radiation a band *loses* is a different quantity in energy-density mode:
each absorbed photon takes its own energy ``E``, so the energy loss weights
the (absorption) cross section by photon energy::

    du_i/dt ⊃ -c * u_i * <σ>_i^(E),   <σ>_i^(E) = ∫ σ E n dE / ∫ E n dE

In photon-density mode the band loses photons at the reaction rate, so the
photon-number weighting applies to both.

Energy units
------------
Band edges, the photon-energy symbol ``ph_E`` and the cross-section tables are
all in **eV**, and the band-average cross sections are energy-unit-free ratios.
The one quantity carrying a net energy dimension, the band-average photon
energy ``eavg`` (``<E>_i``), is converted from eV to **erg** (via
``astropy.units``, ``u.eV.to(u.erg)``) so the rate coefficients and
radiation-moment ODEs are consistent with the CGS solver; the ``radeden``
field is therefore an energy density in erg/cm³.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, List, Tuple

import numpy as np
import sympy as sp
from astropy import units as u

from ...common._integrators import smart_integrate
from ...core.reaction.reaction_props import GroupedPhotoReactionProps
from ._photochemistry import Photochemistry, XsecProps
from .background_field import BackgroundField

if TYPE_CHECKING:
    from ...core.network import Network
    from ...core.reaction import Reaction
    from ...physics import RadiationProps


class RadiationGroup:
    """
    A single frequency band in the radiation field discretisation.

    Each group spans the photon-energy interval ``[lower, upper]`` (in eV).
    Per-reaction band data is not stored here: each photo reaction holds one
    :class:`~jaff.core.reaction.reaction_props.GroupedPhotoReactionProps` per
    band in ``reaction.photo.rgrouped``.

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
    profile_idx : float
        Spectral index *α* of this band.
    ph_E : sympy.Symbol
        The photon-energy symbol (eV), :attr:`NetworkSymbols.ph_E`.
    nph_profile : sympy.Expr
        Photon-number spectral profile ``E^(profile_idx - 2)`` of this band.
    energy_profile : sympy.Expr
        Energy-density spectral profile ``E * nph_profile`` =
        ``E^(profile_idx - 1)`` of this band.
    photden : sympy.Expr or float
        ``∫ nph_profile dE`` over the band, the normalisation for band
        averages.
    eden : sympy.Expr or float
        ``∫ energy_profile dE`` over the band (eV), the normalisation for
        energy-weighted band averages.
    eavg : float
        Photon-number-weighted average energy of this band, in **erg** (the
        band-edge integral is in eV and converted via
        ``astropy.units``, ``u.eV.to(u.erg)``), so that dividing rate/ODE
        terms by it stays CGS-consistent.  Computed on construction and
        shared across all reactions in the band.
    emissivity : sympy.Expr or None
        The network ``globalEmissivity`` integrated over this band, in the units of
        the band's density equation (erg cm⁻³ s⁻¹ in ``"u"`` mode, photons
        cm⁻³ s⁻¹ in ``"nph"`` mode); see :meth:`Radiation.integrate_over_groups`.
        ``None`` without a ``globalEmissivity`` function.
    absorption_coefficient : sympy.Expr or None
        The network ``globalAbsorptionCoefficient`` averaged over this band
        (cm⁻¹), energy-weighted in ``"u"`` mode and photon-number weighted in
        ``"nph"`` mode; see :meth:`Radiation.average_over_groups`.  ``None``
        without a ``globalAbsorptionCoefficient`` function.
    """

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
        from ...core.network import NetworkSymbols

        self.ph_E: sp.Symbol = NetworkSymbols.ph_E
        self.profile_idx: float = profile_idx
        self.nph_profile: sp.Expr = self.ph_E ** (self.profile_idx - 2)
        self.energy_profile: sp.Expr = self.ph_E * self.nph_profile
        # ∫ n(E) dE over the band — used as normalisation for averages.
        self.photden = smart_integrate(
            self.nph_profile, self.ph_E, (self.lower, self.upper)
        )
        # ∫ E n(E) dE over the band — normalisation for energy-weighted averages.
        self.eden = smart_integrate(
            self.energy_profile, self.ph_E, (self.lower, self.upper)
        )
        # Compute the band-average photon energy once per band (shared
        # across all reactions): <E>_i = ∫ E n(E) dE / ∫ n(E) dE
        self.eavg: float = (self.eden / self.photden) * u.eV.to(u.erg)
        # Set by the network once the species are known.
        self.emissivity: sp.Expr | None = None
        self.absorption_coefficient: sp.Expr | None = None

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
    ph_E : sympy.Symbol
        The photon-energy symbol (eV), :attr:`NetworkSymbols.ph_E`, used in the
        symbolic profiles.
    nph_profile : sympy.Piecewise
        Piecewise photon-number spectral profile: ``E^(α_i - 2)`` in band *i*
        (see :meth:`_piecewise_profile`).
    energy_profile : sympy.Piecewise
        Piecewise energy-density spectral profile: ``E^(α_i - 1)`` in band *i*.
    photden_tot : sympy.Expr or float
        Integral of ``nph_profile`` over the full band range (the sum of the
        per-band ``grp.photden``), used to normalise the full-spectrum
        cross section ``reaction.photo.xsec.total``.
    emissivity : sympy.Expr or None
        The network ``.jfunc`` ``globalEmissivity``: the *volume* emissivity,
        i.e. radiation energy emitted per unit volume, time and photon energy
        (erg cm⁻³ s⁻¹ eV⁻¹), integrated over solid angle, as an expression in
        ``ph_E``.  ``None`` without one.  Its
        per-group integrals are ``grp.emissivity``; their sum is
        :attr:`total_emissivity`.
    absorption_coefficient : sympy.Expr or None
        The network ``.jfunc`` ``globalAbsorptionCoefficient``: the fraction of
        the radiation absorbed per unit path length (cm⁻¹) as a function of
        ``ph_E`` (eV).  ``None`` without one.  Its per-group averages are
        ``grp.absorption_coefficient``.
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
        self.override_custom_rate: bool = props.override_custom_rate
        self.override_custom_cross_section: bool = props.override_custom_cross_section
        # Whether the global absorption coefficient also attenuates the flux.
        self.absorption_attenuates_flux: bool = props.absorption_attenuates_flux

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
        self.ph_E: sp.Symbol = self.groups[0].ph_E
        self.nph_profile: sp.Expr = self._piecewise_profile("nph_profile")
        self.energy_profile: sp.Expr = self._piecewise_profile("energy_profile")

        # ∫ n(E) dE over the full range = sum of the per-band integrals.
        self.photden_tot = sum(grp.photden for grp in self.groups)

        # Set by the network once the species are known.
        self.emissivity: sp.Expr | None = None
        self.absorption_coefficient: sp.Expr | None = None

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
            Piecewise expression in ``ph_E`` selecting the band containing ``ph_E``
            (``lower <= E < upper``); the last band is the catch-all.
        """
        pieces = [(getattr(grp, attr), self.ph_E < grp.upper) for grp in self.groups[:-1]]
        pieces.append((getattr(self.groups[-1], attr), True))

        return sp.Piecewise(*pieces)

    def get_reaction_rate_props(
        self, reaction: Reaction
    ) -> Tuple[sp.Expr, float, List[GroupedPhotoReactionProps]]:
        """
        Compute and store symbolic band-averaged rate coefficients for a reaction.

        Reads the tabulated photo cross section σ(E) from
        ``reaction.photo.xsec``, then for each frequency band:

        1. Computes the photon-number-weighted band-average cross section
           ``<σ>_i = ∫ σ n dE / ∫ n dE``.
        2. Computes the band-average photon energy
           ``<E>_i = ∫ E n dE / ∫ n dE`` (stored once per band).
        3. Assembles the symbolic rate coefficient
           ``k_i = c * den[i] * <σ>_i`` (photon-density mode) or
           ``k_i = c * den[i] * <σ>_i / <E>_i`` (energy-density mode).
        4. Assembles the band loss coefficient ``c * den[i] * <σ>_i`` that
           :func:`~jaff.physics._equations.get_sradodes` removes from the band,
           with the absorption cross section energy-weighted in
           energy-density mode (each absorbed photon takes its own energy).
        5. Collects the loss coefficient (as ``k``), the band-average cross
           section and the cross-section fraction into one
           :class:`GroupedPhotoReactionProps` per band.

        After iterating over all bands the total symbolic rate coefficient
        (sum over bands, in units of s⁻¹ or cm³ s⁻¹ depending on reaction
        type) is returned as ``k_tot``.

        If the reaction has no cross section (no ``photo.xsec``, or no
        ``photodecay`` array / ``photodecay_expr``) a ``ValueError`` is raised.

        Parameters
        ----------
        reaction : Reaction
            The photochemical reaction to process.

        Returns
        -------
        tuple[sympy.Expr, float, list[GroupedPhotoReactionProps]]
            ``(k_tot, xsec_tot, grouped)``: the total rate coefficient summed
            over bands, the cross section averaged over the whole spectrum
            (cm²), and one entry per band in ``self.groups`` order.  The
            caller stores them on ``reaction.rate_coefficient``,
            ``reaction.photo.xsec.total`` and ``reaction.photo.rgrouped``.

        Notes
        -----
        The photon-number spectrum used for averaging is the piecewise power
        law ``n(E) ∝ E^(α_i - 2)`` in band *i* (``α_i = grp.profile_idx``).
        For ``α = 1`` this gives a flat energy spectrum; for ``α = 2`` a flat
        photon spectrum.

        Cross-section integrals (``∫ σ n dE``) use
        :func:`~jaff.common._integrators.smart_integrate`: tabulated ``(E, σ)``
        arrays go through ``arr_integrate`` and Verner expressions through
        ``sym_integrate``.  The remaining analytic integrals over ``n(E)`` and
        ``E n(E)`` use the same function, which falls back to numerical
        quadrature when SymPy cannot find a closed form.
        """
        photo = reaction.photo
        if photo is None or photo.xsec is None:
            raise ValueError(
                f"Photo reaction {reaction} has no cross section (photo.xsec is None)"
            )

        xsec: XsecProps = photo.xsec

        # Each reaction carries a single decay-channel cross section, either as
        # tabulated arrays (NORAD / Leiden) or as a SymPy expression (Verner);
        # only tabulated sources can carry a separate absorption cross section.
        absorption: np.ndarray | None = None
        if xsec.photodecay_expr is not None:
            E: np.ndarray | sp.Symbol = self.ph_E
            decay: np.ndarray | sp.Expr = xsec.photodecay_expr
        else:
            if xsec.photodecay is None:
                raise ValueError(
                    f"Photodecay cross-section cannot be None for reaction: {reaction}"
                )

            assert isinstance(xsec.photon_energy, np.ndarray)
            E = xsec.photon_energy  # photon energy array in eV
            decay = xsec.photodecay

            if xsec.process.absorption:
                absorption = xsec.photo_absorption

        k_tot: sp.Expr = sp.Float(
            0.0
        )  # Accumulates total rate coefficient over all bands

        # Total cross section averaged over the full spectrum (cm²),
        # stored on the reaction for later reference
        profile = self.get_photden_profile(E) if isinstance(E, np.ndarray) else None
        xsec_tot: float = self.get_total_average_cross_section(decay, E, profile)

        # Shielding factor, the same for every band: a user-fixed "value", or
        # the shielding model evaluated once (it may write per-reaction tables).
        shielding: sp.Expr | float | None = None
        if photo.shielding is not None:
            shielding = (
                photo.shielding.props["value"]
                if "value" in photo.shielding.props
                else Photochemistry.shielding(reaction, self.network)
            )

        # Reset any back-references from a previous run before repopulating.
        grouped: List[GroupedPhotoReactionProps] = []

        for grp in self.groups:
            # Photon-number-weighted average cross section in the band:
            # <σ>_i = ∫ σ(E) n(E) dE / ∫ n(E) dE
            pr_xsec_avg: float = self.get_group_average_cross_section(
                decay, E, grp, profile
            )
            rad_xsec_avg: float = (
                self.get_group_average_cross_section(absorption, E, grp, profile)
                if absorption is not None
                else pr_xsec_avg
            )
            loss_xsec_avg = (
                self.get_group_average_cross_section(
                    decay if absorption is None else absorption,
                    E,
                    grp,
                    profile,
                    energy_weighted=True,
                )
                if self.mode == "u"
                else rad_xsec_avg
            )

            # Photon-counting rate coefficient: k_rate = c · den[i] · <σ>_i
            # (units: s⁻¹ for photon-density mode, cm³ s⁻¹ for two-body)
            k_rate: sp.Expr = self.c * self.den[grp.index] * rad_xsec_avg
            # Band loss coefficient: what get_sradodes removes from den[i].
            k: sp.Expr = self.c * self.den[grp.index] * loss_xsec_avg
            if shielding is not None:
                k_rate *= shielding
                k *= shielding

            rgroup: GroupedPhotoReactionProps = GroupedPhotoReactionProps(
                k=k,
                xsec=rad_xsec_avg,
                xsec_frac=rad_xsec_avg / xsec_tot
                if xsec_tot != 0.0
                else 0.0,  # fraction of total cross section
            )
            grouped.append(rgroup)

            # In energy-density mode, convert from "per eV" to "per photon"
            # by dividing by the band-average energy <E>_i.
            k_tot += (
                k_rate
                * (
                    1.0
                    if not xsec.process.absorption
                    else (pr_xsec_avg / rad_xsec_avg if rad_xsec_avg != 0.0 else 0.0)
                )
                / (grp.eavg if self.mode == "u" else 1)
            )

        return k_tot, xsec_tot, grouped

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

    def get_total_average_cross_section(
        self,
        xsec: np.ndarray | sp.Expr,
        ph_E: np.ndarray | sp.Symbol,
        ph_profile: np.ndarray | None = None,
    ) -> float:
        """Photon-number-weighted mean of *xsec* over the whole spectrum.

        ``<σ> = ∫ σ(E) n(E) dE / ∫ n(E) dE`` from ``bands[0]`` to ``bands[-1]``,
        with the piecewise power-law ``n(E)`` of the bands.

        Parameters
        ----------
        xsec : numpy.ndarray or sympy.Expr
            Cross section (cm²): tabulated on the *ph_E* grid, or an
            expression in the symbol *ph_E*.
        ph_E : numpy.ndarray or sympy.Symbol
            Photon-energy grid (eV) of a tabulated *xsec*, or the photon-energy
            symbol of a symbolic one.
        ph_profile : numpy.ndarray or None, optional
            :meth:`get_photden_profile` already evaluated on a tabulated
            *ph_E* grid, to reuse across calls; ``None`` (default) evaluates it.

        Returns
        -------
        float
            The spectrum-averaged cross section (cm²).
        """
        profile: np.ndarray | sp.Expr = (
            ph_profile
            if ph_profile is not None
            else self.get_photden_profile(ph_E)
            if isinstance(ph_E, np.ndarray)
            else self.nph_profile
        )

        return (
            smart_integrate(
                xsec * profile,
                ph_E,
                (self.bands[0], self.bands[-1]),
            )
            / self.photden_tot
        )

    def get_group_average_cross_section(
        self,
        xsec: np.ndarray | sp.Expr,
        ph_E: np.ndarray | sp.Symbol,
        group: RadiationGroup,
        ph_profile: np.ndarray | None = None,
        energy_weighted: bool = False,
    ) -> float:
        """Photon-number- or energy-weighted mean of *xsec* over one band.

        ``<σ>_i = ∫ σ(E) n(E) dE / ∫ n(E) dE`` over ``[group.lower, group.upper]``,
        or, with *energy_weighted*, ``∫ σ(E) E n(E) dE / ∫ E n(E) dE``.  The
        photon-number weighting gives reaction rates (they count photons); the
        energy weighting gives the energy a band loses in ``"u"`` mode (each
        absorbed photon takes its own energy ``E``).

        Parameters
        ----------
        xsec : numpy.ndarray or sympy.Expr
            Cross section (cm²): tabulated on the *ph_E* grid, or an
            expression in the symbol *ph_E*.
        ph_E : numpy.ndarray or sympy.Symbol
            Photon-energy grid (eV) of a tabulated *xsec*, or the photon-energy
            symbol of a symbolic one.
        group : RadiationGroup
            The band to average over.
        ph_profile : numpy.ndarray or None, optional
            :meth:`get_photden_profile` already evaluated on a tabulated
            *ph_E* grid, to reuse across calls; ``None`` (default) evaluates it.
        energy_weighted : bool, optional
            Weight by the energy profile ``E n(E)`` instead of ``n(E)``
            (default ``False``).

        Returns
        -------
        float
            The band-averaged cross section (cm²).
        """
        profile: np.ndarray | sp.Expr = (
            ph_profile
            if ph_profile is not None
            else self.get_photden_profile(ph_E)
            if isinstance(ph_E, np.ndarray)
            else group.nph_profile
        )
        norm = group.photden
        if energy_weighted:
            profile = profile * ph_E
            norm = group.eden

        return smart_integrate(xsec * profile, ph_E, (group.lower, group.upper)) / norm

    def integrate_over_groups(self, spectrum: sp.Expr) -> List[sp.Expr]:
        """Integrate an emitted spectrum over each group, in density-equation units.

        *spectrum* is emitted energy per unit photon energy (erg/eV, times
        whatever it is "per": a reaction event for a ``spectralYield``, a
        volume and time for the ``globalEmissivity``).  Each group's integral is in
        the units of that group's density variable:

        - ``"u"`` mode: ``∫ S(E) dE`` over the group, the energy added (erg);
        - ``"nph"`` mode: ``∫ S(E) / E dE`` over the group divided by the
          eV-to-erg factor (one photon carries ``E`` eV), the photons added.

        Parameters
        ----------
        spectrum : sympy.Expr
            The spectrum as an expression in :attr:`ph_E` (eV).  Other symbols
            (e.g. ``tgas``) stay symbolic (see
            :func:`~jaff.common._integrators.smart_integrate`).

        Returns
        -------
        list[sympy.Expr]
            One entry per group in ``self.groups`` order.

        Raises
        ------
        ValueError
            If a spectrum with symbols other than ``ph_E`` has no closed-form
            integral over a group.
        """
        if self.mode == "u":
            return [
                sp.sympify(smart_integrate(spectrum, self.ph_E, (grp.lower, grp.upper)))
                for grp in self.groups
            ]

        return [
            sp.sympify(
                smart_integrate(spectrum / self.ph_E, self.ph_E, (grp.lower, grp.upper))
            )
            / u.eV.to(u.erg)
            for grp in self.groups
        ]

    def average_over_groups(self, spectrum: sp.Expr) -> List[sp.Expr]:
        """Average a per-photon-energy coefficient over each group.

        *spectrum* is a quantity defined at each photon energy, not per eV (e.g.
        an absorption coefficient in cm⁻¹), so it is averaged rather than
        integrated, weighted by the group's assumed spectrum: by the energy
        profile ``E n(E)`` in ``"u"`` mode (the energy density loses energy) and
        by the photon profile ``n(E)`` in ``"nph"`` mode (the photon density
        loses photons).

        Parameters
        ----------
        spectrum : sympy.Expr
            The coefficient as an expression in :attr:`ph_E` (eV).  Other
            symbols (e.g. ``tgas``) stay symbolic (see
            :func:`~jaff.common._integrators.smart_integrate`).

        Returns
        -------
        list[sympy.Expr]
            One weighted mean per group in ``self.groups`` order.

        Raises
        ------
        ValueError
            If a spectrum with symbols other than ``ph_E`` has no closed-form
            integral over a group.
        """
        return [
            sp.sympify(
                self.get_group_average_cross_section(
                    spectrum, self.ph_E, grp, energy_weighted=self.mode == "u"
                )
            )
            for grp in self.groups
        ]

    def group_property(self, name: str) -> List[Any]:
        """Value of the :class:`RadiationGroup` attribute *name* for every group.

        Parameters
        ----------
        name : str
            Attribute of :class:`RadiationGroup`, e.g. ``"profile_idx"``,
            ``"lower"`` or ``"upper"``.

        Returns
        -------
        list
            One value per group, in ``self.groups`` order.  The ``upper`` of the
            last group is ``sympy.oo`` when the last band edge is ``"inf"``.

        Raises
        ------
        AttributeError
            If :class:`RadiationGroup` has no attribute *name*.
        """
        return [getattr(grp, name) for grp in self.groups]

    @property
    def total_emissivity(self) -> sp.Expr | None:
        """Sum of the group emissivities, or ``None`` without a ``globalEmissivity``.

        The groups are contiguous, so this is the emissivity integrated over the
        whole band range, in the units of the density equations (erg cm⁻³ s⁻¹
        in ``"u"`` mode, photons cm⁻³ s⁻¹ in ``"nph"`` mode).
        """
        if self.emissivity is None:
            return None

        return sum((grp.emissivity for grp in self.groups), sp.S.Zero)  # type: ignore[arg-type]
