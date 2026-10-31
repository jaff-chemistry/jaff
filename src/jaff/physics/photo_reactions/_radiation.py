# ABOUTME: Radiation and RadiationGroup: the radiation bands and the band-averaged
# ABOUTME: symbolic photo-reaction rate coefficients
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
The photon number spectrum is assumed to follow a power law in photon energy:

```text
n(E) ∝ E^(α - 2)
```

where ``α`` is the spectral index ``profile_idx`` (from
``RadiationConfig.profile_index``).  ``α`` may be a single global value or one
value per band, in which case the spectrum is a piecewise power law (a
histogram, intentionally discontinuous at band edges): each band ``i`` uses its
own ``α_i`` and is normalised independently.  The energy-integrated version
(energy density per unit energy interval) is therefore:

```text
u(E) = E * n(E) ∝ E^(α - 1)
```

This form is used when computing band-average cross sections and average
photon energies.

Rate coefficient derivation
----------------------------
For a reaction with tabulated cross section σ(E) the *photon-flux-weighted*
average cross section in band *i* is:

```text
<σ>_i = ∫_{E_lo}^{E_hi} σ(E) n(E) dE  /  ∫_{E_lo}^{E_hi} n(E) dE
```

The symbolic rate coefficient stored in the radiation density variable
``den[i]`` (either energy density *u_i* in erg/cm³ or photon density *n_i*
in cm⁻³) is:

```text
k_i = c * den[i] * <σ>_i          (photon density mode)
k_i = c * den[i] * <σ>_i / <E>_i  (energy density mode)
```

where *c* is the speed of light and ``<E>_i`` is the band-average photon
energy (``eavg``).  The rate counts photons, so ``<σ>_i`` is photon-number
weighted.

The radiation a band *loses* is a different quantity in energy-density mode:
each absorbed photon takes its own energy ``E``, so the energy loss weights
the (absorption) cross section by photon energy:

```text
du_i/dt ⊃ -c * u_i * <σ>_i^(E),   <σ>_i^(E) = ∫ σ E n dE / ∫ E n dE
```

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

from typing import TYPE_CHECKING, ClassVar, List

import numpy as np
import sympy as sp
from astropy import units as u

from ...common._integrators import smart_integrate
from ...io import JaffLogger
from ...io._logger import jaff_progress
from .background_field import BackgroundField

if TYPE_CHECKING:
    from ...core.network import Network
    from .config import RadiationConfig


class RadiationGroup:
    """
    A single frequency band in the radiation field discretisation.

    Each group spans the photon-energy interval ``[lower, upper]`` (in eV) and
    holds that band's spectral profile and its integrals, the normalisations of
    the band averages.  :class:`Radiation` builds one group per band
    (``net.radiation.groups``); you normally only read them.

    Per-reaction band data is not stored here: each photo reaction holds one
    :class:`~jaff.physics.PhotoBand` per band in ``reaction.photo.bands``, each
    holding its group; ``reaction.photo.band_xsecs`` shows them as a table.

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
        without a ``globalAbsorptionCoefficient`` function.  It is removed from
        the band's density equation as ``c * <α>_i * den[i]``.

    Examples
    --------
    >>> from jaff import Network
    >>> net = Network(
    ...     "h_photoionization",
    ...     radiation={"bands": [1, 13.6, 100, "inf"]},
    ... )
    >>> grp = net.radiation.groups[1]
    >>> grp, grp.sym, grp.lower, grp.upper
    (Rad_group(1, band=(13.6, 100.0)), photden[1], 13.6, 100.0)
    >>> grp.nph_profile, grp.energy_profile
    (ph_E**(-2.0), ph_E**(-1.0))
    >>> grp.photden  # ∫ E^-2 dE from 13.6 to 100 eV
    0.06352941176470589
    >>> net.radiation.groups[2].upper  # open band (100 eV, oo)
    oo
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

        Raises
        ------
        ValueError
            If *lower* is not an ``int`` or ``float``.
        """
        self.index: int = index
        self.sym: sp.Basic = sym
        self.lower: float | int = lower
        self.upper: float | int | sp.Basic = upper

        if not isinstance(self.lower, (int, float)):
            raise ValueError(
                f"Radiation group lower bound must be a float or int. Found {self.lower}"
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
        return f"Rad_group({self.index}, band={(self.lower, self.upper)})"

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

    The radiation field is split into contiguous photon-energy *bands*.  On
    construction the band-edge list is parsed (replacing the string ``"inf"``
    with ``sympy.oo``) and a :class:`RadiationGroup` object is created for each
    consecutive pair of edges.

    :class:`~jaff.core.Network` builds it from its radiation config, given as
    the ``radiation`` argument (a dict or a
    :class:`~jaff.physics.photo_reactions.config.RadiationConfig`) or as the
    ``[network.radiation]`` table of a config file, and stores it as
    ``net.radiation`` (``None`` when radiation is off).  Every photo reaction
    then assembles its rate coefficient and per-band data from the band
    integrals offered here (see
    :meth:`~jaff.core.reaction.PhotoReaction.band_rates`).

    Parameters
    ----------
    network : Network
        The network the radiation field belongs to; retained on
        ``self.network`` and forwarded to shielding look-ups.
    config : RadiationConfig
        The ``[network.radiation]`` config: band edges, spectral index, mode,
        speed of light, background field and absorption flag.

    Attributes
    ----------
    network : Network
        The network passed to the constructor.
    config : RadiationConfig
        The ``[network.radiation]`` settings (read-only).  ``config.mode`` is
        ``"nph"`` (field tracked as photon number density, cm⁻³) or ``"u"``
        (energy density, erg cm⁻³): it names the density variable
        (``"photden"`` vs. ``"radeden"``) and sets the normalisation of the
        rate coefficients.  ``config.absorption_attenuates_flux`` says whether
        the ``globalAbsorptionCoefficient`` also attenuates the flux equations
        (``-c * <α>_i * rflux[i]``); see :meth:`odes`.
    bands : list of (int, float, or sympy.Basic)
        Band-edge list from ``config.bands``, ``"inf"`` replaced by
        ``sympy.oo``.
    c : float or sympy.Symbol
        Speed of light in cm/s (CGS), used in rate-coefficient expressions as
        ``k = c * den * <σ>``.  A string value in ``config.c`` (e.g. ``"c_hat"``
        for a reduced speed of light) becomes a SymPy symbol.
    background_field : BackgroundField
        Background radiation field built from ``config.background_field``.
    nbands : int
        Number of bands (``len(bands) - 1``).
    den : sympy.IndexedBase
        Symbolic radiation-density variable, shape ``(nbands,)``, named
        ``"radeden"`` in energy-density mode or ``"photden"`` otherwise.
    flux : sympy.IndexedBase
        Symbolic radiation-flux variable ``rflux``, shape ``(nbands,)``: the
        flux moment of each band, in the units of :attr:`den` times a speed.
    groups : list of RadiationGroup
        One :class:`RadiationGroup` per band, in ascending energy order.  Each
        carries its own spectral index *α_i* (``grp.profile_idx``, from
        ``config.profile_index``; a scalar applies to every band).  Typical
        values: 1 (flat energy spectrum), 2 (flat photon spectrum).
    ph_E : sympy.Symbol
        The photon-energy symbol (eV), :attr:`NetworkSymbols.ph_E`, used in the
        symbolic profiles.
    nph_profile : sympy.Piecewise
        Piecewise photon-number spectral profile: ``E^(α_i - 2)`` in band *i*.
        A point on an interior edge belongs to the band above it; energies
        below ``bands[0]`` use the first band's profile and energies above
        ``bands[-1]`` the last band's.  When every band has the same index the
        pieces collapse into one power law.
    energy_profile : sympy.Piecewise
        Piecewise energy-density spectral profile: ``E^(α_i - 1)`` in band *i*
        (same edge rules).
    photden_tot : sympy.Expr or float
        Integral of ``nph_profile`` over the full band range (the sum of the
        per-band ``grp.photden``), used to normalise the full-spectrum
        cross section ``reaction.photo.xsec_total``.
    emissivity : sympy.Expr or None
        The network ``.jfunc`` ``globalEmissivity``: the *volume* emissivity,
        i.e. radiation energy emitted per unit volume, time and photon energy
        (erg cm⁻³ s⁻¹ eV⁻¹), integrated over solid angle, as an expression in
        ``ph_E``.  ``None`` without one.  Its
        per-group integrals are ``grp.emissivity``; their sum is
        :attr:`total_emissivity`.  Set by :meth:`compute`.
    absorption_coefficient : sympy.Expr or None
        The network ``.jfunc`` ``globalAbsorptionCoefficient``: the fraction of
        the radiation absorbed per unit path length (cm⁻¹) as a function of
        ``ph_E`` (eV).  ``None`` without one.  Its per-group averages are
        ``grp.absorption_coefficient``.  Set by :meth:`compute`.

    Notes
    -----
    Within band *i* the photon-number spectrum is a power law with that band's
    own spectral index ``α_i`` (``config.profile_index``: a scalar for every
    band, or one value per band):

    ```
    n(E) ∝ E^(α_i - 2),    u(E) = E n(E) ∝ E^(α_i - 1)
    ```

    Each band is normalised on its own, so the spectrum is a histogram,
    discontinuous at the band edges by design.  Common values are ``α = 1``
    (flat energy spectrum) and ``α = 2`` (flat photon spectrum).

    The rate coefficient of a photo reaction is the sum over bands of

    ```
    k_i = c * den[i] * <σ>_i            ("nph" mode)
    k_i = c * den[i] * <σ>_i / <E>_i    ("u" mode)
    ```

    with ``<σ>_i`` the photon-number-weighted band-average cross section and
    ``<E>_i`` the band-average photon energy (``grp.eavg``, erg), i.e.
    ``c * photon_density(group) * <σ>_i`` (:meth:`photon_density`).  Each photo
    reaction assembles it (:meth:`~jaff.core.reaction.PhotoReaction.band_rates`)
    and stores it as ``reaction.photo.band_rate``, ``reaction.photo.xsec_total``
    and ``reaction.photo.bands``.  :meth:`odes` builds the
    radiation-moment equations from them, laid out by :meth:`ordered_index`.

    Examples
    --------
    >>> import numpy as np
    >>> from jaff import Network
    >>> net = Network(
    ...     "h_photoionization",
    ...     radiation={
    ...         "bands": [1, 13.6, 100, "inf"],
    ...         "profile_index": [0, 1, 0],
    ...     },
    ... )
    >>> rad = net.radiation
    >>> rad.bands, rad.nbands, rad.config.mode, rad.den
    ([1.0, 13.6, 100.0, oo], 3, 'nph', photden)
    >>> rad.groups
    [Rad_group(0, band=(1.0, 13.6)), Rad_group(1, band=(13.6, 100.0)), Rad_group(2, band=(100.0, oo))]
    >>> [grp.profile_idx for grp in rad.groups]
    [0.0, 1.0, 0.0]
    >>> rad.nph_profile
    Piecewise((ph_E**(-2.0), ph_E < 13.6), (ph_E**(-1.0), ph_E < 100.0), (ph_E**(-2.0), True))
    >>> float(rad.photden_tot)  # (1 - 1/13.6) + ln(100/13.6) + 1/100
    2.931570981481379
    """

    ORDERS: ClassVar[tuple[int, ...]] = (0, 1, 2, 3)
    """The layouts of the radiation equations (see :meth:`ordered_index`)."""

    def __init__(
        self,
        network: Network,
        config: RadiationConfig,
    ):
        """Parse band edges and construct one :class:`RadiationGroup` per band.

        Parameters
        ----------
        network : Network
            The network the radiation field belongs to; stored on
            ``self.network`` and used for shielding look-ups.
        config : RadiationConfig
            The ``[network.radiation]`` config: ``bands`` (photon-energy band
            edges in eV, ``"inf"`` allowed last), ``profile_index`` (spectral
            index *α* for ``n(E) ∝ E^(α-2)``; a scalar or one value per band),
            ``mode`` (``"nph"`` photon number density or ``"u"`` energy
            density), ``c`` (speed of light in cm/s, or a string converted to a
            symbol), ``background_field`` and ``absorption_attenuates_flux``.

        Raises
        ------
        RuntimeError
            If the band-averaged photon energy diverges at an open band edge
            (energy-density mode).
        """
        self.logger = JaffLogger().get_logger()
        self.network: Network = network
        self.config: RadiationConfig = config
        profile_indices = self._profile_indices(config)
        self.bands: list[int | float | sp.Basic] = [
            sp.oo if edge == "inf" else edge for edge in config.bands
        ]
        self._check_convergence(profile_indices)
        # Speed of light (cm/s) for k = c * σ * n(E) expressions
        self.c: float | sp.Symbol = (
            sp.symbols(config.c) if isinstance(config.c, str) else config.c
        )
        self.background_field = BackgroundField(config.background_field)

        self.nbands: int = len(self.bands) - 1
        # Symbolic radiation density variable: energy density (erg/cm³) or
        # photon number density (cm⁻³), depending on the mode.
        self.den = sp.IndexedBase(
            "radeden" if self.config.mode == "u" else "photden", shape=(self.nbands,)
        )
        # Symbolic radiation flux variable (the first moment), one per band.
        self.flux = sp.IndexedBase("rflux", shape=(self.nbands,))
        self.groups: list[RadiationGroup] = [
            RadiationGroup(
                lower=lower,  # type: ignore
                upper=self.bands[i + 1],
                index=i,
                sym=self.den[i],
                profile_idx=profile_indices[i],
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

    @staticmethod
    def _profile_indices(config: RadiationConfig) -> list[float]:
        """Return the spectral index *α_i* of each band (a scalar applies to all)."""
        nbands = len(config.bands) - 1
        if isinstance(config.profile_index, list):
            return [float(alpha) for alpha in config.profile_index]

        return [float(config.profile_index)] * nbands

    def _check_convergence(self, alphas: list[float]) -> None:
        """Check the band-averaged photon energy converges at the open band edges.

        Only the first band reaches ``E → 0`` and only the last reaches
        ``E → ∞``, so each edge is checked against its own band's index.

        Raises
        ------
        RuntimeError
            In energy-density mode, if ``<E>`` diverges at a zero lower edge or
            an infinite upper edge.
        """
        bands = self.bands
        alpha_lo, alpha_hi = alphas[0], alphas[-1]
        starts_at_zero = isinstance(bands[0], (float, int)) and float(bands[0]) == 0.0
        ends_at_inf = bands[-1] == sp.oo

        if self.config.mode == "u":
            # The average-energy numerator ∫ E·n(E) dE is the integral of the
            # energy-density spectrum u(E) = E·n(E) ∝ E^(α-1), i.e. a power law
            # with index pl_index = (α-2) + 1 = α - 1.
            # pl_index + 1 = α: α < 0 diverges at E → 0, α > 0 diverges at
            # E → ∞, and α == 0 (pl_index == -1) log-diverges at both ends.
            if starts_at_zero and alpha_lo <= 0.0:
                raise RuntimeError(
                    f"The integral for average energy will diverge since the radiation band starts from bands[0]: {bands[0]}\n"
                    "Please try a non-zero value"
                )
            if ends_at_inf and alpha_hi >= 0.0:
                raise RuntimeError(
                    f'The integral for average energy will diverge since the radiation band ends at bands[{len(bands) - 1}]: "inf"\n'
                    "Please try a non-infinite value or change the profile_index"
                )

        if (
            alpha_lo <= 1.0
            and isinstance(bands[0], (float, int))
            and float(bands[0]) < 1.0
        ):
            self.logger.warning(
                f"Radiation band starts at bands[0]={bands[0]} eV with "
                f"profile_index={alpha_lo} (first band): the photon-number "
                "normalisation integral ∫E^(α-2)dE is lower-edge divergent "
                "(exponent ≤ -1) and near E→0 becomes ill-conditioned, which "
                "can yield a negative/garbage photon density. Use a non-zero "
                "bands[0] ≳ 1 eV or a larger profile_index."
            )

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

    def photon_density(self, group: RadiationGroup) -> sp.Expr:
        """Photon number density of *group* (cm⁻³), from its :attr:`den` symbol.

        ``den[i]`` itself in ``"nph"`` mode; in ``"u"`` mode the energy density
        divided by the band-average photon energy, ``den[i] / <E>_i``.

        Parameters
        ----------
        group : RadiationGroup
            The band.

        Returns
        -------
        sympy.Expr
            The photon density.
        """
        density = self.den[group.index]
        if self.config.mode == "u":
            return density / group.eavg

        return density

    def energy_density(self, group: RadiationGroup) -> sp.Expr:
        """Radiation energy density of *group* (erg cm⁻³), from its :attr:`den` symbol.

        ``den[i]`` itself in ``"u"`` mode; in ``"nph"`` mode the photon density
        times the band-average photon energy, ``den[i] · <E>_i``.

        Parameters
        ----------
        group : RadiationGroup
            The band.

        Returns
        -------
        sympy.Expr
            The energy density.
        """
        density = self.den[group.index]
        if self.config.mode == "nph":
            return density * group.eavg

        return density

    def loss_cross_section(
        self,
        xsec: np.ndarray | sp.Expr,
        ph_E: np.ndarray | sp.Symbol,
        group: RadiationGroup,
        ph_profile: np.ndarray | None,
        photon_average: float | sp.Expr,
    ) -> float | sp.Expr:
        """Band-average cross section that sets how fast an absorber drains *group*.

        Each absorbed photon removes one photon (``"nph"`` mode: the
        photon-number-weighted *photon_average*) or its own energy (``"u"``
        mode: the energy-weighted average of *xsec*, see
        :meth:`get_group_average_cross_section`).

        Parameters
        ----------
        xsec : numpy.ndarray or sympy.Expr
            The absorbing cross section (cm²), as for
            :meth:`get_group_average_cross_section`.
        ph_E : numpy.ndarray or sympy.Symbol
            Its photon-energy grid (eV) or symbol.
        group : RadiationGroup
            The band.
        ph_profile : numpy.ndarray or None
            :meth:`get_photden_profile` on a tabulated *ph_E*, or ``None``.
        photon_average : float or sympy.Expr
            The photon-number-weighted band average of *xsec*, already computed.

        Returns
        -------
        float or sympy.Expr
            The loss cross section (cm²).
        """
        if self.config.mode != "u":
            return photon_average

        return self.get_group_average_cross_section(
            xsec, ph_E, group, ph_profile, energy_weighted=True
        )

    def odes(self, order: int = 0) -> List[sp.Expr]:
        """Symbolic right-hand sides of the radiation-moment equations of every band.

        Each band *i* has two moments, the density ``den[i]`` (``radeden`` in
        erg cm⁻³, or ``photden`` in cm⁻³) and the flux ``rflux[i]``
        (:attr:`flux`).  Its density equation holds:

        - the losses to photo reactions (each band's loss coefficient
          ``bands[i].k`` times the reaction's reactant densities); photo
          reactions whose rate comes from a ``.jfunc`` ``chemRateCoeff<N>`` are not
          coupled to the bands and are skipped;
        - the absorption by the medium ``-c · <α>_i · den[i]`` (a
          ``globalAbsorptionCoefficient``);
        - the radiation emitted by reactions with a spectral yield
          (``reaction.rad.group_yields`` times the reaction rate) and the band's
          ``globalEmissivity``.

        Its flux equation is the density losses with ``den[i] → rflux[i]``;
        emission is isotropic and feeds the density only, and the absorption
        enters it only with ``absorption_attenuates_flux``.  With dust, the dust
        attenuation (:meth:`~jaff.physics.Dust.attenuation`) is subtracted from
        both.  Transport terms (the flux divergence, the radiation pressure) are
        not included.

        Parameters
        ----------
        order : {0, 1, 2, 3}, optional
            Layout of the result (see :meth:`ordered_index`):

            - ``0`` (default): ``[den_0, flux_0, den_1, flux_1, ...]``;
            - ``1``: ``[flux_0, den_0, flux_1, den_1, ...]``;
            - ``2``: ``[den_0, den_1, ..., flux_0, flux_1, ...]``;
            - ``3``: ``[flux_0, flux_1, ..., den_0, den_1, ...]``.

        Returns
        -------
        list of sympy.Expr
            ``2 * nbands`` expressions, laid out as *order* says.

        Raises
        ------
        ValueError
            If *order* is not one of ``0``, ``1``, ``2``, ``3``.
        RuntimeError
            If a photo reaction without a custom rate has no per-band data.

        Examples
        --------
        >>> from jaff import Network
        >>> net = Network(
        ...     "h_photoionization",
        ...     radiation={"bands": [1, 13.6, 100, "inf"]},
        ... )
        >>> net.radiation.odes()  # [den_0, flux_0, den_1, flux_1, den_2, flux_2]
        [0, 0, -5.85154082319921e-8*nden[0]*photden[1], -5.85154082319921e-8*nden[0]*rflux[1], -1.1420465158197e-10*nden[0]*photden[2], -1.1420465158197e-10*nden[0]*rflux[2]]
        """
        if order not in self.ORDERS:
            raise ValueError("Invalid order: Supported orders are 0, 1, 2, 3")

        # Obtain the flux-moment equation from the density-moment equation:
        # replace each density symbol den[i] with the flux rflux[i].
        flux_map = {g.sym: self.flux[g.index] for g in self.groups}
        losses = self._photo_losses()
        emissions = self._emissions()
        dust = self.network.dust

        odes: List[sp.Expr] = [sp.Float(0.0) for _ in range(2 * self.nbands)]
        for group in self.groups:
            rate = losses[group.index]

            # Absorption by the medium: -c <α>_i den_i.  Added before the den → rflux
            # substitution it also attenuates the flux; afterwards only the density.
            absorption: sp.Expr = sp.Float(0.0)
            if group.absorption_coefficient is not None:
                absorption = -self.c * group.absorption_coefficient * group.sym
            if self.config.absorption_attenuates_flux:
                rate += absorption
                absorption = sp.Float(0.0)

            flux = rate.xreplace(flux_map)
            rate += absorption

            # Emission is isotropic, so it feeds the density equation only.
            rate += emissions[group.index]
            if group.emissivity is not None:
                rate += group.emissivity

            if dust is not None:
                density_loss, flux_loss = dust.attenuation(group)
                rate -= density_loss
                flux -= flux_loss

            ei, fi = self.ordered_index(group.index, order)
            odes[ei] = rate
            odes[fi] = flux

        return odes

    def _photo_losses(self) -> List[sp.Expr]:
        """Per band, the radiation photo reactions remove (a negative rate)."""
        net = self.network
        nden = net.symbols.ndens
        losses: List[sp.Expr] = [sp.Float(0.0) for _ in self.groups]
        for reaction in jaff_progress.track(
            net.reactions.with_type("photo"), description="Generating radiation equations"
        ):
            photo = reaction.photo
            # A custom-rate reaction is not coupled to the bands.
            if photo is not None and photo.custom_rad_rate:
                continue

            # With radiation, every other photo reaction has its bands (or the
            # build raised), so missing data is a broken invariant.
            if photo is None or not photo.bands:
                raise RuntimeError(
                    f"Photo reaction {reaction} has no per-band radiation data"
                )

            for band in photo.bands:
                # k times every reactant density (mass-action kinetics).
                loss = band.k
                for reactant in reaction.reactants.core:
                    loss *= nden[net.species[str(reactant)].index]
                losses[band.group.index] -= loss

        return losses

    def _emissions(self) -> List[sp.Expr]:
        """Per band, the radiation emitted by reactions with a spectral yield."""
        emissions: List[sp.Expr] = [sp.Float(0.0) for _ in self.groups]
        for reaction in self.network.reactions.radiation_emitters():
            yields = reaction.rad.group_yields
            for group, group_yield in zip(self.groups, yields, strict=True):
                emissions[group.index] += group_yield * reaction.rate

        return emissions

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
        This method is used by :meth:`odes`
        to place each band's pair of moment equations at the correct position
        in the generated ODE array, matching the memory layout expected by
        the downstream numerical integrator.

        Examples
        --------
        >>> from jaff import Network
        >>> net = Network(
        ...     "h_photoionization",
        ...     radiation={"bands": [1, 13.6, 100, "inf"]},
        ... )
        >>> net.radiation.nbands
        3
        >>> [net.radiation.ordered_index(1, order) for order in range(4)]
        [(2, 3), (3, 2), (1, 4), (4, 1)]
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

        Examples
        --------
        >>> import numpy as np
        >>> from jaff import Network
        >>> net = Network(
        ...     "h_photoionization",
        ...     radiation={
        ...         "bands": [1, 13.6, 100, "inf"],
        ...         "profile_index": [0, 1, 0],
        ...     },
        ... )
        >>> net.radiation.get_photden_profile(np.array([5.0, 50.0, 200.0]))
        array([4.0e-02, 2.0e-02, 2.5e-05])
        """
        lowers = np.array([float(grp.lower) for grp in self.groups])
        alphas = np.array([float(grp.profile_idx) for grp in self.groups])
        band_idx = np.clip(
            np.searchsorted(lowers, ph_energy, side="right") - 1, 0, self.nbands - 1
        )

        return ph_energy ** (alphas[band_idx] - 2)

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

    def compute(self) -> None:
        """Integrate the network's ``.jfunc`` emissivity and absorption over the bands.

        Sets :attr:`emissivity` (the ``globalEmissivity`` spectrum) and each
        group's ``emissivity`` (its integral over the band), and
        :attr:`absorption_coefficient` (the ``globalAbsorptionCoefficient``
        spectrum) and each group's ``absorption_coefficient`` (its band
        average).  Each is left ``None`` without its function.  Run by
        :class:`~jaff.core.Network` once its reactions are linked, so an
        ``rc_<N>`` in either function is reaction *N*'s final rate coefficient.
        """
        functions = self.network.functions
        resolve = self.network.reactions.resolve_rate_references

        if "globalemissivity" in functions:
            self.emissivity = functions.value("globalemissivity")
            values = self.integrate_over_groups(self.emissivity)
            for group, value in zip(self.groups, values, strict=True):
                group.emissivity = resolve(value)

        if "globalabsorptioncoefficient" in functions:
            self.absorption_coefficient = functions.value("globalabsorptioncoefficient")
            values = self.average_over_groups(self.absorption_coefficient)
            for group, value in zip(self.groups, values, strict=True):
                group.absorption_coefficient = resolve(value)

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
        if self.config.mode == "u":
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
                    spectrum, self.ph_E, grp, energy_weighted=self.config.mode == "u"
                )
            )
            for grp in self.groups
        ]

    @property
    def total_emissivity(self) -> sp.Expr | None:
        """Sum of the group emissivities, or ``None`` without a ``globalEmissivity``.

        The groups are contiguous, so this is the emissivity integrated over the
        whole band range, in the units of the density equations (erg cm⁻³ s⁻¹
        in ``"u"`` mode, photons cm⁻³ s⁻¹ in ``"nph"`` mode).
        """
        if self.emissivity is None:
            return None

        return sum((grp.emissivity for grp in self.groups), sp.S.Zero)  # type: ignore
