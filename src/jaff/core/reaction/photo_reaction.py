# ABOUTME: PhotoReaction: a photo reaction's component, held as Reaction.photo: its
# ABOUTME: resolved photo settings, cross sections, band rates and shielding tables

from __future__ import annotations

from dataclasses import replace
from typing import TYPE_CHECKING, Any, Dict, List, Mapping, cast

import numpy as np
from astropy import units as u
from sympy import Basic, Expr, Float, Symbol

from ...config.base import dotted
from ...errors import ParserError
from ._helper import to_float_or_none
from .config import PhotoConfig

if TYPE_CHECKING:
    import matplotlib.pyplot as plt
    import pandas as pd

    from ...physics import PhotoBand, Radiation, XsecProps
    from .reaction import Reaction


class PhotoReaction:
    """A photo reaction's component, held as ``Reaction.photo`` (``None`` otherwise).

    Built by :class:`~jaff.core.reaction.Reaction` for every reaction of type
    ``"photo"`` from the reaction's resolved photo settings
    (``reaction.config.photo``).  It decides where the rate comes from (a
    ``.jfunc`` ``chemRate<N>`` or the cross sections), looks up the cross
    sections (a ``photoCrossSection<N>`` function or the databases) and, with
    a radiation field, builds the per-band rate data and the shielding.

    Attributes
    ----------
    reaction : Reaction
        The reaction it belongs to.
    config : PhotoConfig
        The reaction's resolved photo settings (read-only): ``pi_database``,
        ``override_custom_rate``, ``override_custom_cross_section``,
        ``use_proxy_photoreaction`` and ``shielding``.
    xsec : XsecProps or None
        Cross sections of the reaction's single decay channel, or ``None``
        when no database has an entry.
    xsec_total : float or sympy.Expr or None
        The photodecay cross section averaged over the whole radiation
        spectrum (cm²); ``None`` without a radiation field, and for a
        ``custom_rad_rate`` reaction.  An expression when a
        ``photoCrossSection<N>`` uses other symbols (e.g. ``tgas``).
    custom_rad_rate : bool
        ``True`` when the rate comes from a ``.jfunc`` ``chemRate<N>``
        function rather than from the cross sections; such a reaction is not
        coupled to the radiation bands.
    bands : list[PhotoBand] or None
        One entry per radiation band, in ``net.radiation.groups`` order, with
        its ``group``, the band's loss coefficient ``k``, band-average ``xsec``
        and ``xsec_frac``.  Set when a radiation field is configured; ``None``
        otherwise, and for a ``custom_rad_rate`` reaction.  See
        :attr:`band_xsecs`.
    band_rate : sympy.Expr or None
        The sum of the per-band rate coefficients: the reaction's rate
        coefficient when a radiation field sets it.  ``None`` otherwise, and
        for a ``custom_rad_rate`` reaction.
    shielding_factor : sympy.Expr or None
        The factor the reaction's shielding model (``config.shielding``)
        multiplies its band rates by; ``None`` when it is not shielded, and for
        a ``custom_rad_rate`` reaction.
    shielding_tables : dict[str, Any]
        Data tables the reaction's shielding model built for the generated
        code to interpolate over, keyed by file name (e.g.
        ``shielding_<reaction>.hdf5``); written out by
        :meth:`~jaff.core.network.Network.write_shielding_tables`.  Empty when
        the reaction is not shielded or its model needs no table.

    Examples
    --------
    >>> from jaff import Network
    >>> net = Network("h_photoionization")
    >>> photo = net.reactions.photo_reactions()[0].photo
    >>> photo.config.pi_database, photo.custom_rad_rate, photo.shielding_tables
    ('norad', False, {})
    >>> photo.xsec.database
    'norad'
    """

    def __init__(
        self,
        reaction: Reaction,
        config: PhotoConfig,
    ):
        """Store the component's state; :meth:`compute` fills the rest.

        Parameters
        ----------
        reaction : Reaction
            The photo reaction.
        config : PhotoConfig
            Its resolved photo settings (``reaction.config.photo``).
        """
        self.reaction: Reaction = reaction
        self.config: PhotoConfig = config
        self.xsec: XsecProps | None = None
        self.xsec_total: float | Expr | None = None
        self.custom_rad_rate: bool = False
        self.bands: List[PhotoBand] | None = None
        self.band_rate: Expr | None = None
        self.shielding_factor: Expr | None = None
        self.shielding_tables: Dict[str, Any] = {}

    def compute(self) -> None:
        """Resolve the cross sections and, with radiation, the band rates.

        Sets :attr:`xsec` and, when a radiation field sets the rate, the
        shielding (:attr:`shielding_factor`, :attr:`shielding_tables`),
        :attr:`bands` and :attr:`band_rate`.

        Raises
        ------
        ParserError
            If a ``pi_database`` is set for the reaction itself while its
            cross section comes from a ``photoCrossSection<N>`` function,
            radiation needs a cross section and no database has one, or the
            reaction is shielded without a radiation field.
        """
        reaction = self.reaction
        if reaction.function("chemRate") is not None:
            self.custom_rad_rate = not self._overrides_custom_rate()

        self._set_xsec()

        if reaction.net.radiation is None:
            self._check_no_shielding()
            return

        self._set_band_rates()

    def _check_no_shielding(self) -> None:
        """Reject a shielding setting without a radiation field (it would do nothing)."""
        if self.config.shielding is None:
            return

        serialized = self.reaction.serialized
        own = self.reaction.net.reactions.config.own(serialized)
        path = (
            ("network", "reactions", "reaction", serialized, "photo", "shielding")
            if own is not None and own.photo is not None and own.photo.shielding
            else ("network", "reactions", "photo", "shielding")
        )
        raise ParserError(
            f"{dotted(path)}: shielding of {self.reaction} needs a radiation field "
            "([network.radiation])"
        )

    def _overrides_custom_rate(self) -> bool:
        """Whether radiation recomputes the reaction's ``chemRate`` rate.

        From ``config.override_custom_rate``; without radiation there is no
        cross-section rate to use, so the ``chemRate`` is always kept.
        """
        if self.reaction.net.radiation is None:
            return False

        return bool(self.config.override_custom_rate)

    def _overrides_custom_cross_section(self) -> bool:
        """Whether the reaction ignores its ``photoCrossSection``.

        From ``config.override_custom_cross_section``; like
        ``override_custom_rate``, it has no effect without radiation.
        """
        if self.reaction.net.radiation is None:
            return False

        return bool(self.config.override_custom_cross_section)

    def _set_xsec(self) -> None:
        """Set :attr:`xsec` from a ``photoCrossSection<N>`` or the databases."""
        reaction, net = self.reaction, self.reaction.net
        photochemistry = net.reactions.photochemistry

        custom = reaction.function("photoCrossSection")
        if custom is not None and not self._overrides_custom_cross_section():
            self._check_no_rate_reference(custom)
            own = net.reactions.config.own(reaction.serialized)
            if own is not None and own.photo is not None and own.photo.pi_database:
                raise ParserError(
                    f"pi_database cannot be set for {reaction}: its cross section "
                    "comes from a photoCrossSection function"
                )

            self.xsec = photochemistry.get_custom_xsec(reaction, custom)
            return

        # A cross section is mandatory only when radiation sets the rate.
        self.xsec = photochemistry.get_xsec(
            reaction,
            required=net.radiation is not None and not self.custom_rad_rate,
        )

    def _check_no_rate_reference(self, xsec: Expr) -> None:
        """Reject ``rc_<N>`` in the reaction's ``photoCrossSection`` function.

        A cross section is a per-photon property of the absorber; it cannot
        depend on another reaction's rate coefficient.
        """
        references = sorted(
            str(s) for s in xsec.free_symbols if str(s).lower().startswith("rc_")
        )
        if references:
            raise ParserError(
                f"photoCrossSection{self.reaction.index} uses {', '.join(references)}: "
                "a cross section cannot depend on another reaction's rate"
            )

    def _set_band_rates(self) -> None:
        """Build the per-band rate data and the summed rate coefficient."""
        reaction, net = self.reaction, self.reaction.net
        assert net.radiation is not None

        # A custom rate is not driven by the transported field, so the reaction
        # is not coupled to the bands and gets no band data.
        if self.custom_rad_rate:
            reaction.logger.warning(
                f"Photo reaction [cyan]{reaction}[/] (index {reaction.index}) uses "
                f"the custom rate chemRate{reaction.index}: it is not coupled to "
                "the radiation bands and removes no photons. Supply a "
                f"photoCrossSection{reaction.index} instead, or set "
                "override_custom_rate = true to compute its rate from the "
                "database cross sections."
            )
            return

        self._set_shielding()
        self.band_rate, self.xsec_total, self.bands = self.band_rates(net.radiation)

    def _set_shielding(self) -> None:
        """Run the reaction's shielding model; keep its factor and tables."""
        shielding = self.config.shielding
        if shielding is None:
            return

        from ...physics.photo_reactions.shielding import shielding_function

        plugin = shielding_function(shielding.type, self.reaction.serialized)
        result = plugin.get_shielding(self.reaction)
        self.shielding_factor = result.factor
        self.shielding_tables = dict(result.tables)

    def band_rates(
        self, radiation: Radiation
    ) -> tuple[Expr, float | Expr, list[PhotoBand]]:
        """Assemble the reaction's band-averaged rate data in *radiation*.

        From :attr:`xsec` (tabulated arrays, or an expression in ``ph_E`` for
        Verner and ``photoCrossSection`` cross sections), for each band of
        *radiation*:

        1. the photon-number-weighted band-average cross section
           ``<σ>_i = ∫ σ n dE / ∫ n dE``;
        2. the photon-counting rate coefficient ``c · n_i · <σ>_i``, with
           ``n_i`` the band's photon density
           (:meth:`~jaff.physics.Radiation.photon_density`: ``den[i]``, or
           ``den[i] / <E>_i`` in ``"u"`` mode);
        3. the band loss coefficient ``c · den[i] · <σ>_i^loss``
           (:meth:`~jaff.physics.Radiation.loss_cross_section`: energy-weighted
           in ``"u"`` mode, as each absorbed photon takes its own energy);
        4. one :class:`~jaff.physics.PhotoBand` with the loss coefficient (as
           ``k``), ``<σ>_i`` and its ratio to the spectrum average
           (``xsec_frac``).

        Both coefficients carry :attr:`shielding_factor`, if any.
        Nothing is stored; :meth:`compute` stores the result.

        Parameters
        ----------
        radiation : Radiation
            The radiation field.

        Returns
        -------
        tuple[sympy.Expr, float or sympy.Expr, list[PhotoBand]]
            ``(band_rate, xsec_total, bands)``: the rate coefficient summed
            over the bands (s⁻¹), the cross section averaged over the whole
            spectrum (cm²), and one band in ``radiation.groups`` order.

        Raises
        ------
        ValueError
            If there is no cross section (no :attr:`xsec`, or neither a
            ``photodecay`` array nor a ``photodecay_expr``).

        Notes
        -----
        The photon-number spectrum used for averaging is the piecewise power
        law ``n(E) ∝ E^(α_i - 2)`` in band *i* (``α_i = group.profile_idx``).

        A tabulated cross section with a separate photo-absorption table
        (``xsec.process.absorption``) uses the absorption cross section for the
        band data (``k`` and ``xsec``) and rescales each band's rate term by
        the decay-to-absorption ratio of the band averages, so the rate still
        follows the photodecay cross section.  The spectrum average is always
        that of the photodecay cross section.

        Examples
        --------
        >>> from jaff import Network
        >>> net = Network(
        ...     "h_photoionization",
        ...     radiation={"bands": [1, 13.6, 100, "inf"]},
        ... )
        >>> photo = net.reactions.photo_reactions()[0].photo
        >>> band_rate, xsec_total, bands = photo.band_rates(net.radiation)
        >>> band_rate
        5.85154082319921e-8*photden[1] + 1.1420465158197e-10*photden[2]
        >>> float(xsec_total)
        1.2403886126668026e-19
        >>> bands[1].k, float(bands[1].xsec)
        (5.85154082319921e-8*photden[1], 1.9518639202054948e-18)
        """
        from ...physics.photo_reactions._photo_band import PhotoBand

        energy, decay, absorption = self._cross_section_data(radiation.ph_E)
        profile = (
            radiation.get_photden_profile(energy)
            if isinstance(energy, np.ndarray)
            else None
        )
        xsec_total = radiation.get_total_average_cross_section(decay, energy, profile)
        shielding = self.shielding_factor

        band_rate: Expr = Float(0.0)
        bands: list[PhotoBand] = []
        for group in radiation.groups:
            decay_avg = radiation.get_group_average_cross_section(
                decay, energy, group, profile
            )
            band_avg = (
                radiation.get_group_average_cross_section(
                    absorption, energy, group, profile
                )
                if absorption is not None
                else decay_avg
            )
            loss_avg = radiation.loss_cross_section(
                decay if absorption is None else absorption,
                energy,
                group,
                profile,
                band_avg,
            )

            rate = radiation.c * radiation.photon_density(group) * band_avg
            k = radiation.c * radiation.den[group.index] * loss_avg
            if shielding is not None:
                rate *= shielding
                k *= shielding

            bands.append(
                PhotoBand(
                    group=group,
                    k=k,
                    xsec=band_avg,
                    xsec_frac=band_avg / xsec_total if xsec_total != 0.0 else 0.0,
                )
            )
            # With an absorption table, the rate follows the decay cross section.
            if absorption is not None:
                rate *= decay_avg / band_avg if band_avg != 0.0 else 0.0
            band_rate += rate

        return band_rate, xsec_total, bands

    def _cross_section_data(
        self, ph_E: Symbol
    ) -> tuple[np.ndarray | Symbol, np.ndarray | Expr, np.ndarray | None]:
        """``(energy, decay, absorption)`` of :attr:`xsec` for the band integrals.

        *energy* is the tabulated grid (eV), or *ph_E* for a symbolic cross
        section; *absorption* is the photo-absorption table, when the cross
        section has one.
        """
        xsec = self.xsec
        if xsec is None:
            raise ValueError(
                f"Photo reaction {self.reaction} has no cross section (photo.xsec is None)"
            )

        if xsec.photodecay_expr is not None:
            return ph_E, xsec.photodecay_expr, None

        if xsec.photodecay is None:
            raise ValueError(
                f"Photodecay cross-section cannot be None for reaction: {self.reaction}"
            )

        assert isinstance(xsec.photon_energy, np.ndarray)
        absorption = xsec.photo_absorption if xsec.process.absorption else None

        return xsec.photon_energy, xsec.photodecay, absorption

    def expressions(self) -> list[Expr]:
        """The expressions the component holds.

        Its :attr:`band_rate`, :attr:`shielding_factor`, each band's rate coefficient ``k`` and symbolic
        band cross section, and a symbolic :attr:`xsec_total`.

        Returns
        -------
        list[Expr]
            The expressions; empty without radiation.
        """
        exprs: list[Expr] = [] if self.band_rate is None else [self.band_rate]
        if self.shielding_factor is not None:
            exprs.append(self.shielding_factor)
        for band in self.bands or []:
            exprs.append(band.k)
            if isinstance(band.xsec, Expr):
                exprs.append(band.xsec)
        if isinstance(self.xsec_total, Expr):
            exprs.append(self.xsec_total)

        return exprs

    def substitute(self, mapping: Mapping[Basic, Basic]) -> None:
        """Apply *mapping* (``xreplace``) to every expression the component holds.

        Parameters
        ----------
        mapping : Mapping[Basic, Basic]
            Old symbol or expression to its replacement.
        """
        if self.band_rate is not None:
            self.band_rate = cast(Expr, self.band_rate.xreplace(mapping))
        if self.shielding_factor is not None:
            self.shielding_factor = cast(Expr, self.shielding_factor.xreplace(mapping))
        if self.bands is not None:
            self.bands = [_substitute_band(band, mapping) for band in self.bands]
        if isinstance(self.xsec_total, Expr):
            self.xsec_total = cast(Expr, self.xsec_total.xreplace(mapping))

    @property
    def band_xsecs(self) -> pd.DataFrame:
        """Band-averaged cross sections of the reaction, one row per band.

        One row per entry of :attr:`bands` (populated when a radiation field is
        configured, e.g. through the ``radiation`` argument of
        :class:`~jaff.core.Network` or ``[network.radiation]``), with the edges
        and mean energy of its ``group``.  A photo reaction whose rate comes from a
        ``.jfunc`` ``chemRate<N>`` is not coupled to the bands and has no band
        data, so its frame is empty.  Intended as the data source for band bar
        plots (the ``bands`` option of :meth:`plot_xsecs`).

        Returns
        -------
        pandas.DataFrame
            One row per radiation band the reaction contributes to, with
            columns:

            - ``lower`` : lower band edge in eV.
            - ``upper`` : upper band edge in eV (``inf`` for an open top band).
            - ``eavg``  : photon-number-weighted band-average energy in eV.
            - ``xsec``  : photon-number-weighted band-average cross section in
              cm² (``NaN`` when it is not a number, e.g. a
              ``photoCrossSection<N>`` that depends on ``tgas``).
            - ``xsec_frac`` : the band-average cross section divided by the
              spectrum-averaged one (:attr:`xsec_total`); ``0`` when that is zero.
              Not normalised to sum to 1 over the bands.

            Empty (with the columns above) when the reaction contributes to no
            band, e.g. no radiation field is configured.

        Notes
        -----
        The rows are ordered by ascending band index, as :attr:`bands` is.

        Examples
        --------
        >>> net = Network(
        ...     "h_photoionization",
        ...     radiation={"bands": [1, 13.6, 100, "inf"]},
        ... )
        >>> photo = net.reactions.photo_reactions()[0].photo
        >>> photo.band_xsecs.columns.tolist()
        ['lower', 'upper', 'eavg', 'xsec', 'xsec_frac']
        >>> len(photo.band_xsecs)
        3
        """
        import pandas as pd

        columns = ["lower", "upper", "eavg", "xsec", "xsec_frac"]
        rows = [
            {
                "lower": to_float_or_none(band.group.lower),
                "upper": to_float_or_none(band.group.upper),
                "eavg": (
                    None
                    if band.group.eavg is None
                    else to_float_or_none(band.group.eavg * u.erg.to(u.eV))
                ),
                "xsec": to_float_or_none(band.xsec),
                "xsec_frac": to_float_or_none(band.xsec_frac),
            }
            for band in self.bands or []
        ]

        return pd.DataFrame(rows, columns=columns)

    def plot_xsecs(self, **kwargs: Any) -> tuple[plt.Figure, Any] | None:
        """Plot the reaction's cross sections against photon energy or wavelength.

        Forwards to :func:`jaff.plotting.plot_xsecs` with the reaction as the
        only item; every ``plot_xsecs`` keyword is accepted (``grid``,
        ``processes``, ``energy_unit``, ``xsec_unit``, ``xscale``, ``yscale``,
        ``bands``, ``style``, ``palette``, ``save``, ``show``, ``close`` ...).
        The panel is titled with the reaction equation and the curves are
        labelled by process.

        Returns
        -------
        tuple[matplotlib.figure.Figure, numpy.ndarray] or None
            ``(fig, axes)``, or ``None`` (with a log message) if :attr:`xsec`
            is ``None`` or no requested process has data.

        Examples
        --------
        >>> photo = net.reactions.photo_reactions()[0].photo
        >>> photo.plot_xsecs()  # all processes, eV vs Mb
        >>> photo.plot_xsecs(processes="photodecay", yscale="linear")
        >>> photo.plot_xsecs(energy_unit="nm", xsec_unit="cm^2", bands=True)
        """
        from ...plotting import plot_xsecs

        return plot_xsecs(self.reaction, **kwargs)


def _substitute_band(band: PhotoBand, mapping: Mapping[Basic, Basic]) -> PhotoBand:
    """*band* with *mapping* applied (``xreplace``) to its ``k`` and symbolic ``xsec``."""
    xsec = band.xsec
    if isinstance(xsec, Expr):
        xsec = cast(Expr, xsec.xreplace(mapping))

    return replace(band, k=cast(Expr, band.k.xreplace(mapping)), xsec=xsec)
