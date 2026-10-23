# ABOUTME: PhotoReaction: a photo reaction's component, held as Reaction.photo: its
# ABOUTME: resolved photo settings, cross sections, band rates and shielding tables

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Dict, List, cast

from astropy import units as u
from sympy import Expr

from ...errors import ParserError
from ._helper import to_float_or_none
from .config import PhotoConfig

if TYPE_CHECKING:
    import matplotlib.pyplot as plt
    import pandas as pd

    from ...physics import XsecProps
    from .reaction import Reaction
    from .reaction_props import GroupedPhotoReactionProps


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
        when no database has an entry.  ``xsec.total`` holds the cross
        section averaged over the whole radiation spectrum (cm²).
    custom_rad_rate : bool
        ``True`` when the rate comes from a ``.jfunc`` ``chemRate<N>``
        function rather than from the cross sections; such a reaction is not
        coupled to the radiation bands.
    rgrouped : list[GroupedPhotoReactionProps] or None
        One entry per radiation band, in ``net.radiation.groups`` order, with
        that band's loss coefficient ``k``, band-average ``xsec`` and
        ``xsec_frac``.  Set when a radiation field is configured; ``None``
        otherwise, and for a ``custom_rad_rate`` reaction.  See
        :attr:`band_xsecs`.
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
        self.custom_rad_rate: bool = False
        self.rgrouped: List[GroupedPhotoReactionProps] | None = None
        self.shielding_tables: Dict[str, Any] = {}

    def compute(self) -> Expr | None:
        """Resolve the cross sections and, with radiation, the band rates.

        Returns
        -------
        sympy.Expr or None
            The reaction's new rate coefficient (the sum of its per-band rate
            coefficients) when a radiation field sets it; ``None`` when the
            reaction keeps its own rate.

        Raises
        ------
        ParserError
            If a ``pi_database`` is set for the reaction itself while its
            cross section comes from a ``photoCrossSection<N>`` function, or
            radiation needs a cross section and no database has one.
        """
        reaction = self.reaction
        if reaction.aux_funcs.chemrate is not None:
            self.custom_rad_rate = not self._overrides_custom_rate()

        self._set_xsec()

        if reaction.net.radiation is None:
            return None

        return self._set_band_rates()

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
        from ...physics import Photochemistry

        reaction, net = self.reaction, self.reaction.net
        if net._photochemistry is None:
            net._photochemistry = Photochemistry(net)

        custom = reaction.aux_funcs.photocrosssection
        if custom is not None and not self._overrides_custom_cross_section():
            own = net.reactions.config.own(reaction.serialized)
            if own is not None and own.photo is not None and own.photo.pi_database:
                raise ParserError(
                    f"pi_database cannot be set for {reaction}: its cross section "
                    "comes from a photoCrossSection function"
                )

            self.xsec = net._photochemistry.get_custom_xsec(
                reaction, cast(Expr, custom.expr)
            )
            return

        # A cross section is mandatory only when radiation sets the rate.
        self.xsec = net._photochemistry.get_xsec(
            reaction,
            required=net.radiation is not None and not self.custom_rad_rate,
        )

    def _set_band_rates(self) -> Expr | None:
        """Build the per-band rate data; return the summed rate coefficient."""
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
            return None

        k_tot, xsec_tot, grouped = net.radiation.get_reaction_rate_props(reaction)

        assert self.xsec is not None
        self.xsec.total = xsec_tot
        self.rgrouped = grouped

        return k_tot

    @property
    def band_xsecs(self) -> pd.DataFrame:
        """Band-averaged cross sections of the reaction, one row per band.

        Pairs each per-band entry of :attr:`rgrouped` (populated when a
        radiation field is configured, e.g. through the ``radiation``
        argument of :class:`~jaff.core.Network` or ``[network.radiation]``)
        with the :class:`~jaff.physics.photo_reactions._radiation.RadiationGroup`
        it was computed for.  A photo reaction whose rate comes from a
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
              spectrum-averaged one (``xsec.total``); ``0`` when that is zero.
              Not normalised to sum to 1 over the bands.

            Empty (with the columns above) when the reaction contributes to no
            band, e.g. no radiation field is configured.

        Notes
        -----
        The rows are ordered by ascending band index: :attr:`rgrouped` is
        built by iterating ``net.radiation.groups`` in order, so the two are
        paired by position.

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
        radiation = self.reaction.net.radiation
        grouped = self.rgrouped
        groups = radiation.groups if grouped and radiation is not None else []
        rows = [
            {
                "lower": to_float_or_none(group.lower),
                "upper": to_float_or_none(group.upper),
                "eavg": (
                    None
                    if group.eavg is None
                    else to_float_or_none(group.eavg * u.erg.to(u.eV))
                ),
                "xsec": to_float_or_none(props.xsec),
                "xsec_frac": to_float_or_none(props.xsec_frac),
            }
            for group, props in zip(groups, grouped or [], strict=True)
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
