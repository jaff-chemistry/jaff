# ABOUTME: ReactionEnergy (reaction.energy): a reaction's energy change per event and its
# ABOUTME: term of dE/dt; a dE in ph_E is averaged over the band photons absorbed or emitted
"""
The energy a reaction adds to the gas.

Every reaction holds a :class:`ReactionEnergy` as ``reaction.energy``.  Its
:attr:`~ReactionEnergy.per_event` is the energy change one event adds to the gas
(erg; negative cools): the reaction's ``.jfunc`` ``deltaE<N>``, else its
energy-change plugin, else the network file's value, else 0.
:meth:`~ReactionEnergy.change` is the reaction's term of the chemical heating
rate ``dE/dt``.

A ``per_event`` may depend on the photon energy ``ph_E`` (eV): the energy the
event adds when the photon it involves has energy ``ph_E``.  It is then averaged
over the radiation bands with the distribution of those photons: for a photo
reaction the photons it absorbs, per band and weighted as the band rates are
(photodecay cross section times the band's photon spectrum); for a reaction with
a ``spectralYield`` the photons it emits into the bands, ``Y(E)/E``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, List, Mapping, cast

from sympy import Basic, Expr, Float

from ...errors import ParserError

if TYPE_CHECKING:
    from .reaction import Reaction

__all__ = ["ReactionEnergy"]


class ReactionEnergy:
    """A reaction's energy change per event and its term of ``dE/dt``.

    Attributes
    ----------
    reaction : Reaction
        The reaction.
    per_event : sympy.Expr
        Energy one event adds to the gas (erg; negative cools), as given: the
        ``.jfunc`` ``deltaE<N>``, else the energy-change plugin, else the
        network file's value, else 0.  May depend on ``ph_E``.
    photon : str
        What ``ph_E`` in :attr:`per_event` stands for: ``"none"`` (it does not
        appear), ``"absorbed"`` (the photon a photo reaction absorbs) or
        ``"emitted"`` (a photon the reaction's ``spectralYield`` emits).
    band_values : list of sympy.Expr or None
        With ``photon == "absorbed"``: :attr:`per_event` averaged over each
        radiation band (``net.radiation.groups`` order) with the rate's weight;
        else ``None``.
    mean : sympy.Expr or None
        With ``photon == "emitted"``: :attr:`per_event` averaged over the photons
        the reaction emits into the bands (``Radiation.emitted_average``); else
        ``None``.

    Examples
    --------
    >>> from jaff import Network
    >>> net = Network("h_photoionization")
    >>> energy = net.reactions[0].energy
    >>> energy.per_event, energy.photon
    (6.40000000000000e-12, 'none')
    >>> energy.change() == energy.per_event * net.reactions[0].rate.expanded
    True
    """

    def __init__(self, reaction: Reaction, per_event: Expr):
        """Hold *reaction*'s energy change *per_event*; :meth:`compute` does the rest.

        Parameters
        ----------
        reaction : Reaction
            The reaction.
        per_event : sympy.Expr
            Energy one event adds to the gas (erg), as resolved by the reaction.
        """
        self.reaction: Reaction = reaction
        self.per_event: Expr = per_event
        self.photon: str = "none"
        self.band_values: List[Expr] | None = None
        self.mean: Expr | None = None

    def compute(self) -> None:
        """Decide what ``ph_E`` stands for and average :attr:`per_event` over the bands.

        Run once the reaction's ``photo`` and ``rad`` components are computed.

        Raises
        ------
        ParserError
            If :attr:`per_event` depends on ``ph_E`` and the photon is not
            defined: there is no radiation field, the reaction neither absorbs
            (a photo reaction whose rate the bands set) nor emits (a
            ``spectralYield``), or it does both.
        """
        reaction = self.reaction
        if reaction.net.symbols.ph_E not in self.per_event.free_symbols:
            return

        photo = reaction.photo
        absorbs = photo is not None and bool(photo.bands)
        emits = reaction.rad.spectral_yield is not None
        reason = None
        if reaction.net.radiation is None:
            reason = "the network has no radiation field"
        elif absorbs and emits:
            reason = (
                "the reaction absorbs and emits photons (a photo reaction with a "
                "spectralYield), so ph_E could be either photon"
            )
        elif not (absorbs or emits):
            reason = (
                "the reaction neither absorbs band photons (a photo reaction whose "
                "rate the bands set) nor emits any (a spectralYield)"
            )
        if reason is not None:
            raise ParserError(
                f"Reaction {reaction.verbatim} (index {reaction.index}): its energy "
                "change depends on the photon energy ph_E, which is averaged over the "
                f"band photons the reaction absorbs or emits, but {reason}"
            )

        radiation = reaction.net.radiation
        assert radiation is not None
        self.photon = "absorbed" if absorbs else "emitted"
        try:
            if absorbs:
                assert photo is not None
                energy, decay, _ = photo.cross_section_data(radiation.ph_E)
                self.band_values = radiation.absorbed_average(
                    self.per_event, decay, energy
                )
            else:
                spectral_yield = reaction.rad.spectral_yield
                assert spectral_yield is not None
                self.mean = radiation.emitted_average(self.per_event, spectral_yield)
        except ValueError as err:
            raise ParserError(
                f"Reaction {reaction.verbatim} (index {reaction.index}): cannot average "
                f"its energy change over the {self.photon} photons: {err}"
            ) from err

    def change(self, form: str = "expanded") -> Expr:
        """The reaction's term of the chemical heating rate ``dE/dt`` (erg cm⁻³ s⁻¹).

        ``per_event`` times the rate in *form*; with ``photon == "emitted"``,
        :attr:`mean` times the rate; with ``photon == "absorbed"``,
        ``Σ_b rate_b · <dE>_b`` times the reactant densities (the band rate
        terms ``reaction.photo.bands[b].rate``, in either form).

        Parameters
        ----------
        form : str, optional
            ``"expanded"`` (default) or ``"reduced"`` (see :class:`Rate`).

        Returns
        -------
        sympy.Expr
            The term.
        """
        rate = self.reaction.rate
        if self.mean is not None:
            return self.mean * rate.form(form)
        if self.band_values is None:
            return self.per_event * rate.form(form)

        bands = self.reaction.photo.bands if self.reaction.photo is not None else None
        assert bands is not None
        heating = sum(
            (band.rate * value for band, value in zip(bands, self.band_values)),
            Float(0.0),
        )

        return heating * rate.densities

    def expressions(self) -> list[Expr]:
        """The expressions the component holds: :attr:`per_event`, :attr:`band_values`, :attr:`mean`.

        Returns
        -------
        list[Expr]
            The expressions.
        """
        extra = [self.mean] if self.mean is not None else []

        return [self.per_event, *(self.band_values or []), *extra]

    def substitute(self, mapping: Mapping[Basic, Basic]) -> None:
        """Apply *mapping* (``xreplace``) to every expression the component holds.

        Parameters
        ----------
        mapping : Mapping[Basic, Basic]
            Old symbol or expression to its replacement.
        """
        self.per_event = cast(Expr, self.per_event.xreplace(mapping))
        if self.band_values is not None:
            self.band_values = [cast(Expr, v.xreplace(mapping)) for v in self.band_values]
        if self.mean is not None:
            self.mean = cast(Expr, self.mean.xreplace(mapping))
