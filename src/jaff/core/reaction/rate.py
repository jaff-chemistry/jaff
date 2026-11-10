# ABOUTME: Rate: a reaction's rate per unit volume in two forms, the coefficient written
# ABOUTME: out (expanded) or as its krate[r] symbol (reduced), for the equation builders
"""The rate of one reaction, ``reaction.rate``."""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

from sympy import Expr, Indexed, Integer

if TYPE_CHECKING:
    from .reaction import Reaction

RateForm = Literal["expanded", "reduced"]

RATE_FORMS: tuple[str, ...] = ("expanded", "reduced")
"""The forms of a :class:`Rate`: the coefficient written out, or its ``krate`` symbol."""


def check_form(form: str) -> None:
    """Raise unless *form* is one of :data:`RATE_FORMS`.

    Parameters
    ----------
    form : str
        The form to check.

    Raises
    ------
    ValueError
        If *form* is not ``"expanded"`` or ``"reduced"``.
    """
    if form not in RATE_FORMS:
        raise ValueError(
            f"Invalid rate form {form!r}; expected one of: {', '.join(RATE_FORMS)}"
        )


class Rate:
    """Rate per unit volume of one reaction, ``k · Π nden[reactants]``.

    Reached as ``reaction.rate``.  The rate coefficient ``k`` appears in one of two
    forms: :attr:`expanded` writes it out (``reaction.rate_coefficient``),
    :attr:`reduced` uses the reaction's symbol ``krate[r]``, *r* being its position in
    ``net.reactions``.  Equations built from the reduced form, with the definitions
    ``krate[r] = coefficient`` kept apart, evaluate every coefficient once.

    Attributes
    ----------
    reaction : Reaction
        The reaction whose rate this is.

    Examples
    --------
    >>> from jaff import Network
    >>> net = Network("h_photoionization")
    >>> rate = net.reactions["H+.e-__H"].rate
    >>> rate.expanded
    1.65941781598291e-10*nden[1]*nden[2]/tgas**0.7
    >>> rate.reduced
    krate[1]*nden[1]*nden[2]
    """

    def __init__(self, reaction: Reaction) -> None:
        """Hold *reaction*; the forms are built on first use.

        Parameters
        ----------
        reaction : Reaction
            The reaction whose rate this is.
        """
        self.reaction: Reaction = reaction
        self._factors: list[Indexed] | None = None
        self._expanded: tuple[Expr, Expr] | None = None

    @property
    def coefficient(self) -> Expr:
        """The rate coefficient written out, ``reaction.rate_coefficient``."""
        return self.reaction.rate_coefficient

    @property
    def symbol(self) -> Indexed:
        """The coefficient's symbol, ``krate[r]`` with *r* the reaction's catalogue position."""
        return self.reaction.net.symbols.krates[self.reaction.catalogue_index]

    @property
    def _density_factors(self) -> list[Indexed]:
        """``nden[i]`` of each core reactant, once per occurrence."""
        if self._factors is None:
            net = self.reaction.net
            nden = net.symbols.ndens
            self._factors = [
                nden[net.species[s.name].index] for s in self.reaction.reactants.core
            ]

        return self._factors

    @property
    def densities(self) -> Expr:
        """Product of the core reactants' number densities (``1`` without any).

        Pseudo-species such as ``_PHOTON`` or ``_CR`` carry no density; a reactant
        appearing twice counts twice.
        """
        return self._times_densities(Integer(1))

    def _times_densities(self, factor: Expr) -> Expr:
        """*factor* multiplied by one density at a time.

        One at a time, as SymPy then folds the numbers of a partly unevaluated
        coefficient the same way for every rate.
        """
        for density in self._density_factors:
            factor *= density

        return factor

    @property
    def expanded(self) -> Expr:
        """``coefficient · densities``; cached until the rate coefficient changes."""
        coefficient = self.coefficient
        if self._expanded is None or self._expanded[0] is not coefficient:
            self._expanded = (coefficient, self._times_densities(coefficient))

        return self._expanded[1]

    @property
    def reduced(self) -> Expr:
        """``symbol · densities``: the rate with the coefficient as ``krate[r]``."""
        return self._times_densities(self.symbol)

    def form(self, form: str = "expanded") -> Expr:
        """The rate in *form*: :attr:`expanded` or :attr:`reduced`.

        Parameters
        ----------
        form : {"expanded", "reduced"}, optional
            Which form, by default ``"expanded"``.

        Returns
        -------
        Expr

        Raises
        ------
        ValueError
            If *form* is not one of :data:`RATE_FORMS`.
        """
        check_form(form)

        return self.expanded if form == "expanded" else self.reduced
