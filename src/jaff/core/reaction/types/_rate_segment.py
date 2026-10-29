# ABOUTME: RateSegment: one rate expression valid over one temperature range
# ABOUTME: Frozen dataclass of rate, tmin, tmax and the source line

from __future__ import annotations

from dataclasses import dataclass

from sympy import Expr

from ._reaction_source import ReactionSource


@dataclass(frozen=True)
class RateSegment:
    """A reaction's rate expression valid over one temperature range.

    Immutable: change a reaction's rates through its
    :class:`RateSegments` (:meth:`RateSegments.map_rates`, :meth:`RateSegments.add`).

    Attributes
    ----------
    rate : Expr
        SymPy rate-coefficient expression for this range.
    tmin : float or None
        Lower temperature bound in Kelvin (``None`` = unbounded below).
    tmax : float or None
        Upper temperature bound in Kelvin (``None`` = unbounded above).
    source : ReactionSource or None
        The network-file line the rate came from; ``None`` for a rate not
        read from a file.
    """

    rate: Expr
    tmin: float | None
    tmax: float | None
    source: ReactionSource | None = None
