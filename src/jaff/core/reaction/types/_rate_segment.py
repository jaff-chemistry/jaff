# ABOUTME: RateSegment: one rate expression valid over one temperature range
# ABOUTME: Frozen dataclass of rate, tmin and tmax

from dataclasses import dataclass

from sympy import Expr


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
    """

    rate: Expr
    tmin: float | None
    tmax: float | None
