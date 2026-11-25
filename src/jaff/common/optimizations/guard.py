# ABOUTME: guard(): a CSE temporary read only inside Piecewise branches is computed only
# ABOUTME: under their conditions (else 0.0), so it is never evaluated out of its range
"""Guards keeping CSE temporaries inside the branches that read them."""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Sequence
from typing import cast

from sympy import (
    ITE,
    And,
    Basic,
    Expr,
    Float,
    GreaterThan,
    Integer,
    LessThan,
    Not,
    Or,
    Piecewise,
    StrictGreaterThan,
    StrictLessThan,
    false,
    true,
)
from sympy.core.relational import Relational
from sympy.core.sorting import default_sort_key
from sympy.logic.boolalg import Boolean

Temporaries = list[tuple[Basic, Expr]]
Path = frozenset[Boolean]
"""Conditions that all hold where an expression is evaluated."""

Interval = tuple[float, bool, float, bool]
"""``(low, low closed, high, high closed)``."""

Intervals = list[Interval]
"""Disjoint intervals, sorted."""

# Each comparison with the variable on the left, as (is lower bound, closed).
_BOUNDS: dict[type, tuple[bool, bool]] = {
    StrictLessThan: (False, False),
    LessThan: (False, True),
    StrictGreaterThan: (True, False),
    GreaterThan: (True, True),
}


def guard(temporaries: Temporaries, exprs: Sequence[Basic]) -> Temporaries:
    """*temporaries*, each read only inside ``Piecewise`` branches computed only there.

    :func:`sympy.cse` lifts a sub-expression out of the branches using it into an
    unconditional temporary, which is then evaluated everywhere: a rate fit valid
    below 500 K overflows at 5000 K (Python raises, C traps when floating-point
    exceptions are enabled), and the work of the untaken branches is wasted.  Here
    every temporary that *exprs* read only under conditions, directly or through
    other temporaries, becomes ``Piecewise((value, guard), (0.0, True))``: the guard
    is the union of the conditions it is read under.  The ``0.0`` is never read.

    The conditions are those of ``Piecewise`` branches (each needing the earlier
    conditions false), of ``ITE`` and of the operands of ``And`` / ``Or`` after the
    first one printed (generated code short-circuits; printers sort the operands).
    Each is weakened to what it implies about the inputs: a comparison reading a
    temporary counts as true.  A guard therefore reads no temporary (a condition
    may test the very temporary it selects), and may hold in more cases than the
    exact one, never fewer.  Guards on one variable are simplified to intervals
    (``tgas >= 10 & tgas < 500``).

    Parameters
    ----------
    temporaries : list of (Symbol, Expr)
        CSE temporaries in evaluation order.
    exprs : Sequence[Basic]
        The expressions reading them.

    Returns
    -------
    list of (Symbol, Expr)
        The temporaries in the same order; the outputs of *exprs* are unchanged.

    Examples
    --------
    >>> import sympy as sp
    >>> from jaff.common.optimizations import guard
    >>> t, x, c = sp.symbols("tgas x cse0")
    >>> guard([(c, sp.exp(x))], [sp.Piecewise((1, t < 10), (c, t < 500), (2, True))])
    [(cse0, Piecewise((exp(x), (tgas >= 10) & (tgas < 500)), (0.0, True)))]
    """
    names = {symbol for symbol, _ in temporaries}
    paths = _paths(temporaries, exprs, names)
    guarded: Temporaries = []
    for symbol, value in temporaries:
        condition = _guard(paths[symbol], names)
        if condition is not true:
            value = cast(Expr, Piecewise((value, condition), (Float(0.0), True)))
        guarded.append((symbol, value))

    return guarded


def _paths(
    temporaries: Temporaries, exprs: Sequence[Basic], names: set[Basic]
) -> dict[Basic, set[Path]]:
    """Each temporary's paths: the condition sets it is read under."""
    paths: dict[Basic, set[Path]] = defaultdict(set)
    for expr in exprs:
        _walk(expr, frozenset(), names, paths)

    # A temporary is read only by later ones: last to first, its paths are complete.
    for symbol, value in reversed(temporaries):
        reads = paths[symbol]
        for path in [frozenset()] if frozenset() in reads else list(reads):
            _walk(value, path, names, paths)

    return paths


def _walk(
    expr: Basic, path: Path, names: set[Basic], paths: dict[Basic, set[Path]]
) -> None:
    """Record under *path* each temporary *expr* reads, extending it at branches."""
    if expr in names:
        paths[expr].add(path)
        return
    if not expr.args:
        return

    if isinstance(expr, Piecewise):
        earlier: list[Boolean] = []
        for pair in expr.args:
            value, condition = pair.args[0], cast(Boolean, pair.args[1])
            here = path.union(earlier)
            _walk(condition, here, names, paths)
            _walk(value, here if condition is true else here | {condition}, names, paths)
            earlier.append(Not(condition))
        return

    if isinstance(expr, ITE):
        test, if_true, if_false = expr.args
        test = cast(Boolean, test)
        _walk(test, path, names, paths)
        _walk(if_true, path | {test}, names, paths)
        _walk(if_false, path | {Not(test)}, names, paths)
        return

    if isinstance(expr, (And, Or)):
        # In printed order; each operand runs only if the earlier ones did not decide.
        earlier = []
        for operand in sorted(expr.args, key=default_sort_key):
            _walk(operand, path.union(earlier), names, paths)
            earlier.append(operand if isinstance(expr, And) else Not(operand))
        return

    for arg in expr.args:
        _walk(arg, path, names, paths)


def _guard(paths: set[Path], names: set[Basic]) -> Boolean:
    """The union of *paths*, weakened to the inputs; ``true`` if unconditional."""
    weakened = {frozenset(_weaken(c, names) for c in p) - {true} for p in paths}
    if not weakened or frozenset() in weakened:
        return true

    # (A) | (A & B) is A: drop paths containing another.
    kept = [p for p in weakened if not any(q < p for q in weakened)]

    return _simplified(Or(*(And(*p) for p in kept)))


def _weaken(condition: Boolean, names: set[Basic], negated: bool = False) -> Boolean:
    """A condition on the inputs alone implied by *condition* (by its negation if
    *negated*): each comparison reading a temporary becomes ``true``."""
    if not condition.free_symbols & names:
        return Not(condition) if negated else condition

    args = [cast(Boolean, a) for a in condition.args]
    if isinstance(condition, Not):
        return _weaken(args[0], names, not negated)
    if isinstance(condition, (And, Or)):
        parts = [_weaken(a, names, negated) for a in args]
        return And(*parts) if isinstance(condition, And) != negated else Or(*parts)
    if isinstance(condition, ITE):
        test, if_true, if_false = args
        return Or(
            And(_weaken(test, names), _weaken(if_true, names, negated)),
            And(_weaken(test, names, True), _weaken(if_false, names, negated)),
        )

    return true


def _simplified(condition: Boolean) -> Boolean:
    """*condition* as intervals when it bounds one variable by numbers."""
    if len(condition.free_symbols) != 1:
        return condition

    (variable,) = condition.free_symbols
    intervals = _intervals(condition, variable)
    if intervals is None:
        return condition

    return Or(*(_as_condition(i, variable) for i in intervals))


def _intervals(condition: Basic, variable: Basic) -> Intervals | None:
    """Where *condition* holds, or ``None`` if it is not built from comparisons of
    *variable* with numbers by ``And``, ``Or``, ``Not`` and ``ITE``."""
    if condition is true:
        return [(-math.inf, False, math.inf, False)]
    if condition is false:
        return []
    if isinstance(condition, Relational):
        return _comparison(condition, variable)
    if not isinstance(condition, (Not, And, Or, ITE)):
        return None

    parts: list[Intervals] = []
    for arg in condition.args:
        part = _intervals(arg, variable)
        if part is None:
            return None
        parts.append(part)

    return _combine(condition, parts)


def _combine(condition: Basic, parts: list[Intervals]) -> Intervals:
    """The intervals of the boolean *condition* from those of its arguments."""
    if isinstance(condition, Not):
        return _complement(parts[0])
    if isinstance(condition, ITE):
        test, if_true, if_false = parts
        return _union(
            _intersection(test, if_true) + _intersection(_complement(test), if_false)
        )
    if isinstance(condition, Or):
        return _union([i for p in parts for i in p])

    result = parts[0]
    for p in parts[1:]:
        result = _intersection(result, p)

    return result


def _comparison(condition: Relational, variable: Basic) -> Intervals | None:
    """The interval of ``variable < number`` (or ``<=``, ``>``, ``>=``, either side)."""
    if type(condition) not in _BOUNDS:
        return None

    lhs, rhs = condition.args
    if rhs == variable and lhs.is_Number:
        condition, lhs, rhs = condition.reversed, rhs, lhs
    if lhs != variable or not rhs.is_Number:
        return None

    lower, closed = _BOUNDS[type(condition)]
    bound = float(cast(Expr, rhs))
    if lower:
        return [(bound, closed, math.inf, False)]

    return [(-math.inf, False, bound, closed)]


def _union(intervals: Intervals) -> Intervals:
    """*intervals* sorted, empty ones dropped and overlapping or touching ones merged."""
    merged: Intervals = []
    for low, low_closed, high, high_closed in sorted(
        (i for i in intervals if i[0] < i[2] or (i[0] == i[2] and i[1] and i[3])),
        key=lambda i: (i[0], not i[1]),
    ):
        if merged and _meets(merged[-1], low, low_closed):
            last = merged[-1]
            if high > last[2] or (high == last[2] and high_closed):
                merged[-1] = (last[0], last[1], high, high_closed)
        else:
            merged.append((low, low_closed, high, high_closed))

    return merged


def _meets(interval: Interval, low: float, closed: bool) -> bool:
    """Whether an interval starting at *low* overlaps or touches *interval*."""
    return low < interval[2] or (low == interval[2] and (closed or interval[3]))


def _intersection(a: Intervals, b: Intervals) -> Intervals:
    """Where both *a* and *b* hold."""
    pieces: Intervals = []
    for a_low, a_low_closed, a_high, a_high_closed in a:
        for b_low, b_low_closed, b_high, b_high_closed in b:
            low = max(a_low, b_low)
            low_closed = (a_low_closed or a_low != low) and (b_low_closed or b_low != low)
            high = min(a_high, b_high)
            high_closed = (a_high_closed or a_high != high) and (
                b_high_closed or b_high != high
            )
            pieces.append((low, low_closed, high, high_closed))

    return _union(pieces)


def _complement(intervals: Intervals) -> Intervals:
    """Where *intervals* do not hold."""
    pieces: Intervals = []
    low, low_closed = -math.inf, False
    for start, start_closed, end, end_closed in intervals:
        pieces.append((low, low_closed, start, not start_closed))
        low, low_closed = end, not end_closed
    pieces.append((low, low_closed, math.inf, False))

    return _union(pieces)


def _as_condition(interval: Interval, variable: Basic) -> Boolean:
    """``low < variable < high`` (each side closed or open, infinite sides dropped)."""
    low, low_closed, high, high_closed = interval
    bounds: list[Boolean] = []
    if low != -math.inf:
        bounds.append(
            (GreaterThan if low_closed else StrictGreaterThan)(variable, _number(low))
        )
    if high != math.inf:
        bounds.append(
            (LessThan if high_closed else StrictLessThan)(variable, _number(high))
        )

    return And(*bounds)


def _number(value: float) -> Expr:
    """*value* as a SymPy number: an integer if it is one (``tgas < 500``)."""
    return Integer(int(value)) if value == int(value) else Float(value)
