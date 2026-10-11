"""
Numerical integration helpers for SymPy expressions and tabulated data.

This module provides:

* :func:`integrate` -- fixed-sample quadrature (trapezoid or Simpson) via
  :func:`sympy.lambdify`, for use when a quick approximate answer is needed.
* :func:`get_bounds` -- extract breakpoints (discontinuities / kinks) from a
  piecewise SymPy expression so that adaptive integrators can avoid them.
* :func:`sym_integrate` -- adaptive quadrature via :func:`scipy.integrate.quad`
  that handles piecewise expressions, symbolic infinity bounds, and
  multi-decade ranges.
* :func:`arr_integrate` -- trapezoidal integral of tabulated ``y(x)`` arrays.
* :func:`smart_integrate` -- dispatches to :func:`arr_integrate` for arrays and
  to :func:`sym_integrate` for SymPy expressions.
"""

import operator
from collections.abc import Callable
from functools import lru_cache

import numpy as np
from scipy.integrate import quad, simpson, trapezoid
from sympy import (
    And,
    Basic,
    Expr,
    FiniteSet,
    Not,
    Or,
    Piecewise,
    S,
    Symbol,
    lambdify,
    piecewise_fold,
    solve,
)
from sympy.core.relational import (
    Equality,
    GreaterThan,
    LessThan,
    Relational,
    StrictGreaterThan,
    StrictLessThan,
    Unequality,
)

# Maximum number of distinct ``(expr, sym)`` preparations kept by ``_prepare``.
_PREPARE_CACHE_SIZE = 1024


def integrate(
    expr: Basic | Expr,
    sym: Basic,
    bounds: tuple[float | int, float | int],
    integrator: str = "trapezoid",
    spacing: str = "lin",
):
    """
    Numerically integrate a SymPy expression using fixed-sample quadrature.

    The expression is converted to a NumPy function via
    :func:`sympy.lambdify` and then evaluated on a uniform or logarithmically
    spaced grid of 1 000 000 points before applying the chosen quadrature
    rule.

    Parameters
    ----------
    expr : sympy.Basic or sympy.Expr
        The expression to integrate.
    sym : sympy.Basic
        The integration variable.
    bounds : tuple of (float or int, float or int)
        Integration interval ``(lower, upper)``.
    integrator : {"trapezoid", "simpson"}, optional
        Quadrature rule to apply (default ``"trapezoid"``).
    spacing : {"lin", "log"}, optional
        Sample-point spacing: ``"lin"`` for linear (default) or ``"log"`` for
        logarithmic.

    Returns
    -------
    float
        The approximate definite integral.

    Raises
    ------
    ValueError
        If *integrator* or *spacing* is not one of the supported options.
    """
    integrators = {"trapezoid": trapezoid, "simpson": simpson}
    spacings = {"lin": np.linspace, "log": np.logspace}

    if integrator not in integrators:
        raise ValueError(
            f"Invalid integrator specified: {integrator}\n"
            f"Supported integrators are {integrators.keys()}"
        )
    if spacing not in spacings:
        raise ValueError(
            f"Invalid spacing specified: {spacing}\n"
            f"Supported spacings are {spacings.keys()}"
        )

    lower, upper = bounds
    samples = spacings[spacing](lower, upper, 1_000_000)
    func = lambdify(sym, expr, "numpy")

    return integrators[integrator](func(samples), samples)


def get_bounds(expr: Basic, sym: Basic):
    """
    Extract breakpoints of a piecewise SymPy expression as floats.

    When every condition of the folded :class:`~sympy.Piecewise` is a simple
    comparison of *sym* with an integer or float (optionally combined with
    ``And`` / ``Or`` / ``Not``), the active set of each piece is decided by
    evaluating the conditions on the cells between the thresholds (see
    :func:`_get_bounds_fast`).  This yields exactly the same breakpoints as
    the general algorithm :func:`_get_bounds_solve` while skipping SymPy's
    expensive inequality solver, to which any other expression falls back.

    Parameters
    ----------
    expr : sympy.Basic
        A piecewise (or potentially piecewise) SymPy expression.
    sym : sympy.Basic
        The variable with respect to which breakpoints are identified.

    Returns
    -------
    list of float
        Sorted list of distinct breakpoint values.  Empty if no breakpoints
        were found.
    """
    fast = _get_bounds_fast(piecewise_fold(expr), sym)
    if fast is not None:
        return fast

    return _get_bounds_solve(expr, sym)


# Comparison operators for the simple relationals accepted by the fast path.
_REL_OPS = {
    StrictLessThan: operator.lt,
    LessThan: operator.le,
    StrictGreaterThan: operator.gt,
    GreaterThan: operator.ge,
    Equality: operator.eq,
    Unequality: operator.ne,
}


def _compile_condition(cond: Basic, sym: Basic, thresholds: set) -> Callable | None:
    """
    Compile a simple piecewise condition on *sym* into a float predicate.

    Supported conditions are ``True`` / ``False``, relationals
    (``<``, ``<=``, ``>``, ``>=``, ``Eq``, ``Ne``) between *sym* and an
    :class:`~sympy.Integer` or :class:`~sympy.Float` on either side, and
    ``And`` / ``Or`` / ``Not`` combinations of those.  Integers and
    floats convert to Python floats without ambiguity, so comparing in
    floating point decides each condition exactly as SymPy would.

    Parameters
    ----------
    cond : sympy.Basic
        A :class:`~sympy.Piecewise` condition.
    sym : sympy.Basic
        The integration variable.
    thresholds : set
        Receives the float threshold of every relational in *cond*.

    Returns
    -------
    callable or None
        ``pred(x) -> bool`` telling whether *cond* holds at ``sym = x``, or
        ``None`` if *cond* is not of the supported simple form.
    """
    if cond is S.true or cond is S.false:
        value = cond is S.true
        return lambda x: value

    if isinstance(cond, Relational):
        if cond.rhs == sym and cond.lhs != sym:
            cond = cond.reversed
        num = cond.rhs
        op = _REL_OPS.get(type(cond))
        if op is None or cond.lhs != sym or not (num.is_Integer or num.is_Float):
            return None
        t = float(num)
        if not np.isfinite(t):
            return None
        thresholds.add(t)
        return lambda x: op(x, t)

    if isinstance(cond, (And, Or, Not)):
        preds = [_compile_condition(arg, sym, thresholds) for arg in cond.args]
        if any(p is None for p in preds):
            return None
        if isinstance(cond, Not):
            return lambda x: not preds[0](x)
        combine = all if isinstance(cond, And) else any
        return lambda x: combine(p(x) for p in preds)

    return None


def _get_bounds_fast(folded: Basic, sym: Basic) -> list[float] | None:
    """
    Breakpoints of a folded piecewise with simple conditions on *sym*.

    Mirrors strategy 1 of :func:`_get_bounds_solve`
    (:meth:`sympy.Piecewise.as_expr_set_pairs`): pieces are walked in order,
    each piece is active on its condition set minus the sets of earlier
    conditions, and the finite boundary points of the active sets of
    non-zero pieces are collected.

    Instead of SymPy set arithmetic, the real line is split at the sorted
    condition thresholds ``t_1 < ... < t_k`` into the cells
    ``(-inf, t_1), {t_1}, (t_1, t_2), ..., {t_k}, (t_k, inf)``.  Every
    condition is constant on each cell, so it is evaluated once at a
    representative point per cell.  A threshold ``t_i`` lies on the boundary
    of an active set exactly when membership of ``{t_i}`` and its two
    neighbouring open cells is not all equal.

    Parameters
    ----------
    folded : sympy.Basic
        An expression already passed through :func:`sympy.piecewise_fold`.
    sym : sympy.Basic
        The integration variable.

    Returns
    -------
    list of float or None
        Sorted list of distinct breakpoints, or ``None`` when the fast path
        does not apply (not a piecewise, a non-simple condition, or no
        breakpoints found so that strategy 2 of :func:`_get_bounds_solve`
        would be needed).
    """
    if not isinstance(folded, Piecewise):
        return None

    thresholds = set()
    preds = []
    for _, cond in folded.args:
        pred = _compile_condition(cond, sym, thresholds)
        if pred is None:
            return None
        preds.append(pred)

    # Representative point of every cell: open cells at even, thresholds at odd
    ts = sorted(thresholds)
    reps = [ts[0] - 1.0] if ts else [0.0]
    for i, t in enumerate(ts):
        upper = ts[i + 1] if i + 1 < len(ts) else t + 2.0
        reps.extend([t, 0.5 * (t + upper)])

    boundaries = set()
    remaining = [True] * len(reps)
    for (val, _), pred in zip(folded.args, preds):
        active = [r and pred(x) for r, x in zip(remaining, reps)]
        remaining = [r and not a for r, a in zip(remaining, active)]
        if val == 0 or not any(active):
            continue

        for i, t in enumerate(ts):
            left, point, right = active[2 * i : 2 * i + 3]
            if not left == point == right:
                boundaries.add(t)

    if not boundaries:
        return None

    return sorted(list(boundaries))


def _get_bounds_solve(expr: Basic, sym: Basic) -> list[float]:
    """
    Extract breakpoints of a piecewise SymPy expression using SymPy's solvers.

    This is the general (slow) algorithm behind :func:`get_bounds`, used when
    its interval fast path does not apply.  Uses two strategies in order:

    1. **Domain boundaries** -- calls :func:`sympy.piecewise_fold` then
       iterates over ``(value, domain)`` pairs; for each non-zero value, the
       domain's boundary :class:`~sympy.sets.sets.FiniteSet` is collected.
    2. **Relational atoms** -- if no boundaries were found by strategy 1,
       solves each :class:`~sympy.core.relational.Relational` atom for *sym*
       to recover implicit breakpoints.

    Parameters
    ----------
    expr : sympy.Basic
        A piecewise (or potentially piecewise) SymPy expression.
    sym : sympy.Basic
        The variable with respect to which breakpoints are identified.

    Returns
    -------
    list of float
        Sorted list of distinct breakpoint values.  Empty if no breakpoints
        were found.
    """
    folded = piecewise_fold(expr)
    boundaries = set()

    # Strategy 1: read breakpoints directly from domain boundary sets
    if hasattr(folded, "as_expr_set_pairs"):
        for val, domain in folded.as_expr_set_pairs():
            if val == 0:
                continue

            b = domain.boundary
            if not isinstance(b, FiniteSet):
                continue

            for pt in b:
                boundaries.add(float(pt))

    # Strategy 2: solve relational conditions for sym when strategy 1 found nothing
    if not boundaries:
        for rel in folded.atoms(Relational):
            try:
                cp = solve(rel.lhs - rel.rhs, sym)
                for p in cp:
                    if p.is_real:
                        boundaries.add(float(p))
            except Exception:
                continue

    return sorted(list(boundaries))


def sym_integrate(
    expr: Basic, sym: Basic, bounds: tuple[float | int | Basic, float | int | Basic]
) -> float:
    """
    Adaptively integrate a (piecewise) SymPy expression over a possibly symbolic interval.

    Handles several complications that trip up naive quadrature:

    * **Piecewise expressions** -- :func:`sympy.piecewise_fold` is applied and
      breakpoints are extracted via :func:`get_bounds`; they are then passed to
      :func:`scipy.integrate.quad` as the ``points`` argument so the integrator
      can straddle discontinuities.
    * **Symbolic bounds** -- if *lower* or *upper* are SymPy
      :class:`~sympy.core.basic.Basic` objects (e.g. ``-oo`` / ``oo``)
      they are mapped to ``-np.inf`` / ``np.inf``.
    * **Multi-decade ranges** -- when the effective integration domain spans
      more than 0.1 decades in log10 space, logarithmically spaced interior
      sample hints are added to ``points`` to help the adaptive routine resolve
      steep rate curves.

    Parameters
    ----------
    expr : sympy.Basic
        The expression to integrate.
    sym : sympy.Basic
        The integration variable.
    bounds : tuple of (float or int or sympy.Basic, float or int or sympy.Basic)
        Integration interval ``(lower, upper)``.  Symbolic values are
        interpreted as ``±inf``.

    Returns
    -------
    float
        The approximate definite integral.  Returns ``0.0`` immediately if
        the effective integration interval is empty (``a >= b``).

    Notes
    -----
    When the interval contains infinities, :func:`scipy.integrate.quad` is
    called with ``limit=10000`` and no interior ``points`` hint (because
    ``quad`` does not accept ``points`` with infinite bounds).  For finite
    intervals ``limit=200`` is used.
    """
    lower, upper = bounds

    f_num, pts = _prepare(expr, sym)

    # Resolve symbolic bounds to float ±inf
    t_low = float(lower) if not isinstance(lower, Basic) else -np.inf
    t_high = float(upper) if not isinstance(upper, Basic) else np.inf

    # Clip the integration domain to the range covered by piecewise breakpoints
    a, b = t_low, t_high
    if pts:
        a = max(t_low, min(pts))
        b = min(t_high, max(pts))

    if a >= b:
        return 0.0

    sub_points = []

    # Add log-spaced interior hints when the domain spans >0.1 decades
    if a > 0 and not np.isinf(b) and b > a:
        decades = np.log10(b) - np.log10(a)
        if decades > 0.1:
            log_pts = np.logspace(np.log10(a), np.log10(b), int(decades * 5) + 2)
            sub_points.extend(log_pts[1:-1])  # Exclude the boundaries a and b

    # Add piecewise breakpoints that fall strictly inside [a, b]
    internal_bounds = [p for p in pts if a < p < b]
    sub_points.extend(internal_bounds)

    sub_points = sorted(list(set(sub_points)))

    # quad does not accept the ``points`` argument when bounds are infinite
    if np.isinf(a) or np.isinf(b):
        val, _ = quad(f_num, a, b, limit=10000)

        return val

    val, _ = quad(f_num, a, b, points=sub_points, limit=200)

    return val


@lru_cache(maxsize=_PREPARE_CACHE_SIZE)
def _prepare(expr: Basic, sym: Basic) -> tuple:
    """
    Fold, lambdify and extract breakpoints of *expr* once per ``(expr, sym)``.

    :func:`sym_integrate` is typically called many times with the same
    expression and different bounds (e.g. one call per radiation band), and
    the breakpoint extraction dominates its cost, so the preparation is
    memoised.

    Parameters
    ----------
    expr : sympy.Basic
        The expression to integrate.
    sym : sympy.Basic
        The integration variable.

    Returns
    -------
    tuple
        ``(f_num, pts)``: the NumPy-lambdified folded expression and the
        tuple of breakpoints from :func:`get_bounds`.
    """
    folded = piecewise_fold(expr)
    f_num = lambdify(sym, folded, "numpy")
    pts = tuple(get_bounds(folded, sym))

    return f_num, pts


def arr_integrate(
    y: np.ndarray, x: np.ndarray, bounds: tuple[float | int | Basic, float | int | Basic]
) -> float:
    """
    Trapezoidal integral of tabulated data ``y(x)`` over an interval.

    Used for cross-section integrals where ``σ(E)`` is only known as sampled
    ``(E, σ)`` arrays, so a closed form is unavailable.

    Parameters
    ----------
    y : numpy.ndarray
        Sampled integrand values, aligned with ``x``.
    x : numpy.ndarray
        Sample abscissae, assumed sorted ascending.
    bounds : tuple
        ``(lower, upper)`` integration limits.  A non-symbolic bound is used
        as-is; a symbolic (:class:`sympy.Basic`) bound is treated as ``±inf``,
        i.e. the open end of the tabulated range.  Both limits are clamped to
        ``[x[0], x[-1]]``.

    Returns
    -------
    float
        The integral, or ``0.0`` if the clamped interval is empty.

    Notes
    -----
    The endpoints are inserted into the sample grid via linear interpolation
    so the integral covers exactly ``[lower, upper]`` rather than the nearest
    sample points.
    """
    # Assumes data is sorted (ascending in x).
    lower, upper = bounds
    t_low = float(lower) if not isinstance(lower, Basic) else -np.inf
    t_high = float(upper) if not isinstance(upper, Basic) else np.inf

    t_low = max(t_low, x[0])
    t_high = min(t_high, x[-1])
    if t_high <= t_low:
        return 0.0

    i_low = np.searchsorted(x, t_low)
    i_high = np.searchsorted(x, t_high)

    x_seg = x[i_low:i_high]
    y_seg = y[i_low:i_high]

    x_seg = np.r_[t_low, x_seg, t_high]
    y_seg = np.r_[np.interp(t_low, x, y), y_seg, np.interp(t_high, x, y)]

    return np.trapezoid(y_seg, x_seg)


def smart_integrate(
    y: np.ndarray | Basic,
    x: np.ndarray | Basic,
    bounds: tuple[float | int | Basic, float | int | Basic],
) -> float:
    """
    Integrate tabulated data or a SymPy expression over an interval.

    Dispatches on the type of *y*: a :class:`numpy.ndarray` integrand is
    integrated with :func:`arr_integrate` over the sample abscissae *x*;
    anything else is treated as a SymPy expression in the symbol *x* and
    integrated with :func:`sym_integrate`.

    Parameters
    ----------
    y : numpy.ndarray or sympy.Basic
        Sampled integrand values, or a SymPy expression.
    x : numpy.ndarray or sympy.Basic
        Sample abscissae (for array *y*) or the integration symbol.
    bounds : tuple
        ``(lower, upper)`` integration limits; symbolic values mean ``±inf``.

    Returns
    -------
    float
        The definite integral.

    Raises
    ------
    TypeError
        If *y* is an array but *x* is not, or *y* is symbolic but *x* is not a
        :class:`sympy.Symbol`.
    """
    if isinstance(y, np.ndarray):
        if not isinstance(x, np.ndarray):
            raise TypeError(
                f"Array integrand needs ndarray abscissae, got {type(x).__name__}"
            )

        return arr_integrate(y, x, bounds)

    if not isinstance(x, Symbol):
        raise TypeError(
            f"Symbolic integrand needs a sympy Symbol variable, got {type(x).__name__}"
        )

    return sym_integrate(y, x, bounds)


def safe_integrate():
    """
    Placeholder for a future robust integration routine.

    Raises
    ------
    NotImplementedError
        Always.  This function is not yet implemented.
    """
    raise NotImplementedError("Not yet implemented")
