# ABOUTME: derivative(): ∂expr/∂y per top-level Piecewise branch, written with a symbol
# ABOUTME: holding expr's value where the log-derivative optimization applies
"""Derivatives that may reuse the value of the differentiated expression."""

from __future__ import annotations

from typing import cast

from sympy import Basic, Expr, Piecewise, S, diff, sympify

from .optimizations import log_derivative


def derivative(expr: Basic, y: Basic, value: Basic | None = None) -> Expr:
    """``∂expr/∂y``, written with *value* where that is cheaper.

    Each branch of ``expr`` (the expression itself, or each branch of a top-level
    ``Piecewise``) is differentiated on its own: as ``value · ∂ln(b)/∂y``
    (:func:`~jaff.common.optimizations.log_derivative`) where that rewrite applies,
    else plainly.  The result is exact; a ``Piecewise`` free of *y* gives ``0``.

    Parameters
    ----------
    expr : Basic
        The expression to differentiate.
    y : Basic
        The variable, a ``Symbol`` or an ``Indexed`` such as ``nden[2]``.
    value : Basic, optional
        A symbol holding the value of *expr* (e.g. a rate coefficient's
        ``krate[r]``); without it the result is ``diff(expr, y)``.

    Returns
    -------
    Expr

    Examples
    --------
    >>> import sympy as sp
    >>> from jaff.common.derivative import derivative
    >>> tgas, k = sp.Symbol("tgas"), sp.Symbol("k")
    >>> rate = sp.Piecewise(
    ...     (2e-10 * sp.exp(-50 / tgas), tgas < 100), (1e-10 * tgas, True)
    ... )
    >>> derivative(rate, tgas, k)
    Piecewise((50*k/tgas**2, tgas < 100), (1.0e-10, True))
    """
    expr = cast(Expr, sympify(expr))
    if value is None:
        return diff(expr, y)
    if not isinstance(expr, Piecewise):
        return _branch(expr, y, value)

    pieces = [
        (_branch(cast(Expr, branch), y, value), cond)
        for branch, cond in cast(tuple[tuple[Basic, Basic], ...], expr.args)
    ]
    if all(piece == 0 for piece, _ in pieces):
        return S.Zero

    return Piecewise(*pieces)


def _branch(branch: Expr, y: Basic, value: Basic) -> Expr:
    """``∂branch/∂y``, by the log-derivative optimization when it applies."""
    rewritten = log_derivative(branch, y, value)

    return rewritten if rewritten is not None else diff(branch, y)
