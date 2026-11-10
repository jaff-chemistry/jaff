# ABOUTME: derivative(): ∂expr/∂y, as value · ∂ln(expr)/∂y per Piecewise branch when
# ABOUTME: that is cheaper than plain diff and a symbol already holds expr's value
"""Derivatives that reuse the value of the differentiated expression."""

from __future__ import annotations

from typing import cast

from sympy import Basic, Expr, Piecewise, S, count_ops, diff, expand_log, log, sympify


def derivative(expr: Basic, y: Basic, value: Basic | None = None) -> Expr:
    """``∂expr/∂y``, written with *value* where that is cheaper.

    With *value* a symbol that holds ``expr``'s value (e.g. a rate coefficient's
    ``krate[r]``, computed once for the rhs), each branch ``b`` of ``expr`` (the
    expression itself, or each branch of a top-level ``Piecewise``) becomes
    ``value · ∂ln(b)/∂y`` when that has fewer operations (``count_ops``) than
    ``∂b/∂y``: products of powers and exponentials (Arrhenius forms,
    ``10**poly(log T)`` fits) reduce to a few simple terms.  A branch whose
    logarithm keeps a sum (``ln(a + b)``) or that is zero keeps its plain
    derivative.  The result is exact: ``value · ∂ln b/∂y = ∂b/∂y``.

    Parameters
    ----------
    expr : Basic
        The expression to differentiate.
    y : Basic
        The variable, a ``Symbol`` or an ``Indexed`` such as ``nden[2]``.
    value : Basic, optional
        A symbol holding the value of *expr*; without it the result is
        ``diff(expr, y)``.

    Returns
    -------
    Expr

    Examples
    --------
    >>> import sympy as sp
    >>> from jaff.common.optimizations import derivative
    >>> tgas, k = sp.Symbol("tgas"), sp.Symbol("k")
    >>> derivative(2e-10 * tgas**-0.7 * sp.exp(-50 / tgas), tgas, k)
    k*(-0.7/tgas + 50/tgas**2)
    >>> derivative(1e-10 * tgas + 3e-12 * tgas**2, tgas, k)
    6.0e-12*tgas + 1.0e-10
    """
    expr = cast(Expr, sympify(expr))
    if value is None:
        return diff(expr, y)
    if isinstance(expr, Piecewise):
        pieces = [
            (_branch(cast(Expr, branch), y, value), cond)
            for branch, cond in cast(tuple[tuple[Basic, Basic], ...], expr.args)
        ]
        if all(piece == 0 for piece, _ in pieces):
            return S.Zero

        return Piecewise(*pieces)

    return _branch(expr, y, value)


def _branch(branch: Expr, y: Basic, value: Basic) -> Expr:
    """``∂branch/∂y``, as ``value · ∂ln(branch)/∂y`` when that is cheaper."""
    plain = diff(branch, y)
    if plain == 0:
        return plain

    logarithm = expand_log(log(branch), force=True)
    # A logarithm of a sum is left when the branch is not a product of powers and
    # exponentials; its quotient form saves nothing.
    if any(lg.args[0].is_Add for lg in logarithm.atoms(log)):
        return plain

    candidate = value * diff(logarithm, y)

    return candidate if count_ops(candidate) < count_ops(plain) else plain
