# ABOUTME: log_derivative(): ∂expr/∂y as value · ∂ln(expr)/∂y for products of powers
# ABOUTME: and exponentials, when cheaper than diff and a symbol holds expr's value
"""The exp/log rewrite of a derivative: reuse the differentiated expression's value."""

from __future__ import annotations

from sympy import Basic, Expr, diff, expand_log, log, sympify

from .cost import cost


def log_derivative(expr: Basic, y: Basic, value: Basic) -> Expr | None:
    """``∂expr/∂y`` as ``value · ∂ln(expr)/∂y``, or ``None`` when that saves nothing.

    With *value* a symbol that holds ``expr``'s value (e.g. a rate coefficient's
    ``krate[r]``, computed once for the rhs), a product of powers and exponentials
    (Arrhenius forms, ``10**poly(log T)`` fits) has a logarithm that expands into a
    few simple terms, so ``value · ∂ln(expr)/∂y`` avoids evaluating the
    exponentials again.  The rewrite is exact: ``value · ∂ln f/∂y = ∂f/∂y``.

    It does not apply (``None``) when ``expr`` does not depend on *y*, when its
    logarithm keeps a sum (``ln(a + b)``), or when it costs more to evaluate
    (:func:`~jaff.common.optimizations.cost`) than ``diff(expr, y)``; a tie goes to
    the rewrite, since *value* is already computed.  So a power law ``a·T^b``
    becomes ``b·value/T`` rather than a new ``pow``.

    Parameters
    ----------
    expr : Basic
        The expression to differentiate (no top-level ``Piecewise``: rewrite each
        branch, see :func:`jaff.common.derivative.derivative`).
    y : Basic
        The variable, a ``Symbol`` or an ``Indexed`` such as ``nden[2]``.
    value : Basic
        A symbol holding the value of *expr*.

    Returns
    -------
    Expr or None

    Examples
    --------
    >>> import sympy as sp
    >>> from jaff.common.optimizations import log_derivative
    >>> tgas, k = sp.Symbol("tgas"), sp.Symbol("k")
    >>> log_derivative(2e-10 * tgas**-0.7 * sp.exp(-50 / tgas), tgas, k)
    k*(-0.7/tgas + 50/tgas**2)
    >>> log_derivative(1e-10 * tgas + 3e-12 * tgas**2, tgas, k) is None
    True
    """
    expr = sympify(expr)
    plain = diff(expr, y)
    if plain == 0:
        return None

    logarithm = expand_log(log(expr), force=True)
    # A logarithm of a sum is left when expr is not a product of powers and
    # exponentials; its quotient form saves nothing.
    if any(lg.args[0].is_Add for lg in logarithm.atoms(log)):
        return None

    candidate = value * diff(logarithm, y)

    return candidate if cost(candidate) <= cost(plain) else None
