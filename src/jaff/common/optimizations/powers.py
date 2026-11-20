# ABOUTME: Power rules for optimize(): integer_power (x**3 as an unevaluated (x*x)*x) and
# ABOUTME: exp_power (10**x as exp(2.302585...*x)); both avoid a pow() library call
"""Rules rewriting powers so generated code avoids ``pow()``."""

from __future__ import annotations

import math

from sympy import Basic, Expr, Float, Indexed, Mul, Pow, Symbol, UnevaluatedExpr, exp

LIMIT = 8
"""Largest exponent magnitude :func:`integer_power` turns into a product by default."""


def integer_power(expr: Basic, limit: int = LIMIT) -> Basic:
    """``b**n`` (``b`` a symbol, ``2 <= |n| <= limit``) as a product; else *expr*.

    C and C++ print ``x**3`` as ``pow(x, 3)``, a library call: compilers expand only
    ``pow(x, 2)`` unless ``-ffast-math`` is on.  Here the power becomes a product
    built by squaring (``x**8`` is ``((x*x)*(x*x))*((x*x)*(x*x))``, three
    multiplications once the compiler merges the identical ``x*x``) and a negative
    one its reciprocal, ``1.0/(x*x)``.  Products are wrapped in
    :class:`~sympy.core.expr.UnevaluatedExpr`, which SymPy would otherwise turn
    back into powers, so this runs last, right before printing (in
    :func:`~jaff.common.optimizations.optimize`): any later rewrite may undo it.
    Every language prints the result as plain arithmetic.

    The value is unchanged (``doit()`` gives *expr* back); rounding may differ
    from ``pow`` by a few units in the last place.  Powers of compound bases are
    left alone, as printing the base several times would make the code unreadable.

    Parameters
    ----------
    expr : Basic
        One node of an expression.
    limit : int, optional
        Largest exponent magnitude to expand.

    Returns
    -------
    Basic

    Examples
    --------
    >>> import sympy as sp
    >>> from jaff.common.optimizations import integer_power
    >>> x = sp.Symbol("x")
    >>> sp.cxxcode(integer_power(x**3)), sp.cxxcode(integer_power(x**-2))
    ('(x*x)*x', '1.0/(x*x)')
    """
    if not isinstance(expr, Pow) or not isinstance(expr.base, (Symbol, Indexed)):
        return expr

    n = expr.exp
    if not n.is_Integer or not 2 <= abs(int(n)) <= limit:
        return expr

    product = _squares(expr.base, abs(int(n)))
    if n.is_negative:
        # The reciprocal inside the wrapper: Rust prints (x*x).recip(), not x*x.recip().
        return UnevaluatedExpr(Pow(product, -1, evaluate=False))

    return UnevaluatedExpr(product)


def exp_power(expr: Basic) -> Basic:
    """``b**x`` (``b`` a positive number other than 1, ``x`` not a number) as
    ``exp(ln(b)*x)``; else *expr*.

    ``pow(10, x)`` costs about 2.4 times ``exp`` (glibc; it takes a slow path for
    tiny results), and the cooling and rate fits of networks such as GOW are
    ``10**polynomial``.  ``ln(b)`` is a float, so SymPy does not fold the result
    back into a power.  The relative error grows with the exponent: about
    ``1e-14`` for ``|x·ln b| <= 70`` (fits), ``1e-13`` near the overflow limit.

    Parameters
    ----------
    expr : Basic
        One node of an expression.

    Returns
    -------
    Basic

    Examples
    --------
    >>> import sympy as sp
    >>> from jaff.common.optimizations import exp_power
    >>> exp_power(10 ** sp.Symbol("x"))
    exp(2.30258509299405*x)
    """
    if not isinstance(expr, Pow) or expr.exp.is_Number:
        return expr

    base = expr.base
    if not base.is_Number or not base.is_positive or base == 1:
        return expr

    return exp(Float(math.log(float(base))) * expr.exp)


def _squares(base: Expr, n: int) -> Expr:
    """``base**n`` (``n >= 1``) as products by repeated squaring."""
    if n == 1:
        return base

    half = _squares(base, n // 2)
    square = Mul(_wrapped(half), _wrapped(half), evaluate=False)
    if n % 2 == 0:
        return square

    return Mul(UnevaluatedExpr(square), base, evaluate=False)


def _wrapped(expr: Expr) -> Expr:
    """*expr*, wrapped if it is a product, so it prints in parentheses."""
    return UnevaluatedExpr(expr) if isinstance(expr, Mul) else expr
