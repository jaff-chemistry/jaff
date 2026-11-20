# ABOUTME: integer_powers(): small integer powers of symbols as unevaluated products, so
# ABOUTME: generated code multiplies instead of calling pow() (C/C++ only expand x**2)
"""Integer powers printed as products."""

from __future__ import annotations

from sympy import Basic, Expr, Indexed, Mul, Pow, Symbol, UnevaluatedExpr

LIMIT = 8
"""Largest exponent magnitude turned into a product by default."""


def integer_powers(expr: Basic, limit: int = LIMIT) -> Basic:
    """*expr* with each ``b**n`` (``b`` a symbol, ``2 <= |n| <= limit``) as a product.

    C and C++ print ``x**3`` as ``pow(x, 3)``, a library call: compilers expand only
    ``pow(x, 2)`` unless ``-ffast-math`` is on.  Here the power becomes a product
    built by squaring (``x**8`` is ``((x*x)*(x*x))*((x*x)*(x*x))``, three
    multiplications once the compiler merges the identical ``x*x``) and a negative
    one its reciprocal, ``1.0/(x*x)``.  Products are wrapped in
    :class:`~sympy.core.expr.UnevaluatedExpr`, which SymPy would otherwise turn
    back into powers, so apply this last, right before printing: any later
    rewrite may undo it.  Every language prints the result as plain arithmetic.

    The value is unchanged (``doit()`` gives *expr* back); rounding may differ
    from ``pow`` by a few units in the last place.  Powers of compound bases are
    left alone, as printing the base several times would make the code unreadable.

    Parameters
    ----------
    expr : Basic
        The expression.
    limit : int, optional
        Largest exponent magnitude to expand.

    Returns
    -------
    Basic

    Examples
    --------
    >>> import sympy as sp
    >>> from jaff.common.optimizations import integer_powers
    >>> x, y = sp.symbols("x y")
    >>> sp.cxxcode(integer_powers(3 * x**3 + y / x**2))
    'y*1.0/(x*x) + 3*((x*x)*x)'
    """
    return expr.replace(
        lambda e: _expandable(e, limit),
        lambda e: _product(e.base, int(e.exp)),
    )


def _expandable(expr: Basic, limit: int) -> bool:
    """Whether *expr* is a power :func:`integer_powers` rewrites."""
    if not isinstance(expr, Pow) or not isinstance(expr.base, (Symbol, Indexed)):
        return False

    exp = expr.exp

    return bool(exp.is_Integer and 2 <= abs(int(exp)) <= limit)


def _product(base: Expr, n: int) -> Expr:
    """``base**n`` as an unevaluated product (its reciprocal if *n* < 0)."""
    product = _squares(base, abs(n))
    if n < 0:
        # The reciprocal inside the wrapper: Rust prints (x*x).recip(), not x*x.recip().
        return UnevaluatedExpr(Pow(product, -1, evaluate=False))

    return UnevaluatedExpr(product)


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
