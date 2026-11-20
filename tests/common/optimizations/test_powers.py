# ABOUTME: Tests for optimizations.integer_powers: small integer powers of symbols print
# ABOUTME: as products in every language, with the same value

import math

import pytest
import sympy as sp

from jaff.codegen import Language
from jaff.common.optimizations import integer_powers

X, Y = sp.symbols("x y")
NDEN = sp.IndexedBase("nden")


def _cxx(expr):
    return Language("c++").code_gen(expr, strict=False, allow_unknown_functions=True)


def _integer_powers_left(expr):
    """Powers with an integer exponent other than -1 left in *expr*."""
    return [
        p
        for p in sp.preorder_traversal(expr)
        if isinstance(p, sp.Pow) and p.exp.is_Integer and p.exp != -1
    ]


@pytest.mark.parametrize(
    ("expr", "printed"),
    [
        (X**2, "x*x"),
        (X**3, "(x*x)*x"),
        (X**4, "(x*x)*(x*x)"),
        (X**8, "((x*x)*(x*x))*((x*x)*(x*x))"),
        (X**-2, "1.0/(x*x)"),
        (X**-3, "1.0/((x*x)*x)"),
    ],
)
def test_power_prints_as_a_product(expr, printed):
    assert _cxx(integer_powers(expr)) == printed


@pytest.mark.parametrize("expr", [X, X**-1, X**9, X**-9, X**0.5, X**Y, (X + 1) ** 3])
def test_other_powers_are_unchanged(expr):
    assert integer_powers(expr) == expr


def test_indexed_base():
    rewritten = integer_powers(NDEN[3] ** 3)

    assert not _integer_powers_left(rewritten)
    assert rewritten.doit() == NDEN[3] ** 3


def test_limit():
    assert integer_powers(X**5, limit=4) == X**5
    assert not _integer_powers_left(integer_powers(X**5, limit=5))


def test_powers_inside_expressions():
    expr = 3 * X**3 + sp.exp(Y * X**4) - Y / X**2 + sp.Piecewise((X**2, Y < 1), (0, True))
    rewritten = integer_powers(expr)

    assert not _integer_powers_left(rewritten)
    assert rewritten.doit() == expr


@pytest.mark.parametrize("lang", sorted(lang.name for lang in Language.registered()))
def test_every_language_prints_products(lang):
    expr = 3 * X**3 + Y / X**2 + sp.exp(Y * X**4)
    code = Language(lang).code_gen(
        integer_powers(expr), strict=False, allow_unknown_functions=True
    )

    assert "pow" not in code.replace(".powi", "pow")
    assert "**" not in code and "^" not in code


def test_rust_reciprocal_covers_the_whole_product():
    # `x*x.recip()` would be x * (1/x).
    code = Language("rust").code_gen(integer_powers(Y / X**2))

    assert code == "y*(x*x).recip()"


@pytest.mark.parametrize(("x", "y"), [(1.7, 0.3), (-2.5, 4.0), (1e-3, 1e3)])
def test_python_value_is_unchanged(x, y):
    expr = 3 * X**3 + Y / X**2 + sp.exp(-Y * X**-4) - X**8
    code = Language("python").code_gen(integer_powers(expr))

    assert eval(code, {"math": math, "x": x, "y": y}) == pytest.approx(
        float(expr.subs({X: x, Y: y}))
    )
