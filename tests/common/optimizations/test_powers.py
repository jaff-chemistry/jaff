# ABOUTME: Tests for the power rules of optimizations: integer_power (small integer powers
# ABOUTME: of symbols as products) and exp_power (b**x as exp(ln(b)*x)), via optimize

import math

import pytest
import sympy as sp

from jaff.codegen import Language
from jaff.common.optimizations import exp_power, integer_power, optimize

X, Y = sp.symbols("x y")
NDEN = sp.IndexedBase("nden")
LN10 = sp.Float(math.log(10))


def _cxx(expr):
    return Language("c++").code_gen(expr, strict=False, allow_unknown_functions=True)


def _integer_powers_left(expr):
    """Powers with an integer exponent other than -1 left in *expr*."""
    return [
        p
        for p in sp.preorder_traversal(expr)
        if isinstance(p, sp.Pow) and p.exp.is_Integer and p.exp != -1
    ]


# integer_power


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
    assert _cxx(integer_power(expr)) == printed


@pytest.mark.parametrize("expr", [X, X**-1, X**9, X**-9, X**0.5, X**Y, (X + 1) ** 3])
def test_other_powers_are_unchanged(expr):
    assert integer_power(expr) == expr


def test_indexed_base():
    rewritten = integer_power(NDEN[3] ** 3)

    assert not _integer_powers_left(rewritten)
    assert rewritten.doit() == NDEN[3] ** 3


def test_limit():
    assert integer_power(X**5, limit=4) == X**5
    assert not _integer_powers_left(integer_power(X**5, limit=5))


def test_integer_powers_inside_expressions():
    expr = 3 * X**3 + sp.exp(Y * X**4) - Y / X**2 + sp.Piecewise((X**2, Y < 1), (0, True))
    rewritten = optimize(expr)

    assert not _integer_powers_left(rewritten)
    assert rewritten.doit() == expr


@pytest.mark.parametrize("lang", sorted(lang.name for lang in Language.registered()))
def test_every_language_prints_products(lang):
    expr = 3 * X**3 + Y / X**2 + sp.exp(Y * X**4)
    code = Language(lang).code_gen(
        optimize(expr), strict=False, allow_unknown_functions=True
    )

    assert "pow" not in code.replace(".powi", "pow")
    assert "**" not in code and "^" not in code


def test_rust_reciprocal_covers_the_whole_product():
    # `x*x.recip()` would be x * (1/x).
    code = Language("rust").code_gen(optimize(Y / X**2))

    assert code == "y*(x*x).recip()"


# exp_power


@pytest.mark.parametrize(
    ("expr", "want"),
    [
        (10**X, sp.exp(LN10 * X)),
        (sp.Float(10.0) ** X, sp.exp(LN10 * X)),
        (2 ** (X + 1), sp.exp(sp.Float(math.log(2)) * (X + 1))),
        (sp.Rational(1, 2) ** X, sp.exp(sp.Float(math.log(0.5)) * X)),
    ],
)
def test_numeric_base_power_is_an_exponential(expr, want):
    assert exp_power(expr) == want


@pytest.mark.parametrize("expr", [sp.Integer(10) ** 3, X**10, (-2) ** X, X**Y])
def test_other_bases_and_numeric_powers_are_unchanged(expr):
    assert exp_power(expr) == expr


def test_exponential_does_not_turn_back_into_a_power():
    assert "pow" not in _cxx(optimize(3 * 10 ** (Y - 2) + X))


# optimize


def test_rules_compose_in_one_pass():
    # x**3 inside the exponent is rewritten first (bottom-up), then 10**(...).
    assert _cxx(optimize(10 ** (X**3))) == "std::exp(2.3025850929940459*((x*x)*x))"


def test_each_distinct_subtree_is_visited_once():
    seen = []

    def count(node):
        seen.append(node)
        return node

    shared = sp.exp(X**3)
    optimize(shared * Y + shared**2 + sp.sin(shared), rules=(count,))

    assert len(seen) == len(set(seen))


@pytest.mark.parametrize(("x", "y"), [(1.7, 0.3), (-2.5, 4.0), (1e-3, 1e3)])
def test_python_value_is_unchanged(x, y):
    expr = 3 * X**3 + Y / X**2 + sp.exp(-Y * X**-4) - X**8 + 10 ** (-Y / 1000)
    code = Language("python").code_gen(optimize(expr))

    assert eval(code, {"math": math, "x": x, "y": y}) == pytest.approx(
        float(expr.subs({X: x, Y: y})), rel=1e-13
    )
