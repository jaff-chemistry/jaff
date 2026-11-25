# ABOUTME: Tests for common.derivative: ∂expr/∂y per Piecewise branch, using the
# ABOUTME: log-derivative optimization where it applies, always exact

import pytest
import sympy as sp

from jaff.common.derivative import derivative

T = sp.Symbol("tgas")
K = sp.IndexedBase("krate")[0]

ARRHENIUS = 2e-10 * (T / 300) ** -0.7 * sp.exp(-50 / T)
FIT = 10 ** (-9.5 + 0.3 * sp.log(T, 10) - 0.1 * sp.log(T, 10) ** 2)
SUM = 1e-10 * T**0.5 + 3e-12 * sp.exp(-100 / T)
PIECEWISE = sp.Piecewise((ARRHENIUS, T < 100), (0, T < 200), (SUM, True))


def _is_exact(result, expr, point):
    """*result* with K = *expr* equals diff(expr, T) at *point*."""
    got = float(result.xreplace({K: expr}).xreplace(point))
    want = float(sp.diff(expr, T).xreplace(point))

    return got == pytest.approx(want, rel=1e-12)


@pytest.mark.parametrize("expr", [ARRHENIUS, FIT, SUM], ids=["arrhenius", "fit", "sum"])
@pytest.mark.parametrize("t", [30.0, 150.0, 2000.0])
def test_is_exact(expr, t):
    assert _is_exact(derivative(expr, T, K), expr, {T: t})


@pytest.mark.parametrize("t", [30.0, 150.0, 2000.0])
def test_piecewise_is_exact_in_every_branch(t):
    assert _is_exact(derivative(PIECEWISE, T, K), PIECEWISE, {T: t})


def test_reuses_the_value_where_the_optimization_applies():
    result = derivative(ARRHENIUS, T, K)

    assert result.has(K)
    assert not result.has(sp.exp)


def test_piecewise_chooses_per_branch():
    result = derivative(PIECEWISE, T, K)
    (arrhenius, _), (zero, _), (total, _) = result.args

    assert isinstance(result, sp.Piecewise)
    assert arrhenius.has(K)
    assert zero == 0
    assert total == sp.diff(SUM, T)


def test_piecewise_free_of_the_variable_gives_zero():
    expr = sp.Piecewise((sp.Float(1e-10), T < 100), (sp.Float(2e-10), True))

    assert derivative(expr, T, K) == 0


def test_constant_gives_zero():
    assert derivative(sp.Float(2e-10), T, K) == 0


def test_without_a_value_is_plain_diff():
    assert derivative(ARRHENIUS, T) == sp.diff(ARRHENIUS, T)
    assert derivative(PIECEWISE, T) == sp.diff(PIECEWISE, T)
