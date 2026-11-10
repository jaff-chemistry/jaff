# ABOUTME: Tests for optimizations.derivative: ∂expr/∂y written as value · ∂ln(expr)/∂y
# ABOUTME: where cheaper, per Piecewise branch, else plain diff; always exact

import pytest
import sympy as sp

from jaff.common.optimizations import derivative

T = sp.Symbol("tgas")
K = sp.IndexedBase("krate")[0]
PHOTDEN = sp.IndexedBase("photden")

ARRHENIUS = 2e-10 * (T / 300) ** -0.7 * sp.exp(-50 / T)
FIT = 10 ** (-9.5 + 0.3 * sp.log(T, 10) - 0.1 * sp.log(T, 10) ** 2)
SUM = 1e-10 * T**0.5 + 3e-12 * sp.exp(-100 / T)
PIECEWISE = sp.Piecewise((ARRHENIUS, T < 100), (0, T < 200), (SUM, True))


def _value(expr, at):
    return float(expr.xreplace(at))


def _is_exact(result, expr, y, point):
    """*result* with K = *expr* equals diff(expr, y) at *point*."""
    got = result.xreplace({K: expr})

    return _value(got, point) == pytest.approx(_value(sp.diff(expr, y), point), rel=1e-12)


@pytest.mark.parametrize("expr", [ARRHENIUS, FIT, SUM], ids=["arrhenius", "fit", "sum"])
@pytest.mark.parametrize("t", [30.0, 150.0, 2000.0])
def test_is_exact(expr, t):
    assert _is_exact(derivative(expr, T, K), expr, T, {T: t})


@pytest.mark.parametrize("t", [30.0, 150.0, 2000.0])
def test_piecewise_is_exact_in_every_branch(t):
    assert _is_exact(derivative(PIECEWISE, T, K), PIECEWISE, T, {T: t})


def test_arrhenius_reuses_the_value():
    result = derivative(ARRHENIUS, T, K)

    assert result.has(K)
    assert not result.has(sp.exp)
    assert sp.count_ops(result) < sp.count_ops(sp.diff(ARRHENIUS, T))


def test_fit_reuses_the_value():
    result = derivative(FIT, T, K)

    assert result.has(K)
    assert sp.count_ops(result) < sp.count_ops(sp.diff(FIT, T))


def test_sum_keeps_the_plain_derivative():
    assert derivative(SUM, T, K) == sp.diff(SUM, T)


def test_piecewise_chooses_per_branch():
    result = derivative(PIECEWISE, T, K)
    (arrhenius, _), (zero, _), (total, _) = result.args

    assert isinstance(result, sp.Piecewise)
    assert arrhenius.has(K)
    assert zero == 0
    assert total == sp.diff(SUM, T)


def test_linear_term_keeps_the_plain_derivative():
    # A photo rate is linear in the band densities: diff gives a bare number.
    photo = 1e-10 * PHOTDEN[1] + 2e-11 * PHOTDEN[2]

    assert derivative(photo, PHOTDEN[1], K) == pytest.approx(1e-10)


def test_cheaper_plain_derivative_wins():
    # d(3e-10 T)/dT = 3e-10 beats K/T.
    assert derivative(3e-10 * T, T, K) == pytest.approx(3e-10)


def test_constant_gives_zero():
    assert derivative(sp.Float(2e-10), T, K) == 0


def test_without_a_value_is_plain_diff():
    assert derivative(ARRHENIUS, T) == sp.diff(ARRHENIUS, T)
