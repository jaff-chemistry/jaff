# ABOUTME: Tests for optimizations.log_derivative: value · ∂ln(expr)/∂y for products of
# ABOUTME: powers and exponentials when cheaper than diff, else None

import pytest
import sympy as sp

from jaff.common.optimizations import log_derivative

T = sp.Symbol("tgas")
K = sp.IndexedBase("krate")[0]
PHOTDEN = sp.IndexedBase("photden")

ARRHENIUS = 2e-10 * (T / 300) ** -0.7 * sp.exp(-50 / T)
FIT = 10 ** (-9.5 + 0.3 * sp.log(T, 10) - 0.1 * sp.log(T, 10) ** 2)


@pytest.mark.parametrize("expr", [ARRHENIUS, FIT], ids=["arrhenius", "fit"])
def test_applies_to_products_of_powers_and_exponentials(expr):
    result = log_derivative(expr, T, K)

    assert result is not None
    assert result.has(K)
    assert sp.count_ops(result) < sp.count_ops(sp.diff(expr, T))


@pytest.mark.parametrize("expr", [ARRHENIUS, FIT], ids=["arrhenius", "fit"])
@pytest.mark.parametrize("t", [30.0, 2000.0])
def test_is_exact(expr, t):
    got = float(log_derivative(expr, T, K).xreplace({K: expr}).xreplace({T: t}))

    assert got == pytest.approx(float(sp.diff(expr, T).xreplace({T: t})), rel=1e-12)


def test_arrhenius_drops_the_exponential():
    assert not log_derivative(ARRHENIUS, T, K).has(sp.exp)


def test_logarithm_of_a_sum_does_not_apply():
    assert log_derivative(1e-10 * T**0.5 + 3e-12 * sp.exp(-100 / T), T, K) is None


def test_linear_term_does_not_apply():
    # A photo rate is linear in the band densities: diff gives a bare number.
    photo = 1e-10 * PHOTDEN[1] + 2e-11 * PHOTDEN[2]

    assert log_derivative(photo, PHOTDEN[1], K) is None


def test_does_not_apply_when_diff_is_cheaper():
    # d(3e-10 T)/dT = 3e-10 beats K/T.
    assert log_derivative(3e-10 * T, T, K) is None


@pytest.mark.parametrize("expr", [sp.S.Zero, sp.Float(2e-10)], ids=["zero", "constant"])
def test_does_not_apply_without_a_dependence(expr):
    assert log_derivative(expr, T, K) is None
