# ABOUTME: Tests for optimizations.split: split_log (log(c*exp(b)*R) as log(c) + b + log(R))
# ABOUTME: and split_power ((c*exp(b)*R)**e as c**e*exp(e*b)*R**e), with 10**x like exp

import math

import pytest
import sympy as sp

from jaff.common.optimizations import PRE_CSE_RULES, bottom_up, split_log, split_power

T, X, Y = sp.symbols("tgas x y")


def test_constant_times_exponential():
    rewritten = split_log(sp.log(3.52e-9 * sp.exp(-43900 / T)))

    assert rewritten == sp.Float(math.log(3.52e-9)) - 43900 / T


def test_log_of_exponential():
    assert split_log(sp.log(sp.exp(X))) == X


def test_several_exponentials_and_a_rest():
    rewritten = split_log(sp.log(2 * sp.exp(X) * sp.exp(Y) * sp.sqrt(T)))

    assert rewritten == sp.Float(math.log(2)) + X + Y + sp.log(sp.sqrt(T))


@pytest.mark.parametrize(
    "expr",
    [
        sp.log(3 * T),  # no exponential: nothing underflows
        sp.log(10 * T),  # a numeric base without a symbolic exponent is a constant
        sp.log(-2 * sp.exp(X)),  # negative constant: log(-2) would be complex
        sp.log(T),
        sp.exp(X) * 3,
    ],
)
def test_unchanged(expr):
    assert split_log(expr) == expr


def test_value_is_unchanged():
    expr = sp.log(6.67e-12 * sp.sqrt(T) * sp.exp(-1 - 63590 / T))
    rewritten = split_log(expr)

    for t in (500.0, 3000.0):
        assert float(rewritten.subs(T, t)) == pytest.approx(float(expr.subs(T, t)))


def test_printed_python_is_finite_where_the_exponential_underflows():
    expr = sp.log(3.52e-9 * sp.exp(-43900 / T))
    namespace = {"math": math, "tgas": 20.0}

    with pytest.raises(ValueError):  # exp(-2195) underflows to 0: math.log(0)
        eval(sp.pycode(expr), namespace)
    assert eval(sp.pycode(split_log(expr)), namespace) == pytest.approx(
        math.log(3.52e-9) - 43900 / 20.0
    )


def test_pre_cse_rules_reach_nested_logs():
    expr = 1e-10 / (1 - 0.01 * sp.log(3.52e-9 * sp.exp(-43900 / T))) + X
    rewritten = bottom_up(PRE_CSE_RULES)(expr)

    assert not any(isinstance(a, sp.exp) for a in sp.preorder_traversal(rewritten))


LN10 = sp.Float(math.log(10))


def test_log_of_a_numeric_base_power():
    assert split_log(sp.log(10 ** (X - 3))) == (X - 3) * LN10


def test_log_with_a_numeric_base_power_factor():
    rewritten = split_log(sp.log(2 * 10**X * T))

    assert rewritten == sp.Float(math.log(2)) + X * LN10 + sp.log(T)


# split_power


def test_power_of_constant_times_exponential():
    rewritten = split_power((3.52e-9 * sp.exp(-43900 / T)) ** Y)

    assert rewritten == sp.Float(3.52e-9) ** Y * sp.exp(-43900 * Y / T)


def test_power_of_exponential():
    assert split_power(sp.exp(X) ** 0.5) == sp.exp(0.5 * X)


def test_power_with_a_numeric_base_power_factor():
    rewritten = split_power((10**X * T) ** 0.5)

    assert rewritten == sp.exp(0.5 * X * LN10) * T**0.5


@pytest.mark.parametrize(
    "expr",
    [
        (T + 1) ** 0.5,  # not a product
        (3 * T) ** 0.5,  # no exponential: nothing to split off
        (2 * sp.exp(X)) ** 2,  # integer exponent
        (-2 * sp.exp(X)) ** Y,  # negative coefficient: (-2)**y is not real
        T**Y,
    ],
)
def test_power_unchanged(expr):
    assert split_power(expr) == expr


def test_power_value_is_unchanged():
    expr = (1.18e-10 * sp.exp(-69500 / T) * sp.sqrt(T)) ** (1 / (1 + X))
    rewritten = split_power(expr)

    for t in (2000.0, 5000.0):
        point = {T: t, X: 0.3}
        assert float(rewritten.subs(point)) == pytest.approx(float(expr.subs(point)))


def test_derivative_has_no_growing_exponential_after_split():
    # d(B**e)/dT has B'/B: exp(-E/T)'/exp(-E/T) leaves exp(+E/T), overflowing below
    # ~98 K here (SymPy cancels it only when B is just c*exp(-E/T)).
    expr = (1.18e-10 * sp.exp(-69500 / T) * sp.sqrt(T)) ** (1 / (1 + X))
    namespace = {"math": math, "tgas": 90.0, "x": 0.3}

    with pytest.raises(OverflowError):
        eval(sp.pycode(sp.diff(expr, T)), namespace)
    exact = sp.diff(expr, T).subs({T: 90, X: sp.Rational(3, 10)}).evalf(50)
    assert eval(sp.pycode(sp.diff(split_power(expr), T)), namespace) == pytest.approx(
        float(exact), rel=1e-13
    )
