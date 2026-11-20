# ABOUTME: Tests for optimizations.split_log: log(c*exp(b)*R) as log(c) + b + log(R), so
# ABOUTME: an exp that underflows at low temperature no longer makes a log(0)

import math

import pytest
import sympy as sp

from jaff.common.optimizations import PRE_CSE_RULES, bottom_up, split_log

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
