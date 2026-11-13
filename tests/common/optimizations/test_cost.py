# ABOUTME: Tests for optimizations.cost: a weighted operation count of an expression,
# ABOUTME: cheap arithmetic, dearer division, expensive pow/exp/log/function calls

import pytest
import sympy as sp

from jaff.common.optimizations import cost

X, Y, T, K = sp.symbols("x y tgas k")


@pytest.mark.parametrize(
    ("expr", "expected"),
    [
        (X + Y, 1),
        (X * Y, 1),
        (-X, 1),
        (X / Y, 4),
        (X**2, 2),
        (sp.Float(3.0), 0),
        (X, 0),
    ],
    ids=["add", "mul", "neg", "div", "integer-pow", "number", "symbol"],
)
def test_arithmetic(expr, expected):
    assert cost(expr) == expected


@pytest.mark.parametrize(
    "expr",
    [sp.exp(X), sp.log(X), X**-1.8, sp.sqrt(X), sp.Function("fshield")(X)],
    ids=["exp", "log", "real-pow", "sqrt", "unknown-function"],
)
def test_transcendental_calls_are_expensive(expr):
    assert cost(expr) >= 20


def test_a_power_law_derivative_costs_more_than_its_log_form():
    plain = sp.diff(4.12e-10 * T**-0.8, T)  # -3.3e-10/tgas**1.8: a new pow

    assert cost(-0.8 * K / T) < cost(plain)


def test_costs_add_up():
    assert cost(sp.exp(X) + X * Y) == cost(sp.exp(X)) + cost(X * Y) + 1
