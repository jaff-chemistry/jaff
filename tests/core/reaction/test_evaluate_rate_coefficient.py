# ABOUTME: Tests for evaluating rate coefficients numerically: values for symbols,
# ABOUTME: densities and unknown functions; every missing name reported at once

import math

import numpy as np
import pytest
import sympy as sp

from jaff import Network
from jaff.core.reaction._evaluate import evaluate_expression

NETWORK = (
    "H+ + E -> H  []  3.61e-12*(Tgas/3e2)**(-0.75)\n"
    "H -> H+ + E  []  4.60e-01*crate\n"
    "C -> C+ + E  []  3.39e-10 * exp(-3.761 * av)\n"
    "CO+ + E -> CO  []  2e-14\n"
)
VALUES = {"tgas": 100.0, "crate": 1.3e-17, "av": 2.0}


@pytest.fixture
def net(make_network) -> Network:
    return make_network(NETWORK, name="r.jet")


def _rxn(net, serialized):
    (rxn,) = net.reactions.all(serialized)

    return rxn


# --------------------------------------------------------------------------- #
# Reaction                                                                     #
# --------------------------------------------------------------------------- #


def test_reaction_value(net):
    rxn = _rxn(net, "H+.e-__H")

    assert rxn.evaluate_rate_coefficient(VALUES) == pytest.approx(
        3.61e-12 * (100.0 / 300.0) ** -0.75
    )


def test_constant_rate_needs_no_values(net):
    assert _rxn(net, "CO+.e-__CO").evaluate_rate_coefficient({}) == pytest.approx(2e-14)


def test_unused_values_are_ignored(net):
    rxn = _rxn(net, "C__C+.e-")

    assert rxn.evaluate_rate_coefficient(VALUES | {"chi": 5.0}) == pytest.approx(
        3.39e-10 * math.exp(-3.761 * 2.0)
    )


def test_missing_value_is_named(net):
    with pytest.raises(ValueError, match="av"):
        _rxn(net, "C__C+.e-").evaluate_rate_coefficient({"tgas": 100.0})


def test_result_is_a_float(net):
    value = _rxn(net, "H__H+.e-").evaluate_rate_coefficient(VALUES)

    assert type(value) is float
    assert value == pytest.approx(4.60e-01 * 1.3e-17)


# --------------------------------------------------------------------------- #
# Reactions                                                                    #
# --------------------------------------------------------------------------- #


def test_reactions_values_in_catalogue_order(net):
    values = net.reactions.evaluate_rate_coefficients(VALUES)
    expected = [r.evaluate_rate_coefficient(VALUES) for r in net.reactions]

    assert isinstance(values, np.ndarray)
    assert values.shape == (len(net.reactions),)
    assert values == pytest.approx(expected)


def test_reactions_report_every_missing_name_at_once(net):
    with pytest.raises(ValueError) as info:
        net.reactions.evaluate_rate_coefficients({"tgas": 100.0})

    assert "av" in str(info.value) and "crate" in str(info.value)


# --------------------------------------------------------------------------- #
# Expressions: densities and unknown functions                                 #
# --------------------------------------------------------------------------- #

T = sp.Symbol("tgas")
NDEN = sp.IndexedBase("nden")
F = sp.Function("fshield")


def test_densities_from_an_array():
    expr = 2 * NDEN[0] + NDEN[2]

    assert evaluate_expression(expr, {"nden": [1.0, 5.0, 3.0]}) == pytest.approx(5.0)


def test_missing_density_array_is_named():
    with pytest.raises(ValueError, match="nden"):
        evaluate_expression(NDEN[1] * T, {"tgas": 1.0})


def test_function_as_a_callable_gets_numeric_arguments():
    seen = []

    def fshield(x, y):
        seen.append((x, y))
        return 0.5

    expr = 1e-10 * F(T, 2 * T)
    value = evaluate_expression(expr, {"tgas": 10.0}, {"fshield": fshield})

    assert value == pytest.approx(5e-11)
    assert seen == [(pytest.approx(10.0), pytest.approx(20.0))]


def test_function_as_a_constant():
    value = evaluate_expression(1e-10 * F(T), {"tgas": 10.0}, {"fshield": 0.25})

    assert value == pytest.approx(2.5e-11)


def test_nested_functions_are_evaluated_inside_out():
    g = sp.Function("g")
    value = evaluate_expression(
        F(g(T)), {"tgas": 3.0}, {"g": lambda x: 2 * x, "fshield": lambda x: x + 1}
    )

    assert value == pytest.approx(7.0)


def test_missing_function_is_named():
    with pytest.raises(ValueError, match="fshield"):
        evaluate_expression(F(T), {"tgas": 1.0})


def test_piecewise_picks_the_branch():
    expr = sp.Piecewise((1.0, T < 100), (2.0, True))

    assert evaluate_expression(expr, {"tgas": 50.0}) == pytest.approx(1.0)
    assert evaluate_expression(expr, {"tgas": 500.0}) == pytest.approx(2.0)


def test_gow_with_densities_and_a_shielding_function():
    net = Network("GOW")
    names = {
        "tgas": 50.0,
        "av": 1.0,
        "chi": 1.0,
        "crate": 2e-16,
        "d2g": 0.01,
        "vdisp": 1e5,
        "ncol_C": 1e14,
        "ncol_CO": 1e14,
        "ncol_H2": 1e20,
        "nden": np.full(len(net.species), 1.0),
    }
    values = net.reactions.evaluate_rate_coefficients(names, {"fshield_co_interp2d": 0.1})

    assert np.all(np.isfinite(values))
    assert np.all(values >= 0)
