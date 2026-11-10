# ABOUTME: Tests for net.analysis.evaluate: given values merged over the analysis
# ABOUTME: tables, densities, functions, missing names, and the shape of the result

import math

import numpy as np
import pytest
import sympy as sp

from jaff import Network

NETWORK = (
    "H+ + E -> H  []  3.61e-12*(Tgas/3e2)**(-0.75)\n"
    "H -> H+ + E  []  4.60e-01*crate\n"
    "C -> C+ + E  []  3.39e-10 * exp(-3.761 * av)\n"
    "CO+ + E -> CO  []  2e-14\n"
)
PARAMETERS = {"tgas": 100.0, "crate": 1.3e-17, "av": 2.0}

FUNCTIONS = "C -> C+ + E  []  3.39e-10*fshield(ncol_C, av)\nCO -> C + O  []  1e-10*chi\n"


@pytest.fixture
def net(make_network) -> Network:
    return make_network(NETWORK, name="r.jet")


@pytest.fixture
def fnet(make_network) -> Network:
    return make_network(FUNCTIONS, name="f.jet", funcfile=False)


def _rxn(net, serialized):
    (rxn,) = net.reactions.all(serialized)

    return rxn


def _k(rxn, *args, **kwargs):
    """The rate coefficient of *rxn*, evaluated."""
    return rxn.net.analysis.evaluate(rxn.rate_coefficient, *args, **kwargs)


def _rate(rxn, *args, **kwargs):
    """The rate (coefficient times reactant densities) of *rxn*, evaluated."""
    return rxn.net.analysis.evaluate(rxn.rate.expanded, *args, **kwargs)


# --------------------------------------------------------------------------- #
# Rate coefficients from given values                                          #
# --------------------------------------------------------------------------- #


def test_reaction_value(net):
    rxn = _rxn(net, "H+.e-__H")

    assert _k(rxn, PARAMETERS) == pytest.approx(3.61e-12 * (100.0 / 300.0) ** -0.75)


def test_constant_rate_needs_no_values(net):
    assert _k(_rxn(net, "CO+.e-__CO")) == pytest.approx(2e-14)


def test_unused_values_are_ignored(net):
    rxn = _rxn(net, "C__C+.e-")

    assert _k(rxn, PARAMETERS | {"chi": 5.0}) == pytest.approx(
        3.39e-10 * math.exp(-3.761 * 2.0)
    )


def test_missing_parameter_is_named(net):
    with pytest.raises(ValueError, match="parameter av"):
        _k(_rxn(net, "C__C+.e-"), {"tgas": 100.0})


def test_result_is_a_float(net):
    value = _k(_rxn(net, "H__H+.e-"), PARAMETERS)

    assert type(value) is float
    assert value == pytest.approx(4.60e-01 * 1.3e-17)


def test_nden_in_parameters_points_to_densities(net):
    with pytest.raises(ValueError, match="densities"):
        _k(_rxn(net, "H+.e-__H"), PARAMETERS | {"nden": np.ones(len(net.species))})


# --------------------------------------------------------------------------- #
# net.analysis                                                                 #
# --------------------------------------------------------------------------- #


def test_values_come_from_the_analysis(net):
    net.analysis.enable(parameters=PARAMETERS)
    rxn = _rxn(net, "C__C+.e-")

    assert _k(rxn) == pytest.approx(_k(rxn, PARAMETERS))


def test_given_values_merge_over_the_analysis(net):
    net.analysis.enable(parameters=PARAMETERS)
    rxn = _rxn(net, "H+.e-__H")

    assert _k(rxn, {"tgas": 300.0}) == pytest.approx(3.61e-12)
    assert net.analysis.parameters.get("tgas").value == pytest.approx(100.0)


def test_disabled_analysis_is_ignored(net):
    net.analysis.enable(parameters=PARAMETERS)
    net.analysis.disable()

    with pytest.raises(ValueError, match="tgas"):
        _k(_rxn(net, "H+.e-__H"))


# --------------------------------------------------------------------------- #
# Rates                                                                        #
# --------------------------------------------------------------------------- #


def test_rate_is_the_coefficient_times_the_reactant_densities(net):
    rxn = _rxn(net, "H+.e-__H")
    dens = {"H+": 2.0, "e-": 5.0}

    assert _rate(rxn, PARAMETERS, dens) == pytest.approx(_k(rxn, PARAMETERS) * 10.0)


def test_rate_densities_from_an_array(net):
    rxn = _rxn(net, "H+.e-__H")
    dens = np.full(len(net.species), 3.0)

    assert _rate(rxn, PARAMETERS, dens) == pytest.approx(_k(rxn, PARAMETERS) * 9.0)


def test_missing_density_is_named_by_species(net):
    with pytest.raises(ValueError, match=r"no density for e-, H\+"):
        _rate(_rxn(net, "H+.e-__H"), PARAMETERS)


def test_density_array_of_the_wrong_length_raises(net):
    with pytest.raises(ValueError, match="species"):
        _rate(_rxn(net, "H+.e-__H"), PARAMETERS, [1.0])


def test_rate_densities_from_the_analysis(net):
    net.analysis.enable(parameters=PARAMETERS, densities={"H+": 2.0, "e-": 5.0})
    rxn = _rxn(net, "H+.e-__H")

    assert _rate(rxn) == pytest.approx(_k(rxn) * 10.0)


def test_given_densities_merge_over_the_analysis(net):
    net.analysis.enable(parameters=PARAMETERS, densities={"H+": 2.0, "e-": 5.0})
    rxn = _rxn(net, "H+.e-__H")

    assert _rate(rxn, densities={"e-": 1.0}) == pytest.approx(_k(rxn) * 2.0)


# --------------------------------------------------------------------------- #
# Unknown functions                                                            #
# --------------------------------------------------------------------------- #


def test_function_as_a_constant(fnet):
    rxn = _rxn(fnet, "C__C+.e-")

    assert _k(rxn, {}, functions={"fshield": 0.5}) == pytest.approx(3.39e-10 * 0.5)


def test_function_arguments_are_ignored(fnet):
    # fshield(ncol_c, av): no value for ncol_c or av is needed.
    rxn = _rxn(fnet, "C__C+.e-")

    assert _k(rxn, functions={"fshield": 0.5}) == pytest.approx(1.695e-10)


def test_function_as_an_expression_of_parameters(fnet):
    rxn = _rxn(fnet, "C__C+.e-")
    value = _k(rxn, {"av": 2.0}, functions={"fshield": sp.exp(-sp.Symbol("av"))})

    assert value == pytest.approx(3.39e-10 * math.exp(-2.0))


def test_function_as_a_callable_without_arguments(fnet):
    rxn = _rxn(fnet, "C__C+.e-")

    assert _k(rxn, functions={"fshield": lambda: 0.25}) == pytest.approx(3.39e-10 * 0.25)


def test_a_function_expression_needs_its_parameters(fnet):
    with pytest.raises(ValueError, match="parameter av"):
        _k(_rxn(fnet, "C__C+.e-"), functions={"fshield": sp.exp(-sp.Symbol("av"))})


def test_a_function_value_may_not_use_unknown_functions(fnet):
    with pytest.raises(ValueError, match="g"):
        _k(_rxn(fnet, "C__C+.e-"), functions={"fshield": lambda: sp.Function("g")(1)})


def test_missing_function_is_named(fnet):
    with pytest.raises(ValueError, match="function fshield"):
        _k(_rxn(fnet, "C__C+.e-"))


def test_functions_from_the_analysis(fnet):
    fnet.analysis.enable(functions={"fshield": 0.5})

    assert _k(_rxn(fnet, "C__C+.e-")) == pytest.approx(1.695e-10)


# --------------------------------------------------------------------------- #
# Many expressions                                                             #
# --------------------------------------------------------------------------- #


def test_reactions_values_in_catalogue_order(net):
    values = net.analysis.evaluate(net.reactions.rate_coefficients(), PARAMETERS)
    expected = [_k(r, PARAMETERS) for r in net.reactions]

    assert isinstance(values, np.ndarray)
    assert values.shape == (len(net.reactions),)
    assert values == pytest.approx(expected)


def test_reactions_rates_in_catalogue_order(net):
    dens = np.full(len(net.species), 2.0)
    values = net.analysis.evaluate(net.reactions.rates(), PARAMETERS, dens)
    expected = [_rate(r, PARAMETERS, dens) for r in net.reactions]

    assert values == pytest.approx(expected)


def test_reactions_report_every_missing_name_at_once(net):
    with pytest.raises(ValueError) as info:
        net.analysis.evaluate(net.reactions.rates(), {"tgas": 100.0})

    message = str(info.value)
    assert "parameter av, crate" in message
    assert "no density for" in message and "CO+" in message


def test_reactions_use_the_analysis(net):
    net.analysis.enable(parameters=PARAMETERS)

    assert net.analysis.evaluate(net.reactions.rate_coefficients()) == pytest.approx(
        net.analysis.evaluate(net.reactions.rate_coefficients(), PARAMETERS)
    )


def test_gow_with_densities_and_a_shielding_function():
    net = Network("GOW")
    parameters = {
        "tgas": 50.0,
        "av": 1.0,
        "chi": 1.0,
        "crate": 2e-16,
        "d2g": 0.01,
        "vdisp": 1e5,
        "ncol_C": 1e14,
        "ncol_CO": 1e14,
        "ncol_H2": 1e20,
    }
    values = net.analysis.evaluate(
        net.reactions.rates(),
        parameters,
        np.full(len(net.species), 1.0),
        {"fshield_co_interp2d": 0.1},
    )

    assert np.all(np.isfinite(values))
    assert np.all(values >= 0)


def test_list_gives_an_array(net):
    exprs = [r.rate_coefficient for r in net.reactions]
    values = net.analysis.evaluate(exprs, PARAMETERS)

    assert isinstance(values, np.ndarray)
    assert values == pytest.approx([_k(r, PARAMETERS) for r in net.reactions])


def test_empty_sequence_gives_an_empty_array(net):
    values = net.analysis.evaluate([], PARAMETERS)

    assert isinstance(values, np.ndarray)
    assert values.shape == (0,)


def test_number_gives_a_float(net):
    assert net.analysis.evaluate(sp.Integer(3)) == pytest.approx(3.0)
    assert type(net.analysis.evaluate(2.5)) is float


def test_matrix_keeps_its_shape(net):
    t = sp.Symbol("tgas")
    matrix = sp.Matrix([[t, 2 * t, 0], [1, t**2, 3]])
    values = net.analysis.evaluate(matrix, {"tgas": 2.0})

    assert values.shape == (2, 3)
    assert values == pytest.approx(np.array([[2.0, 4.0, 0.0], [1.0, 4.0, 3.0]]))


def test_sparse_matrix_gives_a_dense_array(net):
    t = sp.Symbol("tgas")
    matrix = sp.ImmutableSparseMatrix(3, 3, {(0, 1): t, (2, 0): 2 * t})
    values = net.analysis.evaluate(matrix, {"tgas": 5.0})

    expected = np.zeros((3, 3))
    expected[0, 1], expected[2, 0] = 5.0, 10.0
    assert values == pytest.approx(expected)


def test_missing_values_in_a_matrix_are_named(net):
    matrix = sp.Matrix([[sp.Symbol("tgas")], [sp.Symbol("av")]])

    with pytest.raises(ValueError, match="parameter av, tgas"):
        net.analysis.evaluate(matrix)


def test_ode_system_rhs_and_jacobian(make_network):
    net = make_network(NETWORK, name="r.jet")
    system = net.ode_system(thermal="dtdt")
    dens = np.full(len(net.species), 2.0)
    rhs = net.analysis.evaluate(system.rhs, PARAMETERS, dens)
    jac = net.analysis.evaluate(system.jacobian, PARAMETERS, dens)

    assert rhs.shape == (len(system.state), 1)
    assert jac.shape == (len(system.state), len(system.state))
    assert np.all(np.isfinite(jac))


def test_reduced_form_asks_for_the_expanded_one(net):
    rxn = _rxn(net, "H+.e-__H")

    with pytest.raises(ValueError, match="krate.*expanded"):
        net.analysis.evaluate(rxn.rate.reduced, PARAMETERS, {"H+": 1.0, "e-": 1.0})
