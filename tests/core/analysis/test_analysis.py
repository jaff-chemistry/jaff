# ABOUTME: Tests for net.analysis: enabling and disabling, the parameter, density
# ABOUTME: and function tables built from the network, set/get/pop and labels

import numpy as np
import pytest
import sympy as sp

from jaff.core.analysis import (
    Analysis,
    AnalysisSpec,
    Densities,
    Entry,
    Functions,
    Parameters,
)

LINES = (
    "H+ + E -> H  []  3.61e-12*(Tgas/3e2)**(-0.75)\n"
    "C -> C+ + E  []  3.39e-10*fshield(ncol_C, av)*exp(-av)\n"
    "CO -> C + O  []  1e-10*chi\n"
)


@pytest.fixture
def net(make_network):
    return make_network(LINES, name="u.jet", funcfile=False)


@pytest.fixture
def on(net):
    net.analysis.enable()

    return net.analysis


# --------------------------------------------------------------------------- #
# Enabling                                                                     #
# --------------------------------------------------------------------------- #


def test_every_network_has_a_disabled_analysis(net):
    assert isinstance(net.analysis, Analysis)
    assert net.analysis is net.analysis
    assert not net.analysis.enabled
    assert net.analysis.parameters is None
    assert net.analysis.densities is None
    assert net.analysis.functions is None
    assert net.analysis.get_spec() is None


def test_enable_builds_the_tables_from_the_network(net, on):
    assert on.enabled
    assert isinstance(on.parameters, Parameters)
    assert isinstance(on.densities, Densities)
    assert isinstance(on.functions, Functions)
    assert on.parameters.names == ("av", "chi", "ncol_c", "tgas")
    assert on.densities.names == tuple(s.name for s in net.species)
    assert on.functions.names == ("fshield",)


def test_new_tables_are_empty(on):
    assert on.parameters.get("tgas") == Entry(None, None)
    assert on.parameters.values() == {}


def test_enable_sets_the_given_entries(net):
    net.analysis.enable(
        parameters={"tgas": 20, "av": {"value": 10, "label": "$A_V$"}},
        densities={"H": 1.0},
        functions={"fshield": 0.5},
    )

    assert net.analysis.parameters.get("tgas") == Entry(20.0, None)
    assert net.analysis.parameters.get("av") == Entry(10.0, "$A_V$")
    assert net.analysis.densities.get("H").value == pytest.approx(1.0)
    assert net.analysis.functions.get("fshield").value == pytest.approx(0.5)


def test_enable_again_merges_over_what_is_stored(net):
    net.analysis.enable(parameters={"tgas": 20, "av": 1})
    net.analysis.enable(parameters={"av": 2})

    assert net.analysis.parameters.values() == {"tgas": 20.0, "av": 2.0}


def test_disable_drops_the_tables(net, on):
    on.parameters.set("tgas", 20)
    net.analysis.disable()

    assert not net.analysis.enabled
    assert net.analysis.parameters is None
    net.analysis.enable()
    assert net.analysis.parameters.values() == {}


def test_densities_from_an_array_in_species_order(net):
    dens = np.arange(len(net.species), dtype=float)
    net.analysis.enable(densities=dens)

    for s in net.species:
        assert net.analysis.densities.get(s.name).value == pytest.approx(s.index)


def test_density_array_of_the_wrong_length_raises(net):
    with pytest.raises(ValueError, match="species"):
        net.analysis.enable(densities=[1.0, 2.0])


# --------------------------------------------------------------------------- #
# set / get / pop                                                              #
# --------------------------------------------------------------------------- #


def test_unknown_name_raises_with_a_hint(on):
    with pytest.raises(KeyError, match="tgas"):
        on.parameters.set("tgass", 1.0)


def test_label_true_gives_the_default_label(on):
    on.parameters.set("tgas", 20, label=True)
    on.densities.set("H", 1, label=True)
    on.functions.set("fshield", 0.5, label=True)

    assert on.parameters.get("tgas").label == "$tgas = {}$"
    assert on.densities.get("H").label == r"$n({\rm H}) = {}$"
    assert on.functions.get("fshield").label == "fshield = {}"


def test_set_without_a_label_keeps_the_label(on):
    on.parameters.set("tgas", 20, label="T")
    on.parameters.set("tgas", 30)

    assert on.parameters.get("tgas") == Entry(30.0, "T")


def test_label_false_removes_the_label(on):
    on.parameters.set("tgas", 20, label="T")
    on.parameters.set("tgas", 20, label=False)

    assert on.parameters.get("tgas").label is None


def test_get_without_a_name_is_a_copy_of_every_entry(on):
    on.parameters.set("tgas", 20)
    entries = on.parameters.get()
    entries["tgas"] = Entry(1.0, None)

    assert set(entries) == set(on.parameters.names)
    assert on.parameters.get("tgas").value == pytest.approx(20.0)


def test_pop_one_clears_it_and_returns_the_old_entry(on):
    on.parameters.set("tgas", 20, label="T")

    assert on.parameters.pop("tgas") == Entry(20.0, "T")
    assert on.parameters.get("tgas") == Entry(None, None)


def test_pop_without_a_name_clears_the_table(on):
    on.parameters.update({"tgas": 20, "av": 1})
    on.densities.set("H", 1.0)
    on.parameters.pop()

    assert on.parameters.values() == {}
    assert on.densities.values() == {"H": 1.0}


def test_analysis_pop_clears_every_table(on):
    on.parameters.set("tgas", 20)
    on.densities.set("H", 1.0)
    on.functions.set("fshield", 0.5)
    on.pop()

    assert on.parameters.values() == on.densities.values() == on.functions.values() == {}
    assert on.enabled


def test_missing_lists_the_unset_names(on):
    on.parameters.set("tgas", 20)

    assert on.parameters.missing(["tgas", "av", "chi"]) == {"av", "chi"}


def test_update_rejects_unknown_entry_fields(on):
    with pytest.raises(ValueError, match="units"):
        on.parameters.update({"tgas": {"value": 1, "units": "K"}})


# --------------------------------------------------------------------------- #
# Values                                                                       #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("value", [np.nan, np.inf, "hot"])
def test_parameters_must_be_finite_numbers(on, value):
    with pytest.raises(ValueError, match="tgas"):
        on.parameters.set("tgas", value)


def test_densities_must_not_be_negative(on):
    with pytest.raises(ValueError, match="H"):
        on.densities.set("H", -1.0)


def test_a_function_can_be_a_float_an_expression_or_a_callable(on):
    av = sp.Symbol("av")
    on.functions.set("fshield", 0.5)
    on.functions.set("fshield", sp.exp(-av))
    on.functions.set("fshield", lambda: 0.25)

    assert callable(on.functions.get("fshield").value)


def test_a_function_expression_may_only_use_parameters(on):
    with pytest.raises(ValueError, match="foo"):
        on.functions.set("fshield", sp.Symbol("foo") * 2)


def test_a_function_expression_may_not_use_unknown_functions(on):
    with pytest.raises(ValueError, match="g"):
        on.functions.set("fshield", sp.Function("g")(sp.Symbol("av")))


# --------------------------------------------------------------------------- #
# Spec                                                                         #
# --------------------------------------------------------------------------- #


def test_spec_is_a_read_only_snapshot(on):
    on.parameters.set("tgas", 20)
    spec = on.get_spec()
    on.parameters.set("tgas", 30)

    assert isinstance(spec, AnalysisSpec)
    assert spec.parameters["tgas"].value == pytest.approx(20.0)
    with pytest.raises(TypeError):
        spec.parameters["tgas"] = Entry(1.0, None)
