# ABOUTME: net.functions: the .jfunc functions, used by calling them with arguments
# ABOUTME: (call) or reading their parameters as network inputs (value)

import pytest
import sympy as sp

from jaff.core.network import NetworkFunctions
from jaff.errors import ParserError

LINES = "H + E -> H-  []  1e-10\nH2 + C -> CH + H  []  1e-10\nCH -> C + H  []  1e-9\n"


@pytest.fixture
def func_network(make_network, tmp_path):
    """Network from *lines* with *jfunc* as its ``.jfunc`` file."""

    def _make(jfunc, lines=LINES, **kwargs):
        (tmp_path / "f.jfunc").write_text(jfunc)
        return make_network(lines, name="f.jet", **kwargs)

    return _make


def _std(net, name):
    return net.symbols.standardize(sp.Symbol(name))


class TestValue:
    def test_parameters_are_network_inputs(self, func_network):
        net = func_network("@function chemRate0(tgas, n_e)\n@return 1e-10*n_e*tgas\n")
        value = net.functions.value("chemRate0")

        assert value == 1e-10 * _std(net, "n_e") * sp.Symbol("tgas")

    def test_name_is_case_insensitive(self, func_network):
        net = func_network("@function chemRate0(n_e)\n@return 2*n_e\n")

        assert net.functions.value("CHEMRATE0") == net.functions.value("chemrate0")

    def test_is_cached(self, func_network):
        net = func_network("@function chemRate0(n_e)\n@return 2*n_e\n")

        assert net.functions.value("chemrate0") is net.functions.value("chemrate0")

    def test_nested_call_is_standardized(self, func_network):
        jfunc = (
            "@function cool(n_CO, tgas)\n@return 1e-24*n_CO*tgas*n_H\n"
            "@function heatingCoolingRate(tgas, n_CO)\n@return -cool(n_CO, tgas)\n"
        )
        net = func_network(jfunc, lines=LINES + "CO -> C + O  []  1e-9\n")
        value = net.functions.value("heatingcoolingrate")
        tgas = sp.Symbol("tgas")

        assert value == -1e-24 * _std(net, "n_CO") * tgas * _std(net, "n_H")


class TestCall:
    def test_arguments_replace_parameters(self, func_network):
        net = func_network("@function f(x)\n@return 3*x*n_H\n")

        assert net.functions.call("f", [sp.Integer(2)]) == 6 * _std(net, "n_H")

    def test_parameter_named_like_a_density_is_a_placeholder(self, func_network):
        net = func_network("@function g(n_C)\n@return 2*n_C\n")
        arg = sp.Symbol("x")

        assert net.functions.call("g", [arg]) == 2 * arg

    def test_wrong_argument_count_raises(self, func_network):
        net = func_network("@function f(x)\n@return 3*x\n")

        with pytest.raises(ParserError, match="argument"):
            net.functions.call("f", [sp.Integer(1), sp.Integer(2)])

    def test_global_calls_helper(self, func_network):
        jfunc = "@function shield(NH2)\n@return 1/(1 + 1e-20*NH2)\n"
        lines = "@var: s = shield(1e21*n_C)\n" + LINES + "H + H -> H2  []  1e-17*s\n"
        net = func_network(jfunc, lines=lines)
        rxn = net.reactions[-1]

        n_c = _std(net, "n_C")
        value = float(rxn.rate_coefficient.subs(n_c, 2.0))

        assert value == pytest.approx(1e-17 / 21.0, rel=1e-12)


class TestTable:
    def test_is_on_the_network(self, func_network):
        net = func_network("@function f(x)\n@return x\n")

        assert isinstance(net.functions, NetworkFunctions)
        assert "F" in net.functions
        assert "g" not in net.functions

    def test_raw_table_is_not_on_the_network(self, func_network):
        net = func_network("@function f(x)\n@return x\n")

        assert not hasattr(net, "aux_funcs")

    def test_unused_helper_with_unknown_species_loads(self, func_network):
        net = func_network("@function unused(x)\n@return x*n_Xe\n")

        assert "unused" in net.functions


class TestReactionFunction:
    def test_reads_the_reaction_function(self, func_network):
        net = func_network("@function deltaE0(n_e)\n@return 2*n_e\n")

        assert net.reactions[0].function("deltaE") == 2 * _std(net, "n_e")
        assert net.reactions[0].dE == 2 * _std(net, "n_e")

    def test_none_without_one(self, func_network):
        net = func_network("")

        assert net.reactions[0].function("deltaE") is None

    def test_aux_funcs_holder_is_gone(self, func_network):
        net = func_network("")

        assert not hasattr(net.reactions[0], "aux_funcs")

    def test_chem_rate_parameters_are_standardized(self, func_network):
        net = func_network("@function chemRate0(tgas, n_H2)\n@return 1e-10*n_H2\n")

        assert net.reactions[0].rate_coefficient == 1e-10 * _std(net, "n_H2")


class TestRateReferences:
    def test_rate_reference_resolves_to_target_rate(self, func_network):
        net = func_network("@function chemRate1(n_e)\n@return 2*rc_2*n_e\n")
        target = net.reactions.by_source_index(2)

        assert net.reactions.by_source_index(1).rate_coefficient == (
            2 * target.rate_coefficient * _std(net, "n_e")
        )

    def test_rate_reference_to_a_later_reaction_is_standardized(self, func_network):
        # rc_1 points at a reaction later in the file whose rate uses n_H.
        lines = "H + E -> H-  []  2*rc_1\nCH -> C + H  []  1e-9*n_H\n"
        net = func_network("", lines=lines)

        free = {str(s) for s in net.reactions[0].rate_coefficient.free_symbols}
        assert not any(name.startswith(("n_", "rc_")) for name in free)

    def test_rate_reference_in_cross_section_raises(self, func_network):
        jfunc = "@function photoCrossSection0()\n@return 1e-18*rc_1\n"
        lines = "H -> H+ + E  []  PHOTO, 13.6\nH+ + E -> H  []  1e-12\n"

        with pytest.raises(ParserError, match="photoCrossSection0"):
            func_network(jfunc, lines=lines)
