# ABOUTME: Reactions.link(): rc_<N> references resolved to reaction N's final rate
# ABOUTME: coefficient in dependency order, with cycle, missing and malformed errors

import pytest
import sympy as sp

from jaff.errors import ParserError

H = "H + E -> H-  []  {}\n"
C = "C + E -> C-  []  {}\n"
OX = "O + E -> O-  []  {}\n"


@pytest.fixture
def link_network(make_network, tmp_path):
    """Network from *lines* with *jfunc* as its ``.jfunc`` file."""

    def _make(lines, jfunc="", **kwargs):
        (tmp_path / "l.jfunc").write_text(jfunc)
        return make_network(lines, name="l.jet", **kwargs)

    return _make


def _rate(net, index):
    return net.reactions.by_source_index(index).rate.coefficient


class TestLink:
    def test_chain_resolves_in_dependency_order(self, link_network):
        net = link_network(H.format("2*rc_1") + C.format("3*rc_2") + OX.format("5.0"))

        assert float(_rate(net, 0)) == pytest.approx(30.0)
        assert float(_rate(net, 1)) == pytest.approx(15.0)

    def test_dE_may_use_its_own_rate(self, link_network):
        net = link_network(H.format("4.0"), "@function deltaE0()\n@return 2*rc_0\n")

        assert float(net.reactions[0].dE) == pytest.approx(8.0)

    def test_dE_uses_the_final_target_rate(self, link_network):
        jfunc = "@function deltaE0()\n@return rc_1\n"
        net = link_network(H.format("1.0") + C.format("3*rc_2") + OX.format("5.0"), jfunc)

        assert float(net.reactions[0].dE) == pytest.approx(15.0)

    def test_group_yield_uses_the_target_rate(self, link_network):
        jfunc = "@function spectralYield0()\n@return 1e-12*rc_1\n"
        net = link_network(
            H.format("1.0") + C.format("2.0"),
            jfunc,
            radiation={"bands": [6.0, 13.6, 100.0], "mode": "u", "profile_index": 0.5},
        )

        yields = [float(g) for g in net.reactions[0].rad.group_yields]
        assert yields == pytest.approx([2e-12 * 7.6, 2e-12 * 86.4], rel=1e-6)

    def test_network_expression_resolves(self, link_network):
        jfunc = "@function heatingCoolingRate()\n@return 7*rc_0\n"
        net = link_network(H.format("2.0"), jfunc)

        assert float(net.thermodynamics.dEdt_extra.volumetric) == pytest.approx(14.0)


class TestErrors:
    def test_cycle_reports_its_path(self, link_network):
        with pytest.raises(ParserError, match="rc_0 → rc_1 → rc_0"):
            link_network(H.format("2*rc_1") + C.format("3*rc_0"))

    def test_self_reference_in_a_rate_is_a_cycle(self, link_network):
        with pytest.raises(ParserError, match="rc_0 → rc_0"):
            link_network(H.format("2*rc_0"))

    def test_missing_target_raises(self, link_network):
        with pytest.raises(ParserError, match="not in the network"):
            link_network(H.format("2*rc_9"))

    def test_malformed_reference_raises(self, link_network):
        with pytest.raises(ParserError, match="rc_x"):
            link_network(H.format("2*rc_x"))


class TestOwnership:
    def test_reactions_resolve_rate_references(self, link_network):
        net = link_network(H.format("2.0"))

        assert net.reactions.resolve_rate_references(3 * sp.Symbol("rc_0")) == 6.0

    def test_symbols_do_not_resolve_rate_references(self, link_network):
        net = link_network(H.format("2.0"))

        assert not hasattr(net.symbols, "resolve_rate_references")

    def test_reaction_substitute_reaches_every_expression(self, link_network):
        jfunc = "@function deltaE0()\n@return 2*x\n"
        net = link_network(H.format("3*x"), jfunc)
        rxn = net.reactions[0]
        x = sp.Symbol("x")

        rxn.substitute({x: sp.Float(1.0)})

        assert all(not e.has(x) for e in rxn.expressions())
        assert float(rxn.rate.coefficient) == pytest.approx(3.0)
        assert float(rxn.dE) == pytest.approx(2.0)
