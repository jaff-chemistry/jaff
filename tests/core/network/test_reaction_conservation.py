# ABOUTME: Conservation (mass/charge) enforcement through Network loading
# ABOUTME: strict=True must abort; non-strict must warn naming the reaction

import logging

import pytest


# H <-> H2 conserves every species (no sink/source), but each reaction breaks
# mass conservation (1 H vs 2 H).  This isolates the per-reaction check from the
# later network-wide sink/source/recombination/isomer diagnostics.
_MASS_VIOLATION = ["H -> H2 [10,1000] 1", "H2 -> H [10,1000] 1"]


class TestStrictTextLoading:
    def test_mass_violation_exits(self, make_network):
        with pytest.raises(SystemExit):
            make_network(_MASS_VIOLATION, strict=True)

    def test_non_strict_loads(self, make_network):
        # strict=False must not abort on a conservation failure.
        net = make_network(_MASS_VIOLATION, strict=False)
        assert net.reactions.count == 2


class TestReactionCheck:
    def test_strict_check_exits(self, make_network):
        net = make_network(_MASS_VIOLATION, strict=False)

        with pytest.raises(SystemExit):
            net.reactions[0].check(strict=True)

    def test_non_strict_check_returns(self, make_network):
        net = make_network(_MASS_VIOLATION, strict=False)

        assert net.reactions[0].check(strict=False) is None


class TestNetworkCheck:
    def test_strict_check_exits(self, make_network):
        net = make_network(_MASS_VIOLATION, strict=False)

        with pytest.raises(SystemExit):
            net.check_conservation(strict=True)

    def test_non_strict_check_returns(self, make_network):
        net = make_network(_MASS_VIOLATION, strict=False)

        assert net.check_conservation(strict=False) is None

    def test_constructor_does_not_check(self, make_network):
        # A hand-built reaction is checked by the network, not on construction.
        from jaff.core.reaction import Reaction

        net = make_network(
            ["H + H -> H2 [10,1000] 1", "H2 -> H + H [10,1000] 1"], strict=True
        )
        h, h2 = net.species["H"], net.species["H2"]
        rxn = Reaction(
            net, [h], [h2], net.reactions[0].rate_coefficient, None, None, "H -> H2", 99
        )

        assert not rxn.check_mass()


class TestWarningNamesReaction:
    def test_warning_includes_reaction_index(self, make_network, caplog):
        # Even non-strict, the warning must identify which reaction failed.
        logging.disable(logging.NOTSET)  # conftest globally disables logging
        try:
            with caplog.at_level(logging.WARNING, logger="JAFF"):
                make_network(_MASS_VIOLATION, strict=False)
        finally:
            logging.disable(logging.CRITICAL)

        mass_msgs = [r.message for r in caplog.records if "Mass not conserved" in r.message]
        assert mass_msgs, "expected a mass-conservation warning"
        # Reaction index (0-based) must appear so the offending reaction is
        # identifiable in a large network.
        assert any("0" in m for m in mass_msgs)
