# ABOUTME: ReactionRadiation: every reaction's radiation component (rxn.rad): its
# ABOUTME: spectral yield, read from the .jfunc, and the yield per radiation group

import pytest

from jaff.core.reaction import ReactionRadiation

RECOMB = "H+ + E -> H  []  2.5e-13\n"
FLAT = "@function spectralYield0()\n@return 1e-12\n"
BANDS = {"bands": [6.0, 13.6, 100.0]}


@pytest.fixture
def emitting(make_network, tmp_path):
    def _make(**kwargs):
        (tmp_path / "y.jfunc").write_text(FLAT)

        return make_network(RECOMB, name="y.jet", **kwargs).reactions[0]

    return _make


class TestComponent:
    def test_every_reaction_has_one(self, make_network):
        rxn = make_network(RECOMB, name="y.jet", funcfile=False).reactions[0]

        assert isinstance(rxn.rad, ReactionRadiation)
        assert rxn.rad.reaction is rxn

    def test_yield_is_the_reaction_function_not_a_copy(self, emitting):
        rxn = emitting()

        assert rxn.rad.spectral_yield is rxn.function("spectralYield")

    def test_no_yield_without_the_function(self, make_network):
        rxn = make_network(RECOMB, name="y.jet", funcfile=False).reactions[0]

        assert rxn.rad.spectral_yield is None
        assert rxn.rad.group_yields is None

    def test_group_yields_need_radiation(self, emitting):
        assert emitting().rad.group_yields is None

    def test_group_yields_one_per_group(self, emitting):
        rxn = emitting(radiation=BANDS)

        assert len(rxn.rad.group_yields) == 2


def test_reaction_radiation_props_is_gone():
    import jaff.core
    import jaff.core.reaction

    assert not hasattr(jaff.core.reaction, "ReactionRadiationProps")
    assert jaff.core.ReactionRadiation is ReactionRadiation
