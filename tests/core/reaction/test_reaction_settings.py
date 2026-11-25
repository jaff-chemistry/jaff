# ABOUTME: A reaction holds its resolved settings as reaction.config (global table with
# ABOUTME: its own entry applied); reaction.photo holds only computed photo state

import pytest

from jaff.core.reaction.config import ReactionConfig

H_PHOTO = "H -> H+ + E  []  PHOTO, 13.6\nH+ + E -> H  []  1e-12\n"
H_KEY = "H._PHOTON__H+.e-"
BANDS = {"bands": [13.6, 20.0, 100.0]}


def _net(make_network, **kwargs):
    return make_network(H_PHOTO, name="h.jet", funcfile=False, **kwargs)


def _photo(net):
    return next(r for r in net.reactions if r.type == "photo")


class TestReactionConfig:
    def test_is_the_resolved_settings(self, make_network):
        net = _net(make_network)
        rxn = _photo(net)

        assert isinstance(rxn.config, ReactionConfig)
        assert rxn.config == net.config.reactions.resolve(H_KEY)

    def test_every_key_is_set(self, make_network):
        rxn = _photo(_net(make_network))

        assert rxn.config.t_cutoff == "clip"
        assert rxn.config.photo is not None
        assert rxn.config.photo.pi_database == "norad"
        assert rxn.config.photo.use_proxy_photoreaction is True

    def test_own_entry_wins_over_global(self, make_network):
        net = _net(
            make_network,
            reactions={
                "t_cutoff": "extrapolate",
                "reaction": {
                    H_KEY: {"t_cutoff": "clip", "photo": {"pi_database": "verner"}}
                },
            },
        )
        rxn = _photo(net)
        other = next(r for r in net.reactions if r.type != "photo")

        assert rxn.config.t_cutoff == "clip"
        assert rxn.config.photo.pi_database == "verner"
        assert other.config.t_cutoff == "extrapolate"

    def test_rate_segments_belong_to_the_reaction(self, make_network):
        rxn = _photo(_net(make_network, reactions={"t_cutoff": "extrapolate"}))

        assert rxn.rate_segments.reaction is rxn
        assert not hasattr(rxn.rate_segments, "mode")

    def test_constructor_takes_no_config(self):
        import inspect

        from jaff.core.reaction import Reaction

        assert "config" not in inspect.signature(Reaction.__init__).parameters


class TestPhotoHoldsComputedStateOnly:
    @pytest.mark.parametrize(
        "name",
        [
            "pi_database",
            "override_custom_rate",
            "override_custom_cross_section",
            "shielding",
        ],
    )
    def test_settings_are_not_copied(self, make_network, name):
        assert not hasattr(_photo(_net(make_network)).photo, name)

    def test_shielding_tables_start_empty(self, make_network):
        assert _photo(_net(make_network)).photo.shielding_tables == {}

    def test_shielding_props_is_gone(self):
        import jaff.core
        import jaff.core.reaction

        assert not hasattr(jaff.core, "ShieldingProps")
        assert not hasattr(jaff.core.reaction, "ShieldingProps")


class TestCustomCrossSectionOverride:
    JFUNC = "@function photoCrossSection0()\n@return 1e-18\n"

    def test_own_override_needs_radiation(self, make_network, tmp_path):
        # Without radiation there is no cross-section rate: neither a global nor
        # a reaction's own override_custom_cross_section has an effect.
        (tmp_path / "h.jfunc").write_text(self.JFUNC)
        net = make_network(
            H_PHOTO,
            name="h.jet",
            reactions={
                "reaction": {H_KEY: {"photo": {"override_custom_cross_section": True}}}
            },
        )

        assert _photo(net).photo.xsec.database == "custom"

    def test_pi_database_with_custom_cross_section_only_if_set_on_reaction(
        self, make_network, tmp_path
    ):
        # The global pi_database default ("norad") must not clash with a
        # photoCrossSection; the reaction's own pi_database does.
        from jaff.errors import ParserError

        (tmp_path / "h.jfunc").write_text(self.JFUNC)
        assert _photo(make_network(H_PHOTO, name="h.jet")).photo.xsec.database == "custom"

        with pytest.raises(ParserError, match="pi_database cannot be set"):
            make_network(
                H_PHOTO,
                name="h.jet",
                reactions={"reaction": {H_KEY: {"photo": {"pi_database": "verner"}}}},
            )
