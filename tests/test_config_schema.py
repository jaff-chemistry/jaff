# ABOUTME: Tests for the pydantic config schema models in jaff.config.schema
# ABOUTME: Defaults, strict key/value validation and per-reaction resolution

import pytest
from pydantic import ValidationError

from jaff.config.schema import (
    PhotoConfig,
    ReactionConfig,
    RadiationConfig,
    ReactionsConfig,
    ShieldingConfig,
)

H2_KEY = "H2._PHOTON__H.H"
CO_KEY = "CO._PHOTON__C.O"


class TestReactionsConfig:
    def test_defaults(self):
        cfg = ReactionsConfig()

        assert cfg.duplicate_policy == "preserve-first"
        assert cfg.t_cutoff == "clip"
        assert cfg.photo == PhotoConfig(
            pi_database="norad",
            override_custom_rate=False,
            override_custom_cross_section=False,
            use_proxy_photoreaction=True,
        )
        assert cfg.reaction == {}

    def test_rejects_unknown_key(self):
        with pytest.raises(ValidationError, match="T_cutoff"):
            ReactionsConfig.model_validate({"T_cutoff": "clip"})

    def test_rejects_unknown_photo_key(self):
        with pytest.raises(ValidationError, match="pi_db"):
            ReactionsConfig.model_validate({"photo": {"pi_db": "norad"}})

    @pytest.mark.parametrize(
        "key, value",
        [
            ("t_cutoff", "clamp"),
            ("duplicate_policy", "keep"),
        ],
    )
    def test_rejects_invalid_value(self, key, value):
        with pytest.raises(ValidationError, match=key):
            ReactionsConfig.model_validate({key: value})

    def test_rejects_string_for_bool(self):
        with pytest.raises(ValidationError, match="override_custom_rate"):
            ReactionsConfig.model_validate({"photo": {"override_custom_rate": "true"}})

    def test_rejects_invalid_pi_database(self):
        with pytest.raises(ValidationError, match="pi_database"):
            ReactionsConfig.model_validate({"photo": {"pi_database": "nist"}})

    def test_choice_values_are_case_insensitive(self):
        cfg = ReactionsConfig.model_validate(
            {
                "t_cutoff": "EXTRAPOLATE",
                "duplicate_policy": "Preserve-Last",
                "photo": {"pi_database": "Verner"},
            }
        )

        assert cfg.t_cutoff == "extrapolate"
        assert cfg.duplicate_policy == "preserve-last"
        assert cfg.photo.pi_database == "verner"

    def test_is_frozen(self):
        cfg = ReactionsConfig()

        with pytest.raises(ValidationError):
            cfg.t_cutoff = "extrapolate"

    def test_per_reaction_entry_rejects_nested_reaction_table(self):
        with pytest.raises(ValidationError, match="reaction"):
            ReactionsConfig.model_validate({"reaction": {"H+.e-__H": {"reaction": {}}}})


class TestReactionConfig:
    def test_defaults_are_unset(self):
        cfg = ReactionConfig()

        assert cfg.duplicate_policy is None
        assert cfg.t_cutoff is None
        assert cfg.photo is None

    def test_accepts_every_global_key(self):
        cfg = ReactionConfig.model_validate(
            {
                "duplicate_policy": "error",
                "t_cutoff": "extrapolate",
                "photo": {
                    "pi_database": "leiden",
                    "override_custom_rate": True,
                    "override_custom_cross_section": True,
                    "use_proxy_photoreaction": False,
                },
            }
        )

        assert cfg.duplicate_policy == "error"
        assert cfg.photo is not None
        assert cfg.photo.use_proxy_photoreaction is False


class TestResolve:
    def test_without_entry_returns_global_settings(self):
        cfg = ReactionsConfig(t_cutoff="extrapolate")

        resolved = cfg.resolve("H+.e-__H")

        assert resolved.t_cutoff == "extrapolate"
        assert resolved.duplicate_policy == "preserve-first"
        assert resolved.photo == cfg.photo

    def test_entry_overrides_only_the_keys_it_sets(self):
        cfg = ReactionsConfig.model_validate(
            {
                "t_cutoff": "extrapolate",
                "duplicate_policy": "error",
                "reaction": {"H+.e-__H": {"t_cutoff": "clip"}},
            }
        )

        resolved = cfg.resolve("H+.e-__H")

        assert resolved.t_cutoff == "clip"
        assert resolved.duplicate_policy == "error"

    def test_photo_entry_overrides_field_by_field(self):
        cfg = ReactionsConfig.model_validate(
            {
                "photo": {"pi_database": "verner", "override_custom_rate": True},
                "reaction": {"H.PHOTON__H+.e-": {"photo": {"pi_database": "leiden"}}},
            }
        )

        resolved = cfg.resolve("H.PHOTON__H+.e-")

        assert resolved.photo is not None
        assert resolved.photo.pi_database == "leiden"
        assert resolved.photo.override_custom_rate is True
        assert resolved.photo.use_proxy_photoreaction is True

    def test_resolved_config_has_no_unset_keys(self):
        # shielding is optional: None means the reaction is not shielded.
        resolved = ReactionsConfig().resolve("H+.e-__H")

        assert resolved.duplicate_policy is not None
        assert resolved.t_cutoff is not None
        assert resolved.photo is not None
        assert None not in resolved.photo.model_dump(exclude={"shielding"}).values()


def _resolved_shielding(cfg: ReactionsConfig, serialized: str):
    photo = cfg.resolve(serialized).photo
    assert photo is not None

    return photo.shielding


class TestShieldingConfig:
    def test_flat_table_splits_type_from_params(self):
        cfg = ShieldingConfig.model_validate({"type": "Leiden", "shielded_by": ["self"]})

        assert cfg.type == "leiden"
        assert cfg.params == {"shielded_by": ["self"]}

    def test_type_defaults_to_leiden(self):
        assert ShieldingConfig.model_validate({"shielded_by": ["self"]}).type == "leiden"

    def test_global_params_validated_by_plugin(self):
        with pytest.raises(ValidationError, match="shielded_by"):
            ReactionsConfig.model_validate({"photo": {"shielding": {"type": "leiden"}}})

    def test_global_rejects_reaction_bound_type(self):
        with pytest.raises(ValidationError, match="hg2015"):
            ReactionsConfig.model_validate({"photo": {"shielding": {"type": "hg2015"}}})

    def test_rejects_unknown_type(self):
        with pytest.raises(ValidationError, match="nope"):
            ReactionsConfig.model_validate({"photo": {"shielding": {"type": "nope"}}})

    def test_per_reaction_accepts_reaction_bound_type(self):
        cfg = ReactionsConfig.model_validate(
            {"reaction": {H2_KEY: {"photo": {"shielding": {"type": "hg2015"}}}}}
        )

        shielding = _resolved_shielding(cfg, H2_KEY)

        assert shielding is not None
        assert shielding.type == "hg2015"

    def test_per_reaction_rejects_type_bound_to_another_reaction(self):
        with pytest.raises(ValidationError, match="CO._PHOTON__C.O"):
            ReactionsConfig.model_validate(
                {"reaction": {CO_KEY: {"photo": {"shielding": {"type": "hg2015"}}}}}
            )

    def test_per_reaction_unknown_param_names_key_and_reaction(self):
        with pytest.raises(ValidationError, match=r"(?s)H2\._PHOTON__H\.H.*foo"):
            ReactionsConfig.model_validate(
                {
                    "reaction": {
                        H2_KEY: {"photo": {"shielding": {"type": "hg2015", "foo": 1}}}
                    }
                }
            )

    def test_per_reaction_shielding_replaces_global(self):
        cfg = ReactionsConfig.model_validate(
            {
                "photo": {
                    "shielding": {"shielded_by": ["self", "H2"], "radiation": "draine"}
                },
                "reaction": {CO_KEY: {"photo": {"shielding": {"shielded_by": ["self"]}}}},
            }
        )

        shielding = _resolved_shielding(cfg, CO_KEY)

        assert shielding is not None
        assert shielding.params == {"shielded_by": ["self"]}

    def test_global_shielding_applies_without_entry(self):
        cfg = ReactionsConfig.model_validate(
            {"photo": {"shielding": {"shielded_by": ["self"]}}}
        )

        assert _resolved_shielding(cfg, CO_KEY) == cfg.photo.shielding


class TestRadiationConfig:
    def test_defaults(self):
        cfg = RadiationConfig(bands=[13.6, "inf"])

        assert cfg.profile_index == 0
        assert cfg.mode == "nph"
        assert cfg.c == pytest.approx(2.99792458e10)
        assert cfg.background_field == "draine"
        assert cfg.absorption_attenuates_flux is True

    def test_bands_are_required(self):
        with pytest.raises(ValidationError, match="bands"):
            RadiationConfig.model_validate({})

    def test_accepts_ints_and_trailing_inf(self):
        cfg = RadiationConfig.model_validate({"bands": [6, 13.6, "inf"]})

        assert cfg.bands == [6, 13.6, "inf"]

    @pytest.mark.parametrize(
        "bands",
        [[13.6], ["inf", 13.6], [6.0, "inf", 100.0], [6.0, "big"]],
        ids=["one-edge", "inf-first", "inf-middle", "other-string"],
    )
    def test_rejects_invalid_bands(self, bands):
        with pytest.raises(ValidationError, match="bands"):
            RadiationConfig.model_validate({"bands": bands})

    def test_accepts_one_profile_index_per_band(self):
        cfg = RadiationConfig.model_validate(
            {"bands": [6.0, 13.6, 100.0], "profile_index": [0, 1.5]}
        )

        assert cfg.profile_index == [0, 1.5]

    def test_rejects_profile_index_length_mismatch(self):
        with pytest.raises(ValidationError, match="profile_index"):
            RadiationConfig.model_validate(
                {"bands": [6.0, 13.6, 100.0], "profile_index": [0, 1, 2]}
            )

    def test_rejects_bool_profile_index(self):
        with pytest.raises(ValidationError, match="profile_index"):
            RadiationConfig.model_validate({"bands": [6.0, 13.6], "profile_index": True})

    def test_choices_are_case_insensitive(self):
        cfg = RadiationConfig.model_validate(
            {"bands": [6.0, 13.6], "mode": "U", "background_field": "Habing"}
        )

        assert cfg.mode == "u"
        assert cfg.background_field == "habing"

    @pytest.mark.parametrize(
        "key, value", [("mode", "energy"), ("background_field", "isrf")]
    )
    def test_rejects_invalid_choice(self, key, value):
        with pytest.raises(ValidationError, match=key):
            RadiationConfig.model_validate({"bands": [6.0, 13.6], key: value})

    @pytest.mark.parametrize("c", ["c_hat", 3e9, 30000000000])
    def test_c_is_number_or_symbol_name(self, c):
        assert RadiationConfig(bands=[6.0, 13.6], c=c).c == c

    @pytest.mark.parametrize("key", ["rsl", "pi_database", "use_proxy_photoreaction"])
    def test_rejects_old_and_moved_keys(self, key):
        with pytest.raises(ValidationError, match=key):
            RadiationConfig.model_validate({"bands": [6.0, 13.6], key: "x"})
