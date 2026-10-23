# ABOUTME: RadiationConfig: [network.radiation] keys, band edges and the
# ABOUTME: per-band profile index

import pytest
from pydantic import ValidationError

from jaff.physics.photo_reactions.config import RadiationConfig


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
