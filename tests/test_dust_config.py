# ABOUTME: DustConfig: [network.dust] keys, tabulated R_V, reductions and the
# ABOUTME: photoelectric-emission thresholds

import pytest
from pydantic import ValidationError

from jaff.physics.dust.config import DustConfig, PhotoelectricEmissionConfig


class TestDustConfig:
    def test_defaults(self):
        cfg = DustConfig()

        assert cfg.rv == pytest.approx(3.1)
        assert cfg.u_reduction == "absorption"
        assert cfg.f_reduction == "transport"
        assert cfg.pe == PhotoelectricEmissionConfig(
            threshold_low=6.0, threshold_high=13.6
        )

    def test_partial_pe_table_keeps_other_default(self):
        cfg = DustConfig.model_validate({"pe": {"threshold_low": 5.0}})

        assert cfg.pe.threshold_low == pytest.approx(5.0)
        assert cfg.pe.threshold_high == pytest.approx(13.6)

    @pytest.mark.parametrize("rv", [3.1, 4.0, 4, 5.5])
    def test_accepts_tabulated_rv(self, rv):
        assert DustConfig(rv=rv).rv == pytest.approx(rv)

    def test_rejects_untabulated_rv(self):
        with pytest.raises(ValidationError, match="rv"):
            DustConfig.model_validate({"rv": 3.0})

    def test_reductions_are_case_insensitive(self):
        cfg = DustConfig.model_validate(
            {"u_reduction": "None", "f_reduction": "Extinction"}
        )

        assert cfg.u_reduction == "none"
        assert cfg.f_reduction == "extinction"

    @pytest.mark.parametrize("key", ["u_reduction", "f_reduction"])
    def test_rejects_invalid_reduction(self, key):
        with pytest.raises(ValidationError, match=key):
            DustConfig.model_validate({key: "emission"})

    def test_rejects_string_threshold(self):
        with pytest.raises(ValidationError, match="threshold_low"):
            DustConfig.model_validate({"pe": {"threshold_low": "6"}})

    @pytest.mark.parametrize("key", ["pe_threshold_low", "pe_threshold_high"])
    def test_rejects_old_flat_keys(self, key):
        with pytest.raises(ValidationError, match=key):
            DustConfig.model_validate({key: 6.0})
