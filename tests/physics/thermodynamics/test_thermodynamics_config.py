# ABOUTME: ThermodynamicsConfig: the EOS model selected by its type and the
# ABOUTME: adiabatic-index bounds

import pytest
from pydantic import ValidationError

from jaff.physics.thermodynamics.config import ThermodynamicsConfig
from jaff.physics.thermodynamics.eos.config import (
    FermiDegenerateEosConfig,
    IdealEosConfig,
    MultiGammaEosConfig,
    RelativisticFermiDegenerateEosConfig,
)


class TestThermodynamicsConfig:
    def test_default_eos_is_ideal(self):
        eos = ThermodynamicsConfig().eos

        assert isinstance(eos, IdealEosConfig)
        assert eos.gamma == pytest.approx(5.0 / 3.0)

    def test_type_defaults_to_ideal(self):
        eos = ThermodynamicsConfig.model_validate({"eos": {"gamma": 1.4}}).eos

        assert isinstance(eos, IdealEosConfig)
        assert eos.gamma == pytest.approx(1.4)

    @pytest.mark.parametrize(
        "table, cls",
        [
            ({"type": "Ideal"}, IdealEosConfig),
            (
                {"type": "multi_gamma", "default_gamma": 1.4, "gamma_map": {"H2": 1.4}},
                MultiGammaEosConfig,
            ),
            ({"type": "fermi_degenerate"}, FermiDegenerateEosConfig),
            (
                {"type": "relativistic_fermi_degenerate"},
                RelativisticFermiDegenerateEosConfig,
            ),
        ],
    )
    def test_type_selects_model(self, table, cls):
        assert isinstance(ThermodynamicsConfig.model_validate({"eos": table}).eos, cls)

    def test_unknown_type_error_lists_choices(self):
        with pytest.raises(ValidationError) as info:
            ThermodynamicsConfig.model_validate({"eos": {"type": "Polytrope"}})

        message = str(info.value)
        assert "unknown EOS type 'polytrope'" in message
        assert (
            "ideal, multi_gamma, fermi_degenerate, relativistic_fermi_degenerate"
            in message
        )
        assert "_eos_type" not in message

    def test_rejects_unknown_type(self):
        with pytest.raises(ValidationError, match="polytrope"):
            ThermodynamicsConfig.model_validate({"eos": {"type": "polytrope"}})

    @pytest.mark.parametrize("missing", ["default_gamma", "gamma_map"])
    def test_multi_gamma_requires_both_keys(self, missing):
        table = {"type": "multi_gamma", "default_gamma": 1.4, "gamma_map": {"H2": 1.4}}
        del table[missing]

        with pytest.raises(ValidationError, match=missing):
            ThermodynamicsConfig.model_validate({"eos": table})

    @pytest.mark.parametrize(
        "table, key",
        [
            ({"type": "ideal", "gamma": 1.0}, "gamma"),
            (
                {"type": "multi_gamma", "default_gamma": 0.9, "gamma_map": {}},
                "default_gamma",
            ),
            (
                {"type": "multi_gamma", "default_gamma": 1.4, "gamma_map": {"H2": 1.0}},
                "gamma_map",
            ),
        ],
    )
    def test_adiabatic_index_must_exceed_one(self, table, key):
        with pytest.raises(ValidationError, match=key):
            ThermodynamicsConfig.model_validate({"eos": table})

    def test_rejects_key_of_another_type(self):
        with pytest.raises(ValidationError, match="gamma"):
            ThermodynamicsConfig.model_validate(
                {"eos": {"type": "fermi_degenerate", "gamma": 1.4}}
            )

    def test_empty_gamma_map_accepted(self):
        eos = ThermodynamicsConfig.model_validate(
            {"eos": {"type": "multi_gamma", "default_gamma": 1.4, "gamma_map": {}}}
        ).eos

        assert isinstance(eos, MultiGammaEosConfig)
        assert eos.gamma_map == {}

    def test_rejects_string_gamma(self):
        with pytest.raises(ValidationError, match="gamma"):
            ThermodynamicsConfig.model_validate({"eos": {"gamma": "5/3"}})

    def test_rejects_non_table_gamma_map(self):
        with pytest.raises(ValidationError, match="gamma_map"):
            ThermodynamicsConfig.model_validate(
                {"eos": {"type": "multi_gamma", "default_gamma": 1.4, "gamma_map": [1.4]}}
            )

    def test_rejects_bool_gamma(self):
        with pytest.raises(ValidationError, match="gamma"):
            ThermodynamicsConfig.model_validate({"eos": {"gamma": True}})
