# ABOUTME: NetworkConfig: [network] top-level keys and its sections
# ABOUTME: (reactions, radiation, dust, thermodynamics)

from pathlib import Path

import pytest
from pydantic import ValidationError

from jaff.core.network.config import NetworkConfig
from jaff.core.reaction.config import ReactionsConfig
from jaff.physics.dust.config import DustConfig
from jaff.physics.thermodynamics.config import ThermodynamicsConfig


class TestNetworkConfig:
    def test_defaults(self):
        cfg = NetworkConfig()

        assert cfg.label is None
        assert cfg.strict is False
        assert cfg.funcfile is True
        assert cfg.expand_nuclei is True
        assert cfg.reactions == ReactionsConfig()
        assert cfg.radiation is None
        assert cfg.dust is None
        assert cfg.thermodynamics == ThermodynamicsConfig()

    def test_sections_are_validated(self):
        cfg = NetworkConfig.model_validate(
            {
                "reactions": {"t_cutoff": "extrapolate"},
                "radiation": {"bands": [6.0, 13.6]},
                "thermodynamics": {"eos": {"gamma": 1.4}},
            }
        )

        assert cfg.reactions.t_cutoff == "extrapolate"
        assert cfg.radiation is not None
        assert cfg.radiation.bands == [6.0, 13.6]
        assert cfg.thermodynamics.eos.type == "ideal"

    def test_empty_dust_table_enables_dust_with_defaults(self):
        assert NetworkConfig.model_validate({"dust": {}}).dust == DustConfig()

    @pytest.mark.parametrize("funcfile", [True, False])
    def test_funcfile_accepts_bool(self, funcfile):
        assert NetworkConfig(funcfile=funcfile).funcfile is funcfile

    @pytest.mark.parametrize("funcfile", ["net.jfunc", Path("net.jfunc")])
    def test_funcfile_path_becomes_path(self, funcfile):
        assert NetworkConfig.model_validate({"funcfile": funcfile}).funcfile == Path(
            "net.jfunc"
        )

    def test_rejects_string_for_bool(self):
        with pytest.raises(ValidationError, match="strict"):
            NetworkConfig.model_validate({"strict": "true"})

    @pytest.mark.parametrize(
        "key", ["errors", "config", "duplicate_policy", "rates", "eos"]
    )
    def test_rejects_old_and_moved_keys(self, key):
        with pytest.raises(ValidationError, match=key):
            NetworkConfig.model_validate({key: {}})
