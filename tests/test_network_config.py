# ABOUTME: NetworkConfig: [network] keys and sections (reactions, radiation, dust, ...),
# ABOUTME: and load_network_config (jaff.toml < config layers < overrides)

from pathlib import Path

import pytest
from pydantic import ValidationError

from jaff.config import NETWORKS_DIR
from jaff.core.network.config import (
    JaffTomlConfig,
    NetworkConfig,
    jaff_toml_for,
    load_network_config,
)
from jaff.core.reaction.config import ReactionsConfig
from jaff.errors import ParserError
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


class TestJaffTomlFor:
    def test_next_to_network_file(self, tmp_path):
        net = tmp_path / "net.jet"
        net.write_text("")

        assert jaff_toml_for(net) == tmp_path / "jaff.toml"

    def test_predefined_network_name(self):
        assert jaff_toml_for(Path("GOW")) == NETWORKS_DIR / "GOW" / "jaff.toml"


@pytest.fixture
def network(tmp_path):
    """A network file in tmp_path/net; write(name, text) writes under tmp_path."""
    (tmp_path / "net").mkdir()
    net = tmp_path / "net" / "net.jet"
    net.write_text("")

    def write(name, text):
        path = tmp_path / name
        path.write_text(text)

        return path

    return net, write


class TestLoadNetworkConfig:
    def test_nothing_gives_defaults(self, network):
        net, _ = network

        assert load_network_config(net) == NetworkConfig()

    def test_jaff_toml_next_to_network_is_applied(self, network):
        net, write = network
        write("net/jaff.toml", '[network]\nlabel = "auto"\n')

        assert load_network_config(net).label == "auto"

    def test_use_jaff_toml_false_skips_it(self, network):
        net, write = network
        write("net/jaff.toml", '[network]\nlabel = "auto"\n')

        assert load_network_config(net, use_jaff_toml=False).label is None

    def test_config_file_layers_over_jaff_toml(self, network):
        net, write = network
        write("net/jaff.toml", '[network]\nlabel = "auto"\nstrict = true\n')
        write("extra.toml", '[network]\nlabel = "extra"\n')

        cfg = load_network_config(net, "extra.toml", cwd=net.parent.parent)

        assert cfg.label == "extra"
        assert cfg.strict is True

    def test_dict_layer_paths_resolve_against_cwd(self, network):
        net, _ = network

        cfg = load_network_config(
            net, {"network": {"funcfile": "aux.jfunc"}}, cwd=net.parent
        )

        assert cfg.funcfile == (net.parent / "aux.jfunc").resolve()

    def test_model_layer_only_overrides_keys_it_sets(self, network):
        net, write = network
        write("net/jaff.toml", "[network]\nstrict = true\n")

        cfg = load_network_config(net, JaffTomlConfig(network=NetworkConfig(label="m")))

        assert cfg.label == "m"
        assert cfg.strict is True

    def test_layers_apply_in_order_then_overrides(self, network):
        net, _ = network

        cfg = load_network_config(
            net,
            [{"network": {"label": "a", "strict": True}}, {"network": {"label": "b"}}],
            overrides={"network": {"strict": False}},
        )

        assert cfg.label == "b"
        assert cfg.strict is False

    def test_model_layer_with_shielding_round_trips(self, network):
        net, _ = network
        model = JaffTomlConfig.model_validate(
            {
                "network": {
                    "reactions": {
                        "photo": {
                            "shielding": {"type": "leiden", "shielded_by": ["self"]}
                        }
                    }
                }
            }
        )

        cfg = load_network_config(net, model)

        assert cfg.reactions.photo.shielding == model.network.reactions.photo.shielding

    def test_invalid_layer_is_parser_error(self, network):
        net, _ = network

        with pytest.raises(ParserError, match=r"network\.rates: unknown key"):
            load_network_config(net, {"network": {"rates": {}}})
