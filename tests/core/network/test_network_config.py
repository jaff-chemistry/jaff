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


class TestNetworkUsesConfig:
    """Network builds net.config: jaff.toml < config layers < keyword arguments."""

    NET = "H + H -> H2 [10,1000] 1e-10\nH2 -> H + H [10,1000] 1e-10\n"

    def _net(self, make_network, **kwargs):
        return make_network(self.NET, name="n.dat", funcfile=False, **kwargs)

    def test_net_config_comes_from_jaff_toml(self, make_network, tmp_path):
        (tmp_path / "jaff.toml").write_text('[network]\nlabel = "from_toml"\n')

        net = self._net(make_network)

        assert isinstance(net.config, NetworkConfig)
        assert net.config.label == "from_toml"
        assert net.label == "from_toml"

    def test_keyword_beats_jaff_toml(self, make_network, tmp_path):
        (tmp_path / "jaff.toml").write_text('[network]\nlabel = "from_toml"\n')

        assert self._net(make_network, label="kw").label == "kw"

    def test_unset_keyword_keeps_jaff_toml(self, make_network, tmp_path):
        (tmp_path / "jaff.toml").write_text("[network]\nexpand_nuclei = false\n")

        net = self._net(make_network)

        assert net.config.expand_nuclei is False
        assert net.spec.expand_nuclei is False

    def test_use_jaff_toml_false_skips_it(self, make_network, tmp_path):
        (tmp_path / "jaff.toml").write_text('[network]\nlabel = "from_toml"\n')

        assert self._net(make_network, use_jaff_toml=False).label == "n"

    def test_config_dict_layer(self, make_network):
        net = self._net(make_network, config={"network": {"label": "from_dict"}})

        assert net.label == "from_dict"

    def test_strict_keyword(self, make_network):
        assert self._net(make_network, strict=True).config.strict is True

    def test_old_jaff_toml_key_is_reported(self, make_network, tmp_path):
        (tmp_path / "jaff.toml").write_text('[network.rates]\nT_cutoff = "clip"\n')

        with pytest.raises(ParserError, match=r"network\.rates: unknown key"):
            self._net(make_network)

    def test_photo_settings_on_non_photo_reaction(self, make_network):
        config = {
            "network": {
                "reactions": {
                    "reaction": {"H.H__H2": {"photo": {"pi_database": "verner"}}}
                }
            }
        }

        with pytest.raises(ParserError, match="only valid for photo reactions"):
            self._net(make_network, config=config)

    def test_settings_for_unknown_reaction_are_reported(self, make_network):
        config = {
            "network": {
                "reactions": {
                    "reaction": {
                        "H.H__H2": {"t_cutoff": "extrapolate"},
                        "He.H__HeH": {"t_cutoff": "clip"},
                        "C.O__CO": {"t_cutoff": "clip"},
                    }
                }
            }
        }

        with pytest.raises(ParserError) as info:
            self._net(make_network, config=config)

        message = str(info.value)
        assert 'network.reactions.reaction."He.H__HeH"' in message
        assert 'network.reactions.reaction."C.O__CO"' in message
        assert "H.H__H2" not in message

    def test_settings_for_known_reactions_pass(self, make_network):
        config = {
            "network": {"reactions": {"reaction": {"H.H__H2": {"t_cutoff": "clip"}}}}
        }

        assert self._net(make_network, config=config).config.reactions.reaction
