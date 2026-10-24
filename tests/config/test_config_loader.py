# ABOUTME: Config loading: TOML layers with paths resolved per file, and deep_merge (later wins,
# ABOUTME: tables merge, lists/scalars replace, _atomic tables replace whole)

import copy
import tomllib
from pathlib import Path

import pytest

from jaff.cli.jaffgen.config import JaffgenTomlConfig
from jaff.config.loader import (
    atomic_tables,
    deep_merge,
    dump_config,
    load_config,
    load_layer,
    resolve_paths,
)
from jaff.core.network.config import JaffTomlConfig
from jaff.errors import ParserError


class TestDeepMerge:
    def test_no_layers_is_empty(self):
        assert deep_merge() == {}

    def test_later_scalar_wins(self):
        assert deep_merge({"a": 1}, {"a": 2}) == {"a": 2}

    def test_disjoint_keys_are_kept(self):
        assert deep_merge({"a": 1}, {"b": 2}) == {"a": 1, "b": 2}

    def test_tables_merge_recursively(self):
        low = {"network": {"label": "x", "reactions": {"t_cutoff": "clip"}}}
        high = {"network": {"reactions": {"duplicate_policy": "error"}}}

        assert deep_merge(low, high) == {
            "network": {
                "label": "x",
                "reactions": {"t_cutoff": "clip", "duplicate_policy": "error"},
            }
        }

    def test_lists_are_replaced(self):
        assert deep_merge({"bands": [6.0, 13.6, 100.0]}, {"bands": [13.6, "inf"]}) == {
            "bands": [13.6, "inf"]
        }

    def test_scalar_replaces_table_and_table_replaces_scalar(self):
        assert deep_merge({"f": {"a": 1}}, {"f": True}) == {"f": True}
        assert deep_merge({"f": True}, {"f": {"a": 1}}) == {"f": {"a": 1}}

    def test_three_layers_apply_in_order(self):
        assert deep_merge({"a": 1, "b": 1}, {"a": 2}, {"a": 3}) == {"a": 3, "b": 1}

    def test_inputs_are_not_modified(self):
        low = {"t": {"a": [1], "b": {"c": 1}}}
        high = {"t": {"b": {"d": 2}}}
        before = copy.deepcopy((low, high))

        merged = deep_merge(low, high)
        merged["t"]["a"].append(2)
        merged["t"]["b"]["e"] = 3

        assert (low, high) == before

    def test_atomic_table_is_replaced_whole(self):
        low = {"photo": {"shielding": {"type": "leiden", "shielded_by": ["self"]}}}
        high = {"photo": {"shielding": {"type": "hg2015"}}}

        merged = deep_merge(low, high, atomic=lambda path: path[-1] == "shielding")

        assert merged == {"photo": {"shielding": {"type": "hg2015"}}}

    def test_atomic_receives_full_key_path(self):
        seen = []

        def atomic(path):
            seen.append(path)
            return False

        deep_merge(
            {"network": {"dust": {"rv": 3.1}}}, {"network": {"dust": {}}}, atomic=atomic
        )

        assert seen == [("network",), ("network", "dust")]


H2_SHIELDING = (
    "network",
    "reactions",
    "reaction",
    "H2._PHOTON__H.H",
    "photo",
    "shielding",
)


class TestAtomicTables:
    @pytest.mark.parametrize(
        "path",
        [
            ("network", "reactions", "photo", "shielding"),
            H2_SHIELDING,
            ("network", "thermodynamics", "eos"),
        ],
    )
    def test_tables_selected_by_type_are_atomic(self, path):
        assert atomic_tables(JaffgenTomlConfig)(path) is True

    @pytest.mark.parametrize(
        "path",
        [
            ("network",),
            ("network", "reactions"),
            ("network", "reactions", "photo"),
            ("network", "reactions", "reaction", "H2._PHOTON__H.H"),
            ("network", "dust"),
            ("network", "dust", "pe"),
            ("network", "radiation"),
            ("network", "thermodynamics"),
            ("jaffgen",),
            ("network", "not_a_key"),
            ("network", "label", "deeper"),
        ],
    )
    def test_other_tables_merge(self, path):
        assert atomic_tables(JaffgenTomlConfig)(path) is False

    def test_works_from_jaff_toml_root(self):
        assert atomic_tables(JaffTomlConfig)(H2_SHIELDING) is True

    def test_eos_type_switch_across_layers(self):
        jaff_toml = {
            "network": {
                "thermodynamics": {
                    "eos": {"type": "multi_gamma", "default_gamma": 1.4, "gamma_map": {}}
                }
            }
        }
        jaffgen_toml = {"network": {"thermodynamics": {"eos": {"gamma": 1.3}}}}

        merged = deep_merge(jaff_toml, jaffgen_toml, atomic=atomic_tables(JaffTomlConfig))

        assert merged["network"]["thermodynamics"]["eos"] == {"gamma": 1.3}


class TestResolvePaths:
    def test_relative_path_resolved_against_base(self, tmp_path):
        layer = resolve_paths(
            {"network": {"funcfile": "aux/net.jfunc"}}, JaffTomlConfig, tmp_path
        )

        assert layer["network"]["funcfile"] == (tmp_path / "aux/net.jfunc").resolve()

    def test_absolute_path_unchanged(self, tmp_path):
        target = tmp_path / "net.jfunc"
        layer = resolve_paths(
            {"network": {"funcfile": str(target)}}, JaffTomlConfig, Path("/")
        )

        assert layer["network"]["funcfile"] == target

    def test_bool_funcfile_unchanged(self, tmp_path):
        layer = resolve_paths({"network": {"funcfile": False}}, JaffTomlConfig, tmp_path)

        assert layer["network"]["funcfile"] is False

    def test_jaffgen_paths_resolved(self, tmp_path):
        layer = resolve_paths(
            {
                "jaffgen": {
                    "input_dir": "tpl",
                    "input_files": ["a.cpp", "sub/b.f90"],
                    "output_dir": "../out",
                    "network_file": "GOW",
                    "network_config": "net/jaff.toml",
                }
            },
            JaffgenTomlConfig,
            tmp_path,
        )

        jg = layer["jaffgen"]
        assert jg["input_dir"] == (tmp_path / "tpl").resolve()
        assert jg["input_files"] == [
            (tmp_path / "a.cpp").resolve(),
            (tmp_path / "sub/b.f90").resolve(),
        ]
        assert jg["output_dir"] == (tmp_path / "../out").resolve()
        assert jg["network_file"] == (tmp_path / "GOW").resolve()
        assert jg["network_config"] == (tmp_path / "net/jaff.toml").resolve()

    def test_non_path_strings_unchanged(self, tmp_path):
        layer = {"jaffgen": {"template": "microphysics"}, "network": {"label": "a/b"}}

        assert resolve_paths(layer, JaffgenTomlConfig, tmp_path) == layer

    def test_unknown_keys_and_tables_unchanged(self, tmp_path):
        layer = {
            "network": {"rsl": "c_hat", "funcfile": 5},
            "jaffgen": {"table": [{"source": {"path": "data.csv"}}]},
        }

        assert resolve_paths(layer, JaffgenTomlConfig, tmp_path) == layer

    def test_input_is_not_modified(self, tmp_path):
        layer = {"network": {"funcfile": "net.jfunc"}}

        resolve_paths(layer, JaffTomlConfig, tmp_path)

        assert layer == {"network": {"funcfile": "net.jfunc"}}


class TestLoadLayer:
    def test_reads_toml_and_resolves_against_its_directory(self, tmp_path):
        cfg = tmp_path / "net" / "jaff.toml"
        cfg.parent.mkdir()
        cfg.write_text(
            '[network]\nfuncfile = "x.jfunc"\n\n[network.reactions]\nt_cutoff = "clip"\n'
        )

        layer = load_layer(cfg, JaffTomlConfig)

        assert layer == {
            "network": {
                "funcfile": (cfg.parent / "x.jfunc").resolve(),
                "reactions": {"t_cutoff": "clip"},
            }
        }

    def test_missing_file(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_layer(tmp_path / "jaff.toml", JaffTomlConfig)

    def test_invalid_toml_names_the_file(self, tmp_path):
        cfg = tmp_path / "jaff.toml"
        cfg.write_text("[network\n")

        with pytest.raises(ParserError, match="jaff.toml"):
            load_layer(cfg, JaffTomlConfig)


class TestLoadConfig:
    def test_no_layers_gives_defaults(self):
        assert load_config([], JaffTomlConfig) == JaffTomlConfig()

    def test_later_layer_wins(self):
        cfg = load_config(
            [
                {"network": {"label": "low", "strict": True}},
                {"network": {"label": "high"}},
            ],
            JaffTomlConfig,
        )

        assert cfg.network.label == "high"
        assert cfg.network.strict is True

    def test_atomic_tables_replace_across_layers(self):
        cfg = load_config(
            [
                {
                    "network": {
                        "thermodynamics": {
                            "eos": {
                                "type": "multi_gamma",
                                "default_gamma": 1.4,
                                "gamma_map": {},
                            }
                        }
                    }
                },
                {"network": {"thermodynamics": {"eos": {"gamma": 1.3}}}},
            ],
            JaffTomlConfig,
        )

        assert cfg.network.thermodynamics.eos.type == "ideal"

    def test_unknown_key_reported_by_path(self):
        with pytest.raises(ParserError, match=r"network\.radiation\.rsl: unknown key"):
            load_config(
                [{"network": {"radiation": {"bands": [6.0, 13.6], "rsl": 1}}}],
                JaffTomlConfig,
            )

    def test_missing_key_reported_by_path(self):
        with pytest.raises(
            ParserError, match=r"network\.radiation\.bands: required key missing"
        ):
            load_config([{"network": {"radiation": {}}}], JaffTomlConfig)

    def test_every_problem_is_listed(self):
        with pytest.raises(ParserError) as info:
            load_config(
                [{"network": {"strict": "yes", "reactions": {"t_cutoff": "clamp"}}}],
                JaffTomlConfig,
            )

        message = str(info.value)
        assert "network.strict:" in message
        assert "network.reactions.t_cutoff:" in message

    def test_reaction_key_is_quoted(self):
        with pytest.raises(ParserError, match=r'reaction\."H\+\.e-__H"\.t_cutoff: '):
            load_config(
                [
                    {
                        "network": {
                            "reactions": {"reaction": {"H+.e-__H": {"t_cutoff": "x"}}}
                        }
                    }
                ],
                JaffTomlConfig,
            )

    def test_validator_message_has_no_pydantic_prefix(self):
        with pytest.raises(ParserError) as info:
            load_config(
                [
                    {
                        "network": {
                            "radiation": {"bands": [6.0, 13.6], "profile_index": [0, 1]}
                        }
                    }
                ],
                JaffTomlConfig,
            )

        message = str(info.value)
        assert "network.radiation: profile_index has 2 entries" in message
        assert "Value error" not in message

    def test_shielding_plugin_errors_use_the_same_format(self):
        layer = {
            "network": {
                "reactions": {
                    "reaction": {
                        "H2._PHOTON__H.H": {
                            "photo": {"shielding": {"type": "hg2015", "min_ncl": 1.0}}
                        }
                    }
                }
            }
        }

        with pytest.raises(ParserError) as info:
            load_config([layer], JaffTomlConfig)

        message = str(info.value)
        assert (
            'network.reactions.reaction."H2._PHOTON__H.H".photo.shielding.min_ncl: unknown key'
            in message
        )
        assert "errors.pydantic.dev" not in message


class TestLoadLayerSections:
    def test_unknown_section_names_the_file(self, tmp_path):
        cfg = tmp_path / "jaff.toml"
        cfg.write_text('[jaffgen]\ntemplate = "microphysics"\n')

        with pytest.raises(ParserError, match=r"jaff\.toml.*jaffgen"):
            load_layer(cfg, JaffTomlConfig)


class TestDumpConfig:
    def _config(self):
        return JaffgenTomlConfig.model_validate(
            {
                "jaffgen": {
                    "output_dir": Path("/abs/out"),
                    "input_files": [Path("/abs/a.cpp")],
                },
                "network": {
                    "radiation": {"bands": [13.6, "inf"], "c": "c_hat"},
                    "reactions": {
                        "reaction": {
                            "CO._PHOTON__C.O": {
                                "photo": {"shielding": {"shielded_by": ["self"]}}
                            }
                        }
                    },
                    "thermodynamics": {"eos": {"gamma": 1.4}},
                },
            }
        )

    def test_round_trips(self):
        config = self._config()

        assert (
            JaffgenTomlConfig.model_validate(tomllib.loads(dump_config(config))) == config
        )

    def test_holds_defaults(self):
        dumped = tomllib.loads(dump_config(self._config()))

        assert dumped["network"]["strict"] is False
        assert dumped["network"]["reactions"]["t_cutoff"] == "clip"
        assert dumped["network"]["thermodynamics"]["eos"]["type"] == "ideal"

    def test_unset_optionals_are_omitted(self):
        dumped = tomllib.loads(dump_config(self._config()))

        assert "label" not in dumped["network"]
        assert "dust" not in dumped["network"]
        assert "template" not in dumped["jaffgen"]

    def test_paths_and_shielding_are_written_as_in_toml(self):
        dumped = tomllib.loads(dump_config(self._config()))
        reaction = dumped["network"]["reactions"]["reaction"]["CO._PHOTON__C.O"]

        assert dumped["jaffgen"]["output_dir"] == "/abs/out"
        assert reaction["photo"]["shielding"] == {
            "type": "leiden",
            "shielded_by": ["self"],
        }


class TestDumpLayout:
    def test_jaffgen_section_comes_first(self):
        text = dump_config(JaffgenTomlConfig())

        assert text.index("[jaffgen]") < text.index("[network]")

    def test_empty_tables_are_left_out(self):
        assert "[network.reactions.reaction]" not in dump_config(JaffgenTomlConfig())

    def test_schema_line(self):
        assert dump_config(
            JaffgenTomlConfig(), schema=Path("/s/jaffgen.schema.json")
        ).startswith("#:schema /s/jaffgen.schema.json\n")
