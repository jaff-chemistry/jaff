# ABOUTME: JaffgenConfig ([jaffgen]) and the per-file models of jaff.toml
# ABOUTME: and jaffgen.toml

from pathlib import Path

import pytest
from pydantic import ValidationError

from jaff.cli.jaffgen.config import JaffgenConfig, JaffgenTomlConfig
from jaff.core.network.config import JaffTomlConfig


class TestJaffgenConfig:
    def test_defaults(self):
        cfg = JaffgenConfig()

        assert cfg.template is None
        assert cfg.input_dir is None
        assert cfg.input_files == []
        assert cfg.output_dir is None
        assert cfg.network_file is None
        assert cfg.network_config is None
        assert cfg.echo is True
        assert cfg.table == []

    def test_paths_become_paths(self):
        cfg = JaffgenConfig.model_validate(
            {
                "input_dir": "templates",
                "input_files": ["a.cpp", "b.f90"],
                "output_dir": "out",
                "network_file": "GOW",
                "network_config": "net/jaff.toml",
            }
        )

        assert cfg.input_dir == Path("templates")
        assert cfg.input_files == [Path("a.cpp"), Path("b.f90")]
        assert cfg.output_dir == Path("out")
        assert cfg.network_file == Path("GOW")
        assert cfg.network_config == Path("net/jaff.toml")

    def test_tables_are_passed_through_unchecked(self):
        table = {"source": {"path": "default"}, "target": {"anything": [1, "x"]}}

        assert JaffgenConfig.model_validate({"table": [table]}).table == [table]

    def test_rejects_string_echo(self):
        with pytest.raises(ValidationError, match="echo"):
            JaffgenConfig.model_validate({"echo": "false"})

    def test_rejects_unknown_key(self):
        with pytest.raises(ValidationError, match="indir"):
            JaffgenConfig.model_validate({"indir": "templates"})


class TestTomlConfigs:
    def test_jaff_toml_holds_network_only(self):
        cfg = JaffTomlConfig.model_validate({"network": {"label": "gow"}})

        assert cfg.network.label == "gow"

    def test_jaff_toml_rejects_jaffgen_section(self):
        with pytest.raises(ValidationError, match="jaffgen"):
            JaffTomlConfig.model_validate({"jaffgen": {"template": "microphysics"}})

    def test_jaffgen_toml_holds_both_sections(self):
        cfg = JaffgenTomlConfig.model_validate(
            {"jaffgen": {"template": "microphysics"}, "network": {"strict": True}}
        )

        assert cfg.jaffgen.template == "microphysics"
        assert cfg.network.strict is True

    def test_empty_files_are_valid(self):
        assert JaffTomlConfig.model_validate({}) == JaffTomlConfig()
        assert JaffgenTomlConfig.model_validate({}) == JaffgenTomlConfig()

    def test_rejects_top_level_table_array(self):
        with pytest.raises(ValidationError, match="table"):
            JaffgenTomlConfig.model_validate({"table": [{}]})
