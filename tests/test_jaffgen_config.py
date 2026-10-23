# ABOUTME: JaffgenConfig ([jaffgen]) and the per-file models of jaff.toml
# ABOUTME: and jaffgen.toml; two-phase load (jaff.toml < jaffgen.toml < CLI)

from pathlib import Path

import pytest
from pydantic import ValidationError

from jaff.cli.jaffgen.config import JaffgenConfig, JaffgenTomlConfig, load_jaffgen_config
from jaff.core.network.config import JaffTomlConfig
from jaff.errors import ParserError


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


@pytest.fixture
def project(tmp_path):
    """A network dir and a run dir; write(name, text) writes a file under tmp_path."""
    (tmp_path / "net").mkdir()
    (tmp_path / "net" / "net.jet").write_text("")
    (tmp_path / "run").mkdir()

    def write(name, text):
        path = tmp_path / name
        path.write_text(text)

        return path

    return tmp_path, write


class TestLoadJaffgenConfig:
    def test_nothing_gives_defaults(self, tmp_path):
        assert load_jaffgen_config(cwd=tmp_path) == JaffgenTomlConfig()

    def test_precedence_jaff_then_jaffgen_then_cli(self, project):
        root, write = project
        write(
            "net/jaff.toml",
            '[network]\nlabel = "jaff"\nstrict = true\nexpand_nuclei = false\n',
        )
        cfg_file = write(
            "run/jaffgen.toml",
            '[jaffgen]\nnetwork_file = "../net/net.jet"\n[network]\nlabel = "jaffgen"\nstrict = false\n',
        )

        cfg = load_jaffgen_config(cfg_file, {"network": {"label": "cli"}}, cwd=root)

        assert cfg.network.label == "cli"
        assert cfg.network.strict is False
        assert cfg.network.expand_nuclei is False

    def test_paths_resolved_per_source(self, project):
        root, write = project
        cfg_file = write(
            "run/jaffgen.toml", '[jaffgen]\noutput_dir = "out"\ninput_dir = "tpl"\n'
        )

        cfg = load_jaffgen_config(
            cfg_file, {"jaffgen": {"output_dir": "cli_out"}}, cwd=root
        )

        assert cfg.jaffgen.input_dir == (root / "run" / "tpl").resolve()
        assert cfg.jaffgen.output_dir == (root / "cli_out").resolve()

    def test_cli_network_file_selects_its_jaff_toml(self, project):
        root, write = project
        (root / "other").mkdir()
        (root / "other" / "o.jet").write_text("")
        write("net/jaff.toml", '[network]\nlabel = "net"\n')
        write("other/jaff.toml", '[network]\nlabel = "other"\n')
        cfg_file = write(
            "run/jaffgen.toml", '[jaffgen]\nnetwork_file = "../net/net.jet"\n'
        )

        cfg = load_jaffgen_config(
            cfg_file, {"jaffgen": {"network_file": "other/o.jet"}}, cwd=root
        )

        assert cfg.network.label == "other"

    def test_network_config_overrides_detection(self, project):
        root, write = project
        write("net/jaff.toml", '[network]\nlabel = "detected"\n')
        write("run/custom.toml", '[network]\nlabel = "custom"\n')
        cfg_file = write(
            "run/jaffgen.toml",
            '[jaffgen]\nnetwork_file = "../net/net.jet"\nnetwork_config = "custom.toml"\n',
        )

        assert load_jaffgen_config(cfg_file, cwd=root).network.label == "custom"

    def test_missing_network_config(self, project):
        root, write = project
        cfg_file = write("run/jaffgen.toml", '[jaffgen]\nnetwork_config = "nope.toml"\n')

        with pytest.raises(FileNotFoundError):
            load_jaffgen_config(cfg_file, cwd=root)

    def test_jaffgen_section_in_jaff_toml_is_rejected(self, project):
        root, write = project
        write("net/jaff.toml", '[jaffgen]\ntemplate = "microphysics"\n')

        with pytest.raises(ParserError, match="jaffgen"):
            load_jaffgen_config(
                None, {"jaffgen": {"network_file": "net/net.jet"}}, cwd=root
            )
