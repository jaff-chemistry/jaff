# ABOUTME: JaffgenConfig ([jaffgen]) and the per-file models of jaff.toml
# ABOUTME: and jaffgen.toml; two-phase load (jaff.toml < jaffgen.toml < CLI)

from pathlib import Path

import pytest
from pydantic import ValidationError

import jaff.cli.jaffgen.config as config_module
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


NET = "H + H -> H2  []  1e-10\nH2 -> H + H  []  1e-10\n"


@pytest.fixture
def project(tmp_path):
    """A network dir and a run dir; write(name, text) writes a file under tmp_path."""
    (tmp_path / "net").mkdir()
    (tmp_path / "net" / "net.jet").write_text(NET)
    (tmp_path / "run").mkdir()

    def write(name, text):
        path = tmp_path / name
        path.write_text(text)

        return path

    return tmp_path, write


def _network(cfg):
    """The Network a jaffgen run builds from *cfg* (its [network] layer as config=)."""
    from jaff import Network

    return Network(cfg.jaffgen.network_file, config=cfg.network, funcfile=False)


class TestLoadJaffgenConfig:
    def test_nothing_gives_defaults(self, tmp_path):
        cfg = load_jaffgen_config(cwd=tmp_path)

        assert cfg.jaffgen == JaffgenConfig()
        assert cfg.network == {}
        assert cfg.config_file is None

    def test_network_layer_is_jaffgen_toml_then_cli(self, project):
        root, write = project
        cfg_file = write(
            "run/jaffgen.toml",
            '[jaffgen]\nnetwork_file = "../net/net.jet"\n'
            '[network]\nlabel = "jaffgen"\nstrict = false\n',
        )

        cfg = load_jaffgen_config(cfg_file, {"network": {"label": "cli"}}, cwd=root)

        assert cfg.network == {"label": "cli", "strict": False}
        assert cfg.config_file == cfg_file

    def test_network_applies_the_jaff_toml_beside_it(self, project):
        root, write = project
        write(
            "net/jaff.toml",
            '[network]\nlabel = "jaff"\nstrict = true\nexpand_nuclei = false\n',
        )
        cfg_file = write(
            "run/jaffgen.toml",
            '[jaffgen]\nnetwork_file = "../net/net.jet"\n'
            '[network]\nlabel = "jaffgen"\nstrict = false\n',
        )

        cfg = load_jaffgen_config(cfg_file, {"network": {"label": "cli"}}, cwd=root)
        net = _network(cfg)

        assert "expand_nuclei" not in cfg.network  # read by Network, not jaffgen
        assert net.config.label == "cli"
        assert net.config.strict is False
        assert net.config.expand_nuclei is False

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

    def test_custom_jaff_toml_is_merged_over_the_one_beside_the_network(self, project):
        root, write = project
        write("net/jaff.toml", '[network]\nlabel = "beside"\nstrict = true\n')
        write("run/custom.toml", '[network]\nlabel = "custom"\n')
        cfg_file = write(
            "run/jaffgen.toml",
            '[jaffgen]\nnetwork_file = "../net/net.jet"\nnetwork_config = "custom.toml"\n',
        )

        cfg = load_jaffgen_config(cfg_file, cwd=root)
        net = _network(cfg)

        assert cfg.network == {"label": "custom"}
        assert net.config.label == "custom"
        assert net.config.strict is True

    def test_custom_jaff_toml_beside_the_network_is_not_read_twice(self, project):
        root, write = project
        write("net/jaff.toml", '[network]\nlabel = "beside"\n')
        cfg_file = write(
            "run/jaffgen.toml",
            '[jaffgen]\nnetwork_file = "../net/net.jet"\n'
            'network_config = "../net/jaff.toml"\n',
        )

        assert load_jaffgen_config(cfg_file, cwd=root).network == {}

    def test_missing_network_config(self, project):
        root, write = project
        cfg_file = write("run/jaffgen.toml", '[jaffgen]\nnetwork_config = "nope.toml"\n')

        with pytest.raises(FileNotFoundError):
            load_jaffgen_config(cfg_file, cwd=root)

    def test_jaffgen_section_in_custom_jaff_toml_is_rejected(self, project):
        root, write = project
        write("run/custom.toml", '[jaffgen]\ntemplate = "microphysics"\n')

        with pytest.raises(ParserError, match="jaffgen"):
            load_jaffgen_config(
                None, {"jaffgen": {"network_config": "run/custom.toml"}}, cwd=root
            )

    def test_jaffgen_section_is_validated(self, project):
        root, _ = project

        with pytest.raises(ParserError, match=r"jaffgen\.echo"):
            load_jaffgen_config(None, {"jaffgen": {"echo": "yes"}}, cwd=root)

    def test_unknown_top_level_table_is_rejected(self, project):
        root, _ = project

        with pytest.raises(ParserError, match="netwrok"):
            load_jaffgen_config(None, {"netwrok": {"strict": True}}, cwd=root)

    def test_network_keys_are_checked_by_the_network(self, project):
        root, _ = project
        cfg = load_jaffgen_config(
            None,
            {"jaffgen": {"network_file": "net/net.jet"}, "network": {"rates": {}}},
            cwd=root,
        )

        with pytest.raises(ParserError, match=r"network\.rates: unknown key"):
            _network(cfg)


BUNDLE = '[network.radiation]\nbands = [13.6, "inf"]\n'


@pytest.fixture
def bundled_template(tmp_path, monkeypatch):
    """A built-in template "bundled" whose jaffgen.toml sets radiation bands.

    Test-owned, so the tests do not depend on what the real templates bundle.
    """
    templates = tmp_path / "templates"
    (templates / "bundled").mkdir(parents=True)
    (templates / "bundled" / "jaffgen.toml").write_text(BUNDLE)
    monkeypatch.setattr(config_module, "TEMPLATES_DIR", templates)
    monkeypatch.setattr(config_module, "predefined_templates", lambda: {"bundled"})

    return "bundled"


class TestBundledTemplateConfig:
    def test_template_bundle_is_the_file_layer(self, project, bundled_template):
        root, _ = project
        cfg = load_jaffgen_config(
            None, {"jaffgen": {"template": bundled_template}}, cwd=root
        )

        assert cfg.config_file is not None
        assert cfg.config_file.name == "jaffgen.toml"
        assert cfg.network["radiation"]["bands"] == [13.6, "inf"]

    def test_cli_key_needing_bundle_keys_is_valid(self, project, bundled_template):
        # bugs.md #46: --set network.radiation.c=1e10 with a template whose
        # bundle sets bands used to fail ("bands: required key missing").
        root, _ = project
        cli = {
            "jaffgen": {"template": bundled_template, "network_file": "net/net.jet"},
            "network": {"radiation": {"c": 1e10}},
        }

        cfg = load_jaffgen_config(None, cli, cwd=root)

        assert cfg.network["radiation"]["bands"] == [13.6, "inf"]
        assert cfg.network["radiation"]["c"] == pytest.approx(1e10)

    def test_explicit_config_wins_over_bundle(self, project, bundled_template):
        root, write = project
        cfg_file = write("run/jaffgen.toml", "[jaffgen]\n")

        cfg = load_jaffgen_config(
            cfg_file, {"jaffgen": {"template": bundled_template}}, cwd=root
        )

        assert cfg.config_file == cfg_file
        assert "radiation" not in cfg.network
