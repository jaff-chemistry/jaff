# ABOUTME: Tests for jaffgen/jaffx CLI wiring: command-line options become a config
# ABOUTME: layer (CLI_KEYS) over jaffgen.toml and jaff.toml, end to end through jaffgen

import shutil
import tomllib
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from jaff import Network

from jaff.cli._helper import cli_layer, command_line_layer, parse_set
from jaff.cli.jaffgen.config import load_jaffgen_config
from jaff.errors import ParserError
from jaff.physics.thermodynamics.eos.config import IdealEosConfig
from tests._paths import FIXTURES_DIR
from tests.jaffgen_run import run_folder

FIXTURES = FIXTURES_DIR
DUP = str(FIXTURES / "duplicate_temp_range.dat")


class TestCliLayer:
    """Each option given sets its config key path; options not given set nothing."""

    def test_options_map_to_key_paths(self):
        args = SimpleNamespace(
            network="net.dat",
            template="microphysics",
            files=["a.cpp"],
            outdir="out",
            label="x",
            strict=True,
            duplicate_policy="error",
        )

        assert cli_layer(args) == {
            "jaffgen": {
                "network_file": "net.dat",
                "template": "microphysics",
                "input_files": ["a.cpp"],
                "output_dir": "out",
            },
            "network": {
                "label": "x",
                "strict": True,
                "reactions": {"duplicate_policy": "error"},
            },
        }

    def test_options_not_given_set_nothing(self):
        assert cli_layer(SimpleNamespace(label=None, funcfile=None)) == {}

    def test_false_is_given(self):
        assert cli_layer(SimpleNamespace(funcfile=False, echo=False)) == {
            "jaffgen": {"echo": False},
            "network": {"funcfile": False},
        }


class TestCliOverConfigFile:
    """The command-line layer wins over jaffgen.toml, which wins over jaff.toml."""

    def test_cli_flag_overrides_config_value(self, tmp_path):
        cfg = tmp_path / "jaffgen.toml"
        cfg.write_text('[network.reactions]\nduplicate_policy = "preserve-last"\n')
        cli = cli_layer(SimpleNamespace(duplicate_policy="error"))

        config = load_jaffgen_config(cfg, cli, cwd=tmp_path)

        assert config.network["reactions"]["duplicate_policy"] == "error"

    def test_no_cli_flag_keeps_config_value(self, tmp_path):
        cfg = tmp_path / "jaffgen.toml"
        cfg.write_text('[network.reactions]\nduplicate_policy = "preserve-last"\n')

        config = load_jaffgen_config(cfg, cli_layer(SimpleNamespace()), cwd=tmp_path)

        assert config.network["reactions"]["duplicate_policy"] == "preserve-last"

    def test_eos_table(self, tmp_path):
        cfg = tmp_path / "jaffgen.toml"
        cfg.write_text('[network.thermodynamics.eos]\ntype = "ideal"\ngamma = 1.4\n')

        layer = load_jaffgen_config(cfg, cwd=tmp_path).network
        eos = Network(DUP, config=layer, funcfile=False).config.thermodynamics.eos

        assert isinstance(eos, IdealEosConfig)
        assert eos.gamma == pytest.approx(1.4)

    @pytest.mark.parametrize("index, expected", [("1.5", 1.5), ("[2, 1.0]", [2, 1.0])])
    def test_profile_index_scalar_or_per_band(self, tmp_path, index, expected):
        cfg = tmp_path / "jaffgen.toml"
        cfg.write_text(
            f"[network.radiation]\nbands = [6, 11.2, 13.6]\nprofile_index = {index}\n"
        )

        radiation = load_jaffgen_config(cfg, cwd=tmp_path).network["radiation"]

        assert radiation["profile_index"] == expected

    def test_non_boolean_echo_is_an_error(self, tmp_path):
        cfg = tmp_path / "jaffgen.toml"
        cfg.write_text('[jaffgen]\necho = "no"\n')

        with pytest.raises(ParserError, match=r"jaffgen\.echo"):
            load_jaffgen_config(cfg, cwd=tmp_path)

    def test_old_keys_are_reported(self, tmp_path):
        cfg = tmp_path / "jaffgen.toml"
        cfg.write_text('[network]\nerrors = true\nconfig = "jaff.toml"\n')

        layer = load_jaffgen_config(cfg, cwd=tmp_path).network
        with pytest.raises(ParserError) as info:
            Network(DUP, config=layer, funcfile=False)

        assert "network.errors: unknown key" in str(info.value)
        assert "network.config: unknown key" in str(info.value)


class TestTyperFlags:
    """The typer command turns flags into the args namespace run_config reads."""

    def _args(self, *flags):
        from typer.testing import CliRunner

        import jaff.cli.jaffgen._engine as engine

        captured = {}

        def _capture(args):
            captured["args"] = args

        with (
            patch.object(engine, "run_config", _capture),
            patch.object(engine, "JaffGen"),
        ):
            CliRunner().invoke(engine.app, ["--network", "n", *flags])

        return captured["args"]

    def test_network_config_flag(self):
        assert self._args("--network-config", "x.toml").network_config == "x.toml"

    def test_network_config_defaults_to_none(self):
        assert self._args().network_config is None

    def test_strict_flag(self):
        assert self._args("--strict").strict is True
        assert self._args("--no-strict").strict is False


class TestJaffxWiring:
    """jaffx forwards its network options to Network."""

    def _args(self, duplicate_policy):
        return SimpleNamespace(
            network=DUP,
            funcfile=False,
            label=None,
            expand_nuclei=None,
            duplicate_policy=duplicate_policy,
        )

    def test_flag_applied(self):
        from jaff.cli.jaffx._engine import JaffX

        jx = JaffX.__new__(JaffX)
        net = jx.get_network(self._args("preserve-last"))
        assert net.config.reactions.duplicate_policy == "preserve-last"

    def test_none_uses_default(self):
        from jaff.cli.jaffx._engine import JaffX

        jx = JaffX.__new__(JaffX)
        net = jx.get_network(self._args(None))
        assert net.config.reactions.duplicate_policy == "preserve-first"


class TestFuncfileWiring:
    """--funcfile false disables aux loading, over a sibling file or a config path."""

    # chemRateCoeff0 overrides reaction 0's rate of 1 with 9 when aux loading is on.
    NETWORK = "H + H -> H2 [10,1000] 1\nH2 -> H + H [10,1000] 1\n"
    JFUNC = "@function chemRateCoeff0(tgas)\n    @return 9\n"
    TEMPLATE = "$$JAFF REPEAT idx, rate IN rates\nk[$idx$] = $rate$\n$$JAFF END\n"

    def _generate(self, tmp_path, *extra):
        from typer.testing import CliRunner

        from jaff.cli.jaffgen._engine import app

        (tmp_path / "net.dat").write_text(self.NETWORK)
        (tmp_path / "rates.py").write_text(self.TEMPLATE)
        # A fresh output directory, so it holds this run's folder alone.
        shutil.rmtree(tmp_path / "out", ignore_errors=True)
        args = ["--network", str(tmp_path / "net.dat")]
        args += ["--files", str(tmp_path / "rates.py")]
        args += ["--outdir", str(tmp_path / "out"), *extra]
        result = CliRunner().invoke(app, args)
        assert result.exit_code == 0, result.output
        return (run_folder(tmp_path / "out") / "rates.py").read_text()

    def test_cli_false_skips_sibling_jfunc(self, tmp_path):
        (tmp_path / "net.jfunc").write_text(self.JFUNC)
        assert "k[0] = 9" in self._generate(tmp_path)
        assert "k[0] = 1" in self._generate(tmp_path, "--funcfile", "false")

    def test_cli_false_skips_config_funcfile(self, tmp_path):
        aux = tmp_path / "aux.jfunc"
        aux.write_text(self.JFUNC)
        cfg = tmp_path / "jaffgen.toml"
        cfg.write_text(f'[network]\nfuncfile = "{aux.as_posix()}"\n')
        # The config is rendered alongside the templates and copied unchanged.
        with_cfg = ("--config", str(cfg))
        assert "k[0] = 9" in self._generate(tmp_path, *with_cfg)
        assert "k[0] = 1" in self._generate(tmp_path, *with_cfg, "--funcfile", "false")


class TestConfigRelativeNetworkPaths:
    """Paths in jaffgen.toml resolve against its directory; CLI paths against the CWD."""

    NETWORK = TestFuncfileWiring.NETWORK
    TEMPLATE = TestFuncfileWiring.TEMPLATE

    @staticmethod
    def _jfunc(rate):
        return f"@function chemRateCoeff0(tgas)\n    @return {rate}\n"

    @pytest.fixture
    def dirs(self, tmp_path, monkeypatch):
        """Config dir with real files; CWD elsewhere holding same-named decoys."""
        cfg_dir, cwd = tmp_path / "cfg", tmp_path / "cwd"
        cfg_dir.mkdir()
        cwd.mkdir()
        (cfg_dir / "net.dat").write_text(self.NETWORK)
        (cfg_dir / "rates.py").write_text(self.TEMPLATE)
        (cfg_dir / "net.jfunc").write_text(self._jfunc(9))
        (cfg_dir / "jaff.toml").write_text('[network.reactions]\nt_cutoff = "clip"\n')
        (cwd / "net.jfunc").write_text(self._jfunc(5))
        (cwd / "jaff.toml").write_text("this is [not valid toml\n")
        monkeypatch.chdir(cwd)
        return cfg_dir, cwd

    def _write_cfg(self, cfg_dir, jaffgen_block="", network_block=""):
        cfg = cfg_dir / "jaffgen.toml"
        cfg.write_text(
            '[jaffgen]\nnetwork_file = "net.dat"\ninput_files = ["rates.py"]\n'
            f"{jaffgen_block}[network]\n{network_block}"
        )
        return cfg

    def test_funcfile_resolves_against_config_dir(self, dirs):
        cfg_dir, _ = dirs
        cfg = self._write_cfg(cfg_dir, network_block='funcfile = "net.jfunc"\n')

        assert load_jaffgen_config(cfg).network["funcfile"] == cfg_dir / "net.jfunc"

    def test_network_config_resolves_against_config_dir(self, dirs):
        cfg_dir, _ = dirs
        cfg = self._write_cfg(cfg_dir, jaffgen_block='network_config = "jaff.toml"\n')

        assert load_jaffgen_config(cfg).jaffgen.network_config == cfg_dir / "jaff.toml"

    def test_boolean_funcfile_is_preserved(self, dirs):
        cfg_dir, _ = dirs
        cfg = self._write_cfg(cfg_dir, network_block="funcfile = false\n")

        assert load_jaffgen_config(cfg).network["funcfile"] is False

    def test_absolute_paths_are_kept(self, dirs):
        cfg_dir, cwd = dirs
        block = f'funcfile = "{(cwd / "net.jfunc").as_posix()}"\n'
        cfg = self._write_cfg(cfg_dir, network_block=block)

        assert load_jaffgen_config(cfg).network["funcfile"] == cwd / "net.jfunc"

    def _generate(self, cfg, *extra):
        from typer.testing import CliRunner

        from jaff.cli.jaffgen._engine import app

        out = cfg.parent / "out"
        args = ["--config", str(cfg), "--outdir", str(out), *extra]
        result = CliRunner().invoke(app, args)
        assert result.exit_code == 0, result.output
        return (run_folder(out) / "rates.py").read_text()

    def test_cli_uses_config_relative_funcfile(self, dirs):
        cfg_dir, _ = dirs
        cfg = self._write_cfg(cfg_dir, network_block='funcfile = "net.jfunc"\n')
        assert "k[0] = 9" in self._generate(cfg)

    def test_cli_uses_config_relative_network_config(self, dirs):
        # The CWD decoy jaff.toml is malformed; loading it would fail the run.
        # funcfile defaults to the sibling scan, which finds cfg/net.jfunc.
        cfg_dir, _ = dirs
        cfg = self._write_cfg(cfg_dir, jaffgen_block='network_config = "jaff.toml"\n')
        assert "k[0] = 9" in self._generate(cfg)

    def test_cli_funcfile_override_stays_cwd_relative(self, dirs):
        cfg_dir, _ = dirs
        cfg = self._write_cfg(cfg_dir, network_block='funcfile = "net.jfunc"\n')
        assert "k[0] = 5" in self._generate(cfg, "--funcfile", "net.jfunc")


class TestEchoWiring:
    """[jaffgen] echo and --echo/--no-echo decide whether directives are echoed."""

    NETWORK = TestFuncfileWiring.NETWORK
    TEMPLATE = TestFuncfileWiring.TEMPLATE

    def _generate(self, tmp_path, *extra, echo=None):
        from typer.testing import CliRunner

        from jaff.cli.jaffgen._engine import app

        (tmp_path / "net.dat").write_text(self.NETWORK)
        (tmp_path / "rates.py").write_text(self.TEMPLATE)
        args = ["--network", str(tmp_path / "net.dat")]
        args += ["--files", str(tmp_path / "rates.py")]
        args += ["--outdir", str(tmp_path / "out"), *extra]
        if echo is not None:
            cfg = tmp_path / "jaffgen.toml"
            cfg.write_text(f"[jaffgen]\necho = {echo}\n")
            args += ["--config", str(cfg)]
        result = CliRunner().invoke(app, args)
        assert result.exit_code == 0, result.output
        return (run_folder(tmp_path / "out") / "rates.py").read_text()

    def test_echo_is_on_by_default(self, tmp_path):
        assert "# $$JAFF REPEAT" in self._generate(tmp_path)

    def test_config_turns_echo_off(self, tmp_path):
        out = self._generate(tmp_path, echo="false")
        assert "$$JAFF" not in out and "k[0] = 1" in out

    def test_cli_overrides_config(self, tmp_path):
        (tmp_path / "a").mkdir()
        (tmp_path / "b").mkdir()

        assert "$$JAFF" not in self._generate(tmp_path / "a", "--no-echo", echo="true")
        assert "# $$JAFF" in self._generate(tmp_path / "b", "--echo", echo="false")


class TestDataTables:
    """[[jaffgen.table]] entries are converted into the output directory."""

    NETWORK = TestFuncfileWiring.NETWORK
    TEMPLATE = TestFuncfileWiring.TEMPLATE

    def test_csv_table_written(self, tmp_path):
        from typer.testing import CliRunner

        from jaff.cli.jaffgen._engine import app

        (tmp_path / "net.dat").write_text(self.NETWORK)
        (tmp_path / "rates.py").write_text(self.TEMPLATE)
        (tmp_path / "data.csv").write_text("# a comment\nT0 N\n1 2\n3 4\n")
        cfg = tmp_path / "jaffgen.toml"
        cfg.write_text(
            '[jaffgen]\nnetwork_file = "net.dat"\ninput_files = ["rates.py"]\n'
            "[[jaffgen.table]]\n"
            '[jaffgen.table.source]\npath = "data.csv"\ndelimiter = " "\n'
            'comment = "#"\ncols = ["T0"]\n'
            '[jaffgen.table.target]\npath = "out.csv"\ndelimiter = ","\n'
        )

        out = tmp_path / "out"
        result = CliRunner().invoke(app, ["--config", str(cfg), "--outdir", str(out)])

        assert result.exit_code == 0, result.output
        assert (run_folder(out) / "out.csv").exists()


class TestParseSet:
    """--set key.path=value: key path and value written as in TOML."""

    @pytest.mark.parametrize(
        "item, expected",
        [
            ("network.strict=true", {"network": {"strict": True}}),
            (
                "network.radiation.bands=[6, 13.6]",
                {"network": {"radiation": {"bands": [6, 13.6]}}},
            ),
            ('network.label="gow"', {"network": {"label": "gow"}}),
            ("network.label=gow", {"network": {"label": "gow"}}),
            ("network.label=a=b", {"network": {"label": "a=b"}}),
            (
                'network.reactions.reaction."H2._PHOTON__H.H".t_cutoff=clip',
                {
                    "network": {
                        "reactions": {
                            "reaction": {"H2._PHOTON__H.H": {"t_cutoff": "clip"}}
                        }
                    }
                },
            ),
        ],
    )
    def test_item(self, item, expected):
        assert parse_set([item]) == expected

    def test_items_merge_later_wins(self):
        items = ["network.label=a", "network.strict=true", "network.label=b"]

        assert parse_set(items) == {"network": {"label": "b", "strict": True}}

    @pytest.mark.parametrize("item", ["network.label", "=1", "network..x=1"])
    def test_invalid_item(self, item):
        with pytest.raises(ParserError, match="--set"):
            parse_set([item])


class TestCommandLineLayer:
    def test_set_beats_named_option(self):
        args = SimpleNamespace(label="flag", set=["network.label=set"])

        assert command_line_layer(args) == {"network": {"label": "set"}}

    def test_without_set(self):
        assert command_line_layer(SimpleNamespace(label="flag")) == {
            "network": {"label": "flag"}
        }


class TestSetEndToEnd:
    def test_jaffgen_set_disables_funcfile(self, tmp_path):
        out = TestFuncfileWiring()
        (tmp_path / "net.jfunc").write_text(out.JFUNC)

        assert "k[0] = 1" in out._generate(tmp_path, "--set", "network.funcfile=false")

    def test_jaffx_set_applies(self):
        from jaff.cli.jaffx._engine import JaffX

        args = SimpleNamespace(
            network=DUP,
            funcfile=False,
            label=None,
            expand_nuclei=None,
            duplicate_policy=None,
            set=['network.reactions.duplicate_policy="preserve-last"'],
        )

        net = JaffX.__new__(JaffX).get_network(args)

        assert net.config.reactions.duplicate_policy == "preserve-last"

    def test_jaffx_rejects_jaffgen_keys(self):
        from jaff.cli.jaffx._engine import JaffX

        args = SimpleNamespace(
            network=DUP,
            funcfile=False,
            label=None,
            expand_nuclei=None,
            duplicate_policy=None,
            set=["jaffgen.echo=false"],
        )

        with pytest.raises(ParserError, match="jaffgen"):
            JaffX.__new__(JaffX).get_network(args)

    def test_typer_collects_repeated_set(self):
        args = TestTyperFlags()._args(
            "--set", "network.label=a", "--set", "network.strict=true"
        )

        assert args.set == ["network.label=a", "network.strict=true"]


class TestConfigToml:
    """jaffgen writes the run's whole config as config.toml, not a jaffgen.toml copy."""

    def test_written_instead_of_jaffgen_copy(self, tmp_path):
        cfg = tmp_path / "jaffgen.toml"
        cfg.write_text('[network.reactions]\nt_cutoff = "extrapolate"\n')
        gen = TestFuncfileWiring()
        gen._generate(tmp_path, "--config", str(cfg), "--label", "from_cli")

        out = run_folder(tmp_path / "out")
        written = tomllib.loads((out / "config.toml").read_text())
        assert not (out / "jaffgen.toml").exists()
        assert written["network"]["label"] == "from_cli"
        assert written["network"]["reactions"]["t_cutoff"] == "extrapolate"
        assert written["jaffgen"]["network_file"] == str((tmp_path / "net.dat").resolve())

    def test_reloads_as_the_same_config(self, tmp_path):
        TestFuncfileWiring()._generate(tmp_path, "--label", "x")
        written = run_folder(tmp_path / "out") / "config.toml"

        reloaded = load_jaffgen_config(written)

        assert reloaded.network["label"] == "x"
        assert reloaded.jaffgen.output_dir == (tmp_path / "out").resolve()


class TestCliErrors:
    """Config and input errors print their message, not a traceback."""

    def test_jaffgen_config_error_is_one_message(self, tmp_path):
        from typer.testing import CliRunner

        from jaff.cli.jaffgen._engine import app

        cfg = tmp_path / "jaffgen.toml"
        cfg.write_text('[network.radiation]\nbands = [6.0, 13.6]\nrsl = "c_hat"\n')

        result = CliRunner().invoke(app, ["--config", str(cfg), "--network", DUP])

        assert result.exit_code == 1
        assert "network.radiation.rsl: unknown key" in result.output
        assert "Traceback" not in result.output
        assert not isinstance(result.exception, ParserError)

    def test_jaffx_config_error_is_one_message(self):
        from typer.testing import CliRunner

        from jaff.cli.jaffx._engine import app

        result = CliRunner().invoke(
            app,
            [
                "get",
                "num-species",
                "--network",
                DUP,
                "--set",
                "network.reactions.t_cutoff=clamp",
            ],
        )

        assert result.exit_code == 1
        assert "network.reactions.t_cutoff" in result.output
        assert not isinstance(result.exception, ParserError)


class TestJaffGenEngine:
    """JaffGen(config) prepares a run from a loaded config; run() writes it."""

    def _config(self, tmp_path, **options):
        from jaff.cli.jaffgen._engine import run_config

        (tmp_path / "net.dat").write_text(TestFuncfileWiring.NETWORK)
        (tmp_path / "rates.py").write_text(TestFuncfileWiring.TEMPLATE)
        args = SimpleNamespace(
            network=str(tmp_path / "net.dat"),
            files=[str(tmp_path / "rates.py")],
            outdir=str(tmp_path / "out"),
            **options,
        )

        return run_config(args)

    def test_run_config_is_the_loaded_config(self, tmp_path):
        from jaff.cli.jaffgen.config import JaffgenRunConfig

        config = self._config(tmp_path, strict=True)

        assert isinstance(config, JaffgenRunConfig)
        assert config.network == {"strict": True}
        assert config.jaffgen.output_dir == (tmp_path / "out").resolve()

    def test_constructor_prepares_without_writing(self, tmp_path):
        from jaff.cli import JaffGen

        config = self._config(tmp_path)
        gen = JaffGen(config)

        assert gen.config is config
        assert gen.net.reactions.count == 2
        assert [f.relpath.name for f in gen.output_files] == ["rates.py"]
        assert gen.output_root == (tmp_path / "out").resolve()
        assert gen.output_dir.parent == gen.output_root
        assert not (tmp_path / "out").exists()

    def test_run_writes_the_files_and_the_config(self, tmp_path):
        from jaff.cli import JaffGen

        gen = JaffGen(self._config(tmp_path))
        gen.run()

        assert (gen.output_dir / "rates.py").exists()
        assert (gen.output_dir / "config.toml").exists()

    def test_holds_no_cli_arguments_or_state(self, tmp_path):
        from jaff.cli import JaffGen

        gen = JaffGen(self._config(tmp_path))

        assert not hasattr(gen, "args")
        assert not hasattr(gen, "state")


class TestDuplicatePolicyChoices:
    def test_cli_choices_are_the_config_values(self):
        from typing import get_args

        from jaff.cli._helper import DuplicatePolicyChoice
        from jaff.core.reaction.config import DuplicatePolicy

        allowed = get_args(get_args(DuplicatePolicy)[0])

        assert [choice.value for choice in DuplicatePolicyChoice] == list(allowed)

    def test_cli_type_has_its_own_name(self):
        import jaff.cli._helper as helper

        assert not hasattr(helper, "DuplicatePolicy")


class TestSetTables:
    """[[jaffgen.table]] paths resolve against the layer that set them."""

    TABLE = 'jaffgen.table=[{source={path="data.csv"}}]'

    def test_set_table_resolves_against_the_working_directory(self, tmp_path):
        config = load_jaffgen_config(None, parse_set([self.TABLE]), cwd=tmp_path)

        assert config.jaffgen.table[0]["source"]["path"] == str(tmp_path / "data.csv")

    def test_file_table_resolves_against_the_file(self, tmp_path):
        (tmp_path / "run").mkdir()
        cfg = tmp_path / "run" / "jaffgen.toml"
        cfg.write_text('[[jaffgen.table]]\nsource = { path = "data.csv" }\n')

        config = load_jaffgen_config(cfg, cwd=tmp_path)

        assert config.jaffgen.table[0]["source"]["path"] == str(
            (tmp_path / "run" / "data.csv").resolve()
        )

    def test_default_source_is_kept(self, tmp_path):
        layer = parse_set(['jaffgen.table=[{source={path="default"}}]'])

        config = load_jaffgen_config(None, layer, cwd=tmp_path)

        assert config.jaffgen.table[0]["source"]["path"] == "default"

    def test_table_needs_no_config_file(self, tmp_path):
        from jaff.cli.jaffgen._config_table import ConfigTable

        source = tmp_path / "data.csv"
        source.write_text("a,b\n1,2\n")
        table = ConfigTable({"source": {"path": str(source)}}, None, tmp_path / "n.jet")

        assert table.source_config["path"] == str(source)


class TestOutputDirRecorded:
    def test_default_output_dir_is_in_config_toml(self, tmp_path, monkeypatch):
        import jaff.cli.jaffgen._engine as engine
        from jaff.cli import JaffGen

        monkeypatch.setattr(engine, "DEFAULT_OUTPUT", tmp_path / "generated")
        (tmp_path / "net.dat").write_text(TestFuncfileWiring.NETWORK)
        (tmp_path / "rates.py").write_text(TestFuncfileWiring.TEMPLATE)
        args = SimpleNamespace(
            network=str(tmp_path / "net.dat"), files=[str(tmp_path / "rates.py")]
        )

        gen = JaffGen(engine.run_config(args))
        gen.run()
        written = tomllib.loads((gen.output_dir / "config.toml").read_text())

        assert written["jaffgen"]["output_dir"] == str(tmp_path / "generated")

    def test_configured_output_dir_is_in_config_toml(self, tmp_path):
        """The base directory is recorded, not the run folder inside it."""
        from jaff.cli import JaffGen

        gen = JaffGen(TestRunFolder._config(tmp_path))
        gen.run()
        written = tomllib.loads((gen.output_dir / "config.toml").read_text())

        assert written["jaffgen"]["output_dir"] == str((tmp_path / "out").resolve())


class TestRunFolder:
    """Each run writes into a new <network>_<template>_<datetime> folder in output_dir."""

    STAMP = "20261031-101500"

    @staticmethod
    def _config(tmp_path, **options):
        from jaff.cli.jaffgen._engine import run_config

        (tmp_path / "net.dat").write_text(TestFuncfileWiring.NETWORK)
        (tmp_path / "rates.py").write_text(TestFuncfileWiring.TEMPLATE)
        options.setdefault("files", [str(tmp_path / "rates.py")])
        args = SimpleNamespace(
            network=str(tmp_path / "net.dat"), outdir=str(tmp_path / "out"), **options
        )

        return run_config(args)

    @pytest.fixture(autouse=True)
    def _fixed_time(self, monkeypatch):
        import jaff.cli.jaffgen._engine as engine

        monkeypatch.setattr(engine, "_timestamp", lambda: self.STAMP)

    def _gen(self, tmp_path, **options):
        from jaff.cli import JaffGen

        return JaffGen(self._config(tmp_path, **options))

    def test_named_after_network_template_and_time(self, tmp_path):
        gen = self._gen(tmp_path, files=None, template="python_solve_ivp")

        assert gen.output_dir == gen.output_root / f"net_python_solve_ivp_{self.STAMP}"

    def test_input_dir_name_without_template(self, tmp_path):
        (tmp_path / "mytemplates").mkdir()
        (tmp_path / "mytemplates" / "rates.py").write_text(TestFuncfileWiring.TEMPLATE)

        gen = self._gen(tmp_path, files=None, indir=str(tmp_path / "mytemplates"))

        assert gen.output_dir.name == f"net_mytemplates_{self.STAMP}"

    def test_custom_with_only_input_files(self, tmp_path):
        assert self._gen(tmp_path).output_dir.name == f"net_custom_{self.STAMP}"

    def test_label_names_the_network(self, tmp_path):
        gen = self._gen(tmp_path, label="my net/v2")

        assert gen.output_dir.name == f"my_net_v2_custom_{self.STAMP}"

    def test_same_second_runs_get_a_suffix(self, tmp_path):
        first = self._gen(tmp_path)
        first.run()
        second = self._gen(tmp_path)
        second.run()
        third = self._gen(tmp_path)

        assert second.output_dir.name == f"net_custom_{self.STAMP}-2"
        assert third.output_dir.name == f"net_custom_{self.STAMP}-3"
        assert (first.output_dir / "rates.py").exists()

    def test_timestamp_format(self):
        import re

        from jaff.cli.jaffgen._engine import _timestamp

        assert re.fullmatch(r"\d{8}-\d{6}", _timestamp())
