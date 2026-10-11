# ABOUTME: Tests for jaffgen/jaffx CLI wiring of network options
# ABOUTME: --network-config resolution, duplicate_policy/funcfile plumbing, [network] block

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

FIXTURES = Path(__file__).parent / "fixtures"
DUP = str(FIXTURES / "duplicate_temp_range.dat")
RXN = "H+.H2__H.H2+"


def _bare_jaffgen(network_config=None, duplicate_policy=None, funcfile=None):
    """A JaffGen instance with a minimal args namespace, no real CLI parse."""
    from jaff.cli.jaffgen._engine import JaffGen
    from jaff.cli.jaffgen._structs import State

    jg = JaffGen.__new__(JaffGen)
    jg.state = State()
    jg.args = SimpleNamespace(
        label=None,
        funcfile=funcfile,
        expand_nuclei=None,
        errors=None,
        network_config=network_config,
        duplicate_policy=duplicate_policy,
        lang=None,
    )
    return jg


class TestNetworkConfigArg:
    """--network-config resolution onto NetworkArgs.config."""

    def test_typer_callback_maps_flag(self):
        from typer.testing import CliRunner

        import jaff.cli.jaffgen._engine as engine

        captured = {}

        def _capture(self, args):
            captured["args"] = args

        with patch.object(engine.JaffGen, "__init__", _capture):
            runner = CliRunner()
            runner.invoke(engine.app, ["--network", "n", "--network-config", "x.toml"])
            assert captured["args"].network_config == "x.toml"
            runner.invoke(engine.app, ["--network", "n"])
            assert captured["args"].network_config is None

    def test_none_leaves_config_unset(self):
        jg = _bare_jaffgen(network_config=None)
        jg.set_network_options()
        assert jg.state.network_args.config is None

    def test_relative_path_resolves_against_cwd(self, tmp_path, monkeypatch):
        (tmp_path / "networks").mkdir()
        rel = Path("networks") / "jaff.toml"
        (tmp_path / rel).write_text('[network.rates]\nT_cutoff = "clip"\n')
        monkeypatch.chdir(tmp_path)

        jg = _bare_jaffgen(network_config=str(rel))
        jg.set_network_options()
        cfg = jg.state.network_args.config
        assert cfg == (tmp_path / rel).resolve() and cfg.is_absolute()

    def test_missing_path_raises(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        jg = _bare_jaffgen(network_config="does_not_exist.toml")
        with pytest.raises(FileNotFoundError):
            jg.set_network_options()

    def test_directory_path_raises(self, tmp_path, monkeypatch):
        (tmp_path / "networks").mkdir()
        monkeypatch.chdir(tmp_path)
        jg = _bare_jaffgen(network_config="networks")
        with pytest.raises(FileNotFoundError):
            jg.set_network_options()


class TestNetworkBlockExtraction:
    """[network.*] parses into the shapes fed to metadata slots."""

    def test_block_shapes(self, tmp_path):
        from jaff.drivers import Toml

        body = (
            '[network]\nlabel = "demo"\n'
            '[network.rates]\nT_cutoff = "extrapolate"\n'
            f'[network.reactions."{RXN}"]\nT_cutoff = "clip"\n'
            f'[network.reactions."{RXN}".shielding]\ntype = "leiden"\n'
        )
        p = tmp_path / "jaffgen.toml"
        p.write_text(body)
        network_cfg = Toml(str(p)).get_key("network") or {}

        assert (network_cfg.get("rates") or {}).get("T_cutoff") == "extrapolate"
        reactions_cfg = network_cfg.get("reactions") or {}
        assert reactions_cfg[RXN]["T_cutoff"] == "clip"
        assert reactions_cfg[RXN]["shielding"]["type"] == "leiden"
        assert network_cfg.get("label") == "demo"  # scalar still reachable


class TestDuplicatePolicyWiring:
    """jaffgen resolves duplicate_policy from CLI flag and jaffgen.toml."""

    def _bare(self, cli_value=None, toml_value=None):
        jg = _bare_jaffgen(duplicate_policy=cli_value)
        if toml_value is not None:
            jg.state.network_args.duplicate_policy = toml_value
        return jg

    def test_cli_flag_sets_network_args(self):
        jg = self._bare(cli_value="preserve-last")
        jg.set_network_options()
        assert jg.state.network_args.duplicate_policy == "preserve-last"

    def test_cli_flag_overrides_config_value(self):
        jg = self._bare(cli_value="error", toml_value="preserve-last")
        jg.set_network_options()
        assert jg.state.network_args.duplicate_policy == "error"

    def test_no_cli_flag_keeps_config_value(self):
        jg = self._bare(cli_value=None, toml_value="preserve-last")
        jg.set_network_options()
        assert jg.state.network_args.duplicate_policy == "preserve-last"

    def test_reads_key_from_config_file(self, tmp_path):
        from jaff.cli.jaffgen._engine import JaffGen
        from jaff.cli.jaffgen._structs import ResolvedPath, State
        from jaff.drivers import Toml

        cfg = tmp_path / "jaffgen.toml"
        cfg.write_text('[network]\nduplicate_policy = "error"\n')

        jg = JaffGen.__new__(JaffGen)
        jg.state = State()
        jg.state.config_dir = ResolvedPath(tmp_path, tmp_path)
        jg.state.config_raw = Toml(cfg)
        jg.set_state_from_config()
        assert jg.state.network_args.duplicate_policy == "error"


class TestJaffxWiring:
    """jaffx forwards the duplicate_policy flag onto NetworkArgs."""

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
        assert net.spec.duplicate_policy == "preserve-last"

    def test_none_uses_default(self):
        from jaff.cli.jaffx._engine import JaffX

        jx = JaffX.__new__(JaffX)
        net = jx.get_network(self._args(None))
        assert net.spec.duplicate_policy == "preserve-first"


class TestRadiationProfileIndexWiring:
    """[network.radiation] profile_index accepts a scalar or a per-band list."""

    def _from_config(self, tmp_path, profile_index):
        from jaff.cli.jaffgen._engine import JaffGen
        from jaff.cli.jaffgen._structs import ResolvedPath, State
        from jaff.drivers import Toml

        cfg = tmp_path / "jaffgen.toml"
        cfg.write_text(
            f"[network.radiation]\nbands = [6, 11.2, 13.6]\n"
            f"profile_index = {profile_index}\n"
        )
        jg = JaffGen.__new__(JaffGen)
        jg.state = State()
        jg.state.config_dir = ResolvedPath(tmp_path, tmp_path)
        jg.state.config_raw = Toml(cfg)
        jg.set_state_from_config()
        return jg.state.network_args

    def test_scalar_profile_index(self, tmp_path):
        assert self._from_config(tmp_path, "1.5").rad_profile_index == 1.5

    def test_list_profile_index(self, tmp_path):
        assert self._from_config(tmp_path, "[2, 1.0]").rad_profile_index == [2, 1.0]


class TestEosWiring:
    """[network.eos] maps onto NetworkArgs.eos and from there to EosProps."""

    def _from_config(self, tmp_path, block):
        from jaff.cli.jaffgen._engine import JaffGen
        from jaff.cli.jaffgen._structs import ResolvedPath, State
        from jaff.drivers import Toml

        cfg = tmp_path / "jaffgen.toml"
        cfg.write_text(block)
        jg = JaffGen.__new__(JaffGen)
        jg.state = State()
        jg.state.config_dir = ResolvedPath(tmp_path, tmp_path)
        jg.state.config_raw = Toml(cfg)
        jg.set_state_from_config()
        return jg.state.network_args

    def test_eos_table_is_stored(self, tmp_path):
        args = self._from_config(tmp_path, '[network.eos]\ntype = "ideal"\ngamma = 1.4\n')
        assert args.eos == {"type": "ideal", "gamma": 1.4}

    def test_missing_eos_table_leaves_none(self, tmp_path):
        assert self._from_config(tmp_path, '[network]\nlabel = "x"\n').eos is None


class TestFuncfileWiring:
    """--funcfile false must survive set_network_options and disable aux loading."""

    # chemRate0 overrides reaction 0's rate of 1 with 9 when aux loading is on.
    NETWORK = "H + H -> H2 [10,1000] 1\nH2 -> H + H [10,1000] 1\n"
    JFUNC = "@function chemRate0(tgas)\n    @return 9\n"
    TEMPLATE = "# $JAFF REPEAT idx, rate IN rates\nk[$idx$] = $rate$\n# $JAFF END\n"

    def _bare(self, cli_value=None, toml_value=None):
        jg = _bare_jaffgen(funcfile=cli_value)
        if toml_value is not None:
            jg.state.network_args.funcfile = toml_value
        return jg

    def test_cli_false_overrides_default(self):
        jg = self._bare(cli_value=False)
        jg.set_network_options()
        assert jg.state.network_args.funcfile is False

    def test_cli_false_overrides_config_path(self):
        jg = self._bare(cli_value=False, toml_value="aux.jfunc")
        jg.set_network_options()
        assert jg.state.network_args.funcfile is False

    def test_no_cli_flag_keeps_config_path(self):
        jg = self._bare(cli_value=None, toml_value="aux.jfunc")
        jg.set_network_options()
        assert jg.state.network_args.funcfile == "aux.jfunc"

    def _generate(self, tmp_path, *extra):
        from typer.testing import CliRunner

        from jaff.cli.jaffgen._engine import app

        (tmp_path / "net.dat").write_text(self.NETWORK)
        (tmp_path / "rates.py").write_text(self.TEMPLATE)
        args = ["--network", str(tmp_path / "net.dat")]
        args += ["--files", str(tmp_path / "rates.py")]
        args += ["--outdir", str(tmp_path / "out"), *extra]
        result = CliRunner().invoke(app, args)
        assert result.exit_code == 0, result.output
        return (tmp_path / "out" / "rates.py").read_text()

    def test_cli_false_skips_sibling_jfunc(self, tmp_path):
        (tmp_path / "net.jfunc").write_text(self.JFUNC)
        assert "k[0] = 9" in self._generate(tmp_path)
        assert "k[0] = 1" in self._generate(tmp_path, "--funcfile", "false")

    def test_cli_false_skips_config_funcfile(self, tmp_path):
        aux = tmp_path / "aux.jfunc"
        aux.write_text(self.JFUNC)
        cfg = tmp_path / "jaffgen.toml"
        cfg.write_text(f'[network]\nfuncfile = "{aux.as_posix()}"\n')
        # The config is rendered alongside the templates; .toml needs --lang.
        with_cfg = ("--config", str(cfg), "--lang", "python")
        assert "k[0] = 9" in self._generate(tmp_path, *with_cfg)
        assert "k[0] = 1" in self._generate(tmp_path, *with_cfg, "--funcfile", "false")


class TestConfigRelativeNetworkPaths:
    """[network] funcfile/config in jaffgen.toml resolve against the config dir."""

    NETWORK = TestFuncfileWiring.NETWORK
    TEMPLATE = TestFuncfileWiring.TEMPLATE

    @staticmethod
    def _jfunc(rate):
        return f"@function chemRate0(tgas)\n    @return {rate}\n"

    @pytest.fixture
    def dirs(self, tmp_path, monkeypatch):
        """Config dir with real files; CWD elsewhere holding same-named decoys."""
        cfg_dir, cwd = tmp_path / "cfg", tmp_path / "cwd"
        cfg_dir.mkdir()
        cwd.mkdir()
        (cfg_dir / "net.dat").write_text(self.NETWORK)
        (cfg_dir / "rates.py").write_text(self.TEMPLATE)
        (cfg_dir / "net.jfunc").write_text(self._jfunc(9))
        (cfg_dir / "jaff.toml").write_text('[network.rates]\nT_cutoff = "clip"\n')
        (cwd / "net.jfunc").write_text(self._jfunc(5))
        (cwd / "jaff.toml").write_text("this is [not valid toml\n")
        monkeypatch.chdir(cwd)
        return cfg_dir, cwd

    def _write_cfg(self, cfg_dir, network_block):
        cfg = cfg_dir / "jaffgen.toml"
        cfg.write_text(
            '[jaffgen]\nnetwork_file = "net.dat"\ninput_files = ["rates.py"]\n'
            f'lang = "python"\n[network]\n{network_block}'
        )
        return cfg

    def _state_from(self, cfg):
        from jaff.cli.jaffgen._engine import JaffGen
        from jaff.cli.jaffgen._structs import ResolvedPath, State
        from jaff.drivers import Toml

        jg = JaffGen.__new__(JaffGen)
        jg.state = State()
        jg.state.config_dir = ResolvedPath(cfg.parent, cfg.parent)
        jg.state.config_raw = Toml(cfg)
        jg.set_state_from_config()
        return jg.state.network_args

    def test_funcfile_resolves_against_config_dir(self, dirs):
        cfg_dir, _ = dirs
        sn = self._state_from(self._write_cfg(cfg_dir, 'funcfile = "net.jfunc"\n'))
        assert Path(sn.funcfile) == cfg_dir / "net.jfunc"

    def test_config_resolves_against_config_dir(self, dirs):
        cfg_dir, _ = dirs
        sn = self._state_from(self._write_cfg(cfg_dir, 'config = "jaff.toml"\n'))
        assert Path(sn.config) == cfg_dir / "jaff.toml"

    def test_boolean_funcfile_is_preserved(self, dirs):
        cfg_dir, _ = dirs
        sn = self._state_from(self._write_cfg(cfg_dir, "funcfile = false\n"))
        assert sn.funcfile is False

    def test_absolute_paths_are_kept(self, dirs):
        cfg_dir, cwd = dirs
        block = (
            f'funcfile = "{(cwd / "net.jfunc").as_posix()}"\n'
            f'config = "{(cwd / "jaff.toml").as_posix()}"\n'
        )
        sn = self._state_from(self._write_cfg(cfg_dir, block))
        assert Path(sn.funcfile) == cwd / "net.jfunc"
        assert Path(sn.config) == cwd / "jaff.toml"

    def _generate(self, cfg, *extra):
        from typer.testing import CliRunner

        from jaff.cli.jaffgen._engine import app

        out = cfg.parent / "out"
        args = ["--config", str(cfg), "--outdir", str(out), *extra]
        result = CliRunner().invoke(app, args)
        assert result.exit_code == 0, result.output
        return (out / "rates.py").read_text()

    def test_cli_uses_config_relative_funcfile(self, dirs):
        cfg_dir, _ = dirs
        cfg = self._write_cfg(cfg_dir, 'funcfile = "net.jfunc"\n')
        assert "k[0] = 9" in self._generate(cfg)

    def test_cli_uses_config_relative_network_config(self, dirs):
        # The CWD decoy jaff.toml is malformed; loading it would fail the run.
        # funcfile defaults to the sibling scan, which finds cfg/net.jfunc.
        cfg_dir, _ = dirs
        cfg = self._write_cfg(cfg_dir, 'config = "jaff.toml"\n')
        assert "k[0] = 9" in self._generate(cfg)

    def test_cli_funcfile_override_stays_cwd_relative(self, dirs):
        cfg_dir, _ = dirs
        cfg = self._write_cfg(cfg_dir, 'funcfile = "net.jfunc"\n')
        assert "k[0] = 5" in self._generate(cfg, "--funcfile", "net.jfunc")
