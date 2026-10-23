# ABOUTME: Tests for photoionization database selection (pi_database): props,
# ABOUTME: resolution/fallback in Photochemistry, network wiring, Verner integration.
import logging
from dataclasses import replace
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest
import sympy as sp

from jaff import Network
from jaff.core.network.config import NetworkConfig
from jaff.cli import JaffGen
from jaff.common._integrators import get_bounds
from jaff.core.network import NetworkSymbols
from jaff.errors import ParserError
from jaff.physics.photo_reactions.config import RadiationConfig
from jaff.physics.photo_reactions._photochemistry import Photochemistry, XsecProcess
from jaff.physics.photo_reactions._radiation import Radiation
from jaff.plotting.adapters.xsecs import _sampled_xsecs
from tests._paths import REPO_ROOT
from tests.xsec_props import stub_photo

BANDS = [6.0, 13.6, 100.0, 1000.0]


class TestReactionAttribute:
    def test_reaction_pi_database_defaults_to_global(self, make_network):
        net = make_network(
            "C -> C+ + E  []  PHOTO, 11.3\nH + H -> H2  []  1e-10\n", name="t.jet"
        )
        r = _by_key(net)
        assert r["C._PHOTON__C+.e-"].config.photo.pi_database == "norad"
        assert r["H.H__H2"].photo is None


def _net(pi_database=None, own=None):
    """Minimal network stand-in: config (no proxy lookups, global pi_database).

    *own* maps a serialized reaction to the ``pi_database`` set for it alone.
    """
    photo = {"use_proxy_photoreaction": False}
    if pi_database is not None:
        photo["pi_database"] = pi_database

    entries = {key: {"photo": {"pi_database": db}} for key, db in (own or {}).items()}
    config = NetworkConfig.model_validate(
        {"reactions": {"photo": photo, "reaction": entries}}
    )

    return SimpleNamespace(config=config, radiation=None)


def _rxn(key, net):
    """Reaction stand-in with its settings resolved from *net*'s config."""
    return SimpleNamespace(
        serialized=key,
        normalized_proxy_reaction_str=lambda: key,
        photo=stub_photo(config=net.config.reactions.resolve(key).photo),
    )


@pytest.fixture
def pc_and_warn():
    """Photochemistry factory whose logger.warning is a mock."""
    logger = logging.getLogger("pi-database-test")

    def _make(net):
        with patch("jaff.physics.photo_reactions._photochemistry.JaffLogger") as jl:
            jl.return_value.get_logger.return_value = logger
            return Photochemistry(net)

    with patch.object(logger, "warning") as warn:
        yield _make, warn


class TestResolution:
    def test_default_norad_without_radiation(self, pc_and_warn):
        make, warn = pc_and_warn
        x = (pc := make(_net())).get_xsec(_rxn("He._PHOTON__He+.e-", pc.net))
        assert x.database == "norad"
        assert x.photodecay is not None and x.photodecay_expr is None
        warn.assert_not_called()

    @pytest.mark.parametrize("db", ["norad", "verner", "leiden"])
    def test_global_choice(self, pc_and_warn, db):
        make, warn = pc_and_warn
        x = (pc := make(_net(db))).get_xsec(_rxn("C._PHOTON__C+.e-", pc.net))
        assert x.database == db
        warn.assert_not_called()

    def test_verner_dict_is_symbolic(self, pc_and_warn):
        make, _ = pc_and_warn
        x = (pc := make(_net("verner"))).get_xsec(_rxn("C._PHOTON__C+.e-", pc.net))
        assert isinstance(x.photodecay_expr, sp.Basic)
        assert x.photodecay_expr.free_symbols == {NetworkSymbols.ph_E}
        assert x.photon_energy is None and x.photodecay is None
        assert x.process == XsecProcess(absorption=False, decay_type="ionization")

    def test_per_reaction_override_beats_global(self, pc_and_warn):
        make, _ = pc_and_warn
        x = (pc := make(_net("verner", own={"C._PHOTON__C+.e-": "leiden"}))).get_xsec(_rxn("C._PHOTON__C+.e-", pc.net))
        assert x.database == "leiden"

    def test_leiden_missing_falls_back_to_norad(self, pc_and_warn):
        make, warn = pc_and_warn
        pc = make(_net("leiden"))
        x = pc.get_xsec(_rxn("P+._PHOTON__P++.e-", pc.net))
        assert x.database == "norad"
        warn.assert_not_called()
        assert pc.fallbacks == [("P+._PHOTON__P++.e-", "leiden", "norad")]

    def test_norad_missing_falls_back_to_verner(self, pc_and_warn):
        make, warn = pc_and_warn
        pc = make(_net())
        x = pc.get_xsec(_rxn("Ca+++++._PHOTON__Ca++++++.e-", pc.net))
        assert x.database == "verner"
        warn.assert_not_called()
        assert pc.fallbacks == [("Ca+++++._PHOTON__Ca++++++.e-", "norad", "verner")]

    def test_molecule_falls_back_to_leiden(self, pc_and_warn):
        make, warn = pc_and_warn
        pc = make(_net())
        x = pc.get_xsec(_rxn("CO._PHOTON__CO+.e-", pc.net))
        assert x.database == "leiden"
        warn.assert_not_called()
        assert pc.fallbacks == [("CO._PHOTON__CO+.e-", "norad", "leiden")]

    def test_norad_never_reports_leiden_absorption(self, pc_and_warn):
        make, _ = pc_and_warn
        x = (pc := make(_net("norad"))).get_xsec(_rxn("C._PHOTON__C+.e-", pc.net))
        assert x.process.absorption is False and x.photo_absorption is None

    def test_missing_everywhere_returns_none_when_not_required(self, pc_and_warn):
        make, _ = pc_and_warn
        assert (pc := make(_net())).get_xsec(_rxn("C9H9._PHOTON__C9H9+.e-", pc.net)) is None

    def test_missing_everywhere_raises_when_required(self, pc_and_warn):
        make, _ = pc_and_warn
        with pytest.raises(ParserError, match="C9H9"):
            (pc := make(_net())).get_xsec(_rxn("C9H9._PHOTON__C9H9+.e-", pc.net), required=True)

    def test_dissociation_ignores_pi_database(self, pc_and_warn):
        make, warn = pc_and_warn
        x = (pc := make(_net("verner"))).get_xsec(_rxn("CO._PHOTON__C.O", pc.net))
        assert x.database == "leiden"
        assert x.process.decay_type == "dissociation"
        warn.assert_not_called()

    def test_override_on_dissociation_raises(self, pc_and_warn):
        make, _ = pc_and_warn
        with pytest.raises(ParserError, match="photoionization"):
            (pc := make(_net(own={"CO._PHOTON__C.O": "verner"}))).get_xsec(_rxn("CO._PHOTON__C.O", pc.net))


class _FakeReaction:
    """Stand-in for the Reaction attributes Radiation touches."""

    def __init__(self, xsecs):
        self.photo = stub_photo(xsec=xsecs)


class TestVernerIntegration:
    def _rad(self, bands=BANDS):
        with patch("jaff.physics.photo_reactions._radiation.BackgroundField"):
            return Radiation(
                None,
                RadiationConfig(bands=list(bands), profile_index=[2.0, 1.0, 0.5], c=1.0),
            )

    @staticmethod
    def _reaction(xsecs):
        return _FakeReaction(xsecs)

    @pytest.mark.parametrize("bands", [BANDS, [6.0, 13.6, 100.0, "inf"]])
    def test_symbolic_matches_dense_tabulation(self, pc_and_warn, bands):
        make, _ = pc_and_warn
        sym = (pc := make(_net("verner"))).get_xsec(_rxn("C._PHOTON__C+.e-", pc.net))
        E = NetworkSymbols.ph_E
        # The fit is zero above E_max, so a dense grid up to E_max suffices.
        e_max = float(get_bounds(sym.photodecay_expr, E)[-1])
        grid = np.geomspace(11.0, e_max, 400_000)
        tab = replace(
            sym,
            photon_energy=grid,
            photodecay=np.asarray(
                sp.lambdify(E, sym.photodecay_expr, "numpy")(grid), dtype=float
            ),
            photodecay_expr=None,
        )
        rad = self._rad(bands)
        r_sym, r_tab = self._reaction(sym), self._reaction(tab)
        _, tot_sym, grouped_sym = rad.get_reaction_rate_props(r_sym)
        _, tot_tab, grouped_tab = rad.get_reaction_rate_props(r_tab)
        for g_sym, g_tab in zip(grouped_sym, grouped_tab, strict=True):
            a, b = float(g_sym.xsec), float(g_tab.xsec)
            assert a == pytest.approx(b, rel=1e-4, abs=1e-30)
        assert float(tot_sym) == pytest.approx(float(tot_tab), rel=1e-4)


PHOTO_JET = (
    "C -> C+ + E  []  PHOTO, 11.3\n"
    "Ca+++++ -> Ca++++++ + E  []  PHOTO, 100\n"
    "CO -> C + O  []  PHOTO, 11\n"
)


def _by_key(net):
    return {r.serialized: r for r in net.reactions}


def _global_photo(**photo):
    """Config layer setting the global photo settings ([network.reactions.photo])."""
    return {"network": {"reactions": {"photo": photo}}}


class TestNetworkWiring:
    def test_global_choice_reaches_reactions(self, make_network):
        net = make_network(
            PHOTO_JET,
            name="t.jet",
            radiation=RadiationConfig(bands=BANDS),
            config=_global_photo(pi_database="verner"),
        )
        r = _by_key(net)
        assert r["C._PHOTON__C+.e-"].photo.xsec.database == "verner"
        assert r["CO._PHOTON__C.O"].photo.xsec.database == "leiden"

    def test_jaff_toml_per_reaction_override(self, make_network, tmp_path):
        cfg = tmp_path / "jaff.toml"
        cfg.write_text(
            '[network.reactions.reaction."C._PHOTON__C+.e-".photo]\npi_database = "Leiden"\n'
        )
        net = make_network(
            PHOTO_JET,
            name="t.jet",
            config=str(cfg),
            radiation=RadiationConfig(bands=BANDS),
        )
        rxn = _by_key(net)["C._PHOTON__C+.e-"]
        assert rxn.config.photo.pi_database == "leiden"
        assert rxn.photo.xsec.database == "leiden"

    def test_override_on_non_photo_reaction_raises(self, make_network, tmp_path):
        cfg = tmp_path / "jaff.toml"
        cfg.write_text(
            '[network.reactions.reaction."H.H__H2".photo]\npi_database = "verner"\n'
        )
        with pytest.raises(ParserError, match="photo"):
            make_network("H + H -> H2  []  1e-10\n", name="t.jet", config=str(cfg))

    def test_invalid_override_value_raises(self, make_network, tmp_path):
        cfg = tmp_path / "jaff.toml"
        cfg.write_text(
            '[network.reactions.reaction."C._PHOTON__C+.e-".photo]\npi_database = "x"\n'
        )
        with pytest.raises(ParserError, match="pi_database"):
            make_network(PHOTO_JET, name="t.jet", config=str(cfg))

    def test_empty_override_value_raises(self, make_network, tmp_path):
        cfg = tmp_path / "jaff.toml"
        cfg.write_text(
            '[network.reactions.reaction."C._PHOTON__C+.e-".photo]\npi_database = ""\n'
        )
        with pytest.raises(ParserError, match="pi_database"):
            make_network(PHOTO_JET, name="t.jet", config=str(cfg))

    def test_override_keyed_by_written_reaction_under_proxy(self, make_network, tmp_path):
        cfg = tmp_path / "jaff.toml"
        cfg.write_text(
            '[network.reactions.reaction."Cx._PHOTON__C+.e-".photo]\npi_database = "verner"\n'
        )
        net = make_network(
            "Cx -> C+ + E  []  PHOTO, 11.3\n",
            name="p.jet",
            config=str(cfg),
            reactions={"photo": {"use_proxy_photoreaction": True}},
            radiation=RadiationConfig(bands=BANDS),
        )
        rxn = _by_key(net)["Cx._PHOTON__C+.e-"]
        assert rxn.photo.xsec.database == "verner"

    def test_missing_xsec_raises_only_with_radiation(self, make_network):
        line = "C9H9 -> C9H9+ + E  []  PHOTO, 9\n"
        net = make_network(line, name="a.jet")  # no radiation: fine
        assert _by_key(net)["C9H9._PHOTON__C9H9+.e-"].photo.xsec is None
        with pytest.raises(ParserError, match="C9H9"):
            make_network(line, name="b.jet", radiation=RadiationConfig(bands=BANDS))


class TestFallbackSummary:
    def test_one_warning_lists_all_fallbacks(self, pc_and_warn):
        make, warn = pc_and_warn
        pc = make(_net())
        for key in (
            "CO._PHOTON__CO+.e-",
            "Ca+++++._PHOTON__Ca++++++.e-",
            "C._PHOTON__C+.e-",
        ):
            pc.get_xsec(_rxn(key, pc.net))
        warn.assert_not_called()
        pc.log_fallback_summary()
        warn.assert_called_once()
        msg = warn.call_args.args[0]
        assert "norad -> leiden: CO._PHOTON__CO+.e-" in msg
        assert "norad -> verner: Ca+++++._PHOTON__Ca++++++.e-" in msg
        assert "C._PHOTON__C+.e-" not in msg  # NORAD had it: no fallback
        assert pc.fallbacks == []  # cleared after logging

    def test_no_fallbacks_no_warning(self, pc_and_warn):
        make, warn = pc_and_warn
        pc = make(_net())
        pc.get_xsec(_rxn("C._PHOTON__C+.e-", pc.net))
        pc.log_fallback_summary()
        warn.assert_not_called()

    def test_network_load_emits_single_summary(self, make_network):
        logger = logging.getLogger("pi-database-net-test")
        lines = PHOTO_JET + "CO -> CO+ + E  []  PHOTO, 14\n"
        with (
            patch("jaff.physics.photo_reactions._photochemistry.JaffLogger") as jl,
            patch.object(logger, "warning") as warn,
        ):
            jl.return_value.get_logger.return_value = logger
            make_network(
                lines, name="t.jet", radiation=RadiationConfig(bands=BANDS)
            )
        warn.assert_called_once()
        msg = warn.call_args.args[0]
        assert "CO._PHOTON__CO+.e-" in msg and "Ca+++++._PHOTON__Ca++++++.e-" in msg


REPO = REPO_ROOT
H_PHOTO = REPO / "networks" / "h_photoionization" / "h_photo.jet"


def _jaffgen(tmp_path, toml):
    cfg = tmp_path / "jaffgen.toml"
    cfg.write_text(toml)
    return JaffGen(
        SimpleNamespace(
            network=str(H_PHOTO),
            config=str(cfg),
            label=None,
            funcfile=None,
            duplicate_policy=None,
            expand_nuclei=None,
            strict=None,
            network_config=None,
            outdir=str(tmp_path / "out"),
            indir=None,
            files=None,
            template="microphysics",
        )
    )


class TestJaffgenWiring:
    RAD = "[network.radiation]\nbands = [6, 13.6, 100]\n"
    PROXY = "[network.reactions.photo]\nuse_proxy_photoreaction = true\n"

    def test_global_pi_database(self, tmp_path):
        gen = _jaffgen(tmp_path, self.RAD + self.PROXY + 'pi_database = "verner"\n')
        assert gen.net.config.reactions.photo.pi_database == "verner"
        (rxn,) = [r for r in gen.net.reactions if r.type == "photo"]
        assert rxn.photo.xsec.database == "verner"

    def test_per_reaction_pi_database(self, tmp_path):
        gen = _jaffgen(
            tmp_path,
            self.RAD
            + self.PROXY
            + '\n[network.reactions.reaction."H._PHOTON__H+.e-".photo]\npi_database = "leiden"\n',
        )
        (rxn,) = [r for r in gen.net.reactions if r.type == "photo"]
        assert rxn.config.photo.pi_database == "leiden"
        assert rxn.photo.xsec.database == "leiden"


class TestJaffRoundTrip:
    def test_override_survives_round_trip(self, make_network, tmp_path):
        cfg = tmp_path / "jaff.toml"
        cfg.write_text(
            '[network.reactions.reaction."C._PHOTON__C+.e-".photo]\npi_database = "verner"\n'
        )
        props = RadiationConfig(bands=BANDS)
        net = make_network(
            PHOTO_JET, name="t.jet", config=str(cfg), radiation=props
        )
        out = tmp_path / "t.jaff"
        net.to_jaff(out)
        back = Network(str(out), radiation=props)
        rxn = _by_key(back)["C._PHOTON__C+.e-"]
        assert rxn.config.photo.pi_database == "verner"
        assert rxn.photo.xsec.database == "verner"
        assert _by_key(back)["Ca+++++._PHOTON__Ca++++++.e-"].config.photo.pi_database == "norad"


class TestPlotSampling:
    def test_expression_sampled_on_threshold_range(self, pc_and_warn):
        make, _ = pc_and_warn
        x = (pc := make(_net("verner"))).get_xsec(_rxn("C._PHOTON__C+.e-", pc.net))
        s = _sampled_xsecs(x)
        E = s.photon_energy
        assert E[0] == pytest.approx(11.26, rel=1e-3)
        assert np.all(np.diff(E) > 0) and np.all(s.photodecay > 0)

    def test_tabulated_passthrough(self, pc_and_warn):
        make, _ = pc_and_warn
        x = (pc := make(_net("norad"))).get_xsec(_rxn("C._PHOTON__C+.e-", pc.net))
        assert _sampled_xsecs(x) is x

    def test_xsec_series_accepts_verner(self, pc_and_warn):
        from types import SimpleNamespace

        from jaff.plotting.adapters import xsec_series

        make, _ = pc_and_warn
        x = (pc := make(_net("verner"))).get_xsec(_rxn("C._PHOTON__C+.e-", pc.net))
        (s,) = xsec_series(
            SimpleNamespace(photo=stub_photo(xsec=x)), processes="photodecay"
        )
        assert s["x"][0] == pytest.approx(11.26, rel=1e-3)
        assert np.all(s["y"] > 0)
