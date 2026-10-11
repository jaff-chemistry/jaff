# ABOUTME: Tests for photoionization database selection (pi_database): props,
# ABOUTME: resolution/fallback in Photochemistry, network wiring, Verner integration.
import logging
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest
import sympy as sp

from jaff import Network
from jaff.common._integrators import get_bounds
from jaff.cli import JaffGen
from jaff.errors import ParserError
from jaff.physics import RadiationProps
from jaff.physics.photo_reactions._photochemistry import Photochemistry
from jaff.physics.photo_reactions._radiation import Radiation
from jaff.plotting._api import _sampled_xsecs

BANDS = [6.0, 13.6, 100.0, 1000.0]


class TestRadiationPropsPiDatabase:
    def test_default_is_norad(self):
        assert RadiationProps(bands=BANDS).pi_database == "norad"

    @pytest.mark.parametrize("value", ["norad", "Leiden", "VERNER"])
    def test_valid_values_lowercased(self, value):
        assert RadiationProps(bands=BANDS, pi_database=value).pi_database == value.lower()

    @pytest.mark.parametrize("value", ["topbase", "", 3, None])
    def test_invalid_values_rejected(self, value):
        with pytest.raises(ParserError, match="pi_database"):
            RadiationProps(bands=BANDS, pi_database=value)

    def test_radiation_stores_pi_database(self):
        with patch("jaff.physics.photo_reactions._radiation.BackgroundField"):
            rad = Radiation(None, RadiationProps(bands=BANDS, pi_database="verner"))
        assert rad.pi_database == "verner"


class TestReactionAttribute:
    def test_reaction_pi_database_defaults_to_none(self, make_network):
        net = make_network("H + H -> H2  []  1e-10\n", name="t.jet")
        assert all(r.pi_database is None for r in net.reactions)


def _net(pi_database=None):
    """Minimal network stand-in: proxy flag + optional radiation."""
    rad = None if pi_database is None else SimpleNamespace(pi_database=pi_database)
    return SimpleNamespace(_use_proxy_photoreaction=False, radiation=rad)


def _rxn(key, pi_database=None):
    return SimpleNamespace(
        serialized=key,
        normalized_proxy_reaction_str=lambda: key,
        pi_database=pi_database,
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
        x = make(_net()).get_xsec(_rxn("He._PHOTON__He+.e-"))
        assert x["database"] == "norad"
        assert x["photodecay"] is not None and x["photodecay_expr"] is None
        warn.assert_not_called()

    @pytest.mark.parametrize("db", ["norad", "verner", "leiden"])
    def test_global_choice(self, pc_and_warn, db):
        make, warn = pc_and_warn
        x = make(_net(db)).get_xsec(_rxn("C._PHOTON__C+.e-"))
        assert x["database"] == db
        warn.assert_not_called()

    def test_verner_dict_is_symbolic(self, pc_and_warn):
        make, _ = pc_and_warn
        x = make(_net("verner")).get_xsec(_rxn("C._PHOTON__C+.e-"))
        assert isinstance(x["photodecay_expr"], sp.Basic)
        assert x["photon_energy"] is None and x["photodecay"] is None
        assert x["_equations"] == {"pa": False, "decay_type": "ionization"}

    def test_per_reaction_override_beats_global(self, pc_and_warn):
        make, _ = pc_and_warn
        x = make(_net("verner")).get_xsec(_rxn("C._PHOTON__C+.e-", "leiden"))
        assert x["database"] == "leiden"

    def test_leiden_missing_falls_back_to_norad(self, pc_and_warn):
        make, warn = pc_and_warn
        pc = make(_net("leiden"))
        x = pc.get_xsec(_rxn("P+._PHOTON__P++.e-"))
        assert x["database"] == "norad"
        warn.assert_not_called()
        assert pc.fallbacks == [("P+._PHOTON__P++.e-", "leiden", "norad")]

    def test_norad_missing_falls_back_to_verner(self, pc_and_warn):
        make, warn = pc_and_warn
        pc = make(_net())
        x = pc.get_xsec(_rxn("Ca+++++._PHOTON__Ca++++++.e-"))
        assert x["database"] == "verner"
        warn.assert_not_called()
        assert pc.fallbacks == [("Ca+++++._PHOTON__Ca++++++.e-", "norad", "verner")]

    def test_molecule_falls_back_to_leiden(self, pc_and_warn):
        make, warn = pc_and_warn
        pc = make(_net())
        x = pc.get_xsec(_rxn("CO._PHOTON__CO+.e-"))
        assert x["database"] == "leiden"
        warn.assert_not_called()
        assert pc.fallbacks == [("CO._PHOTON__CO+.e-", "norad", "leiden")]

    def test_norad_never_reports_leiden_absorption(self, pc_and_warn):
        make, _ = pc_and_warn
        x = make(_net("norad")).get_xsec(_rxn("C._PHOTON__C+.e-"))
        assert x["_equations"]["pa"] is False and x["photo_absorption"] is None

    def test_missing_everywhere_returns_none_when_not_required(self, pc_and_warn):
        make, _ = pc_and_warn
        assert make(_net()).get_xsec(_rxn("C9H9._PHOTON__C9H9+.e-")) is None

    def test_missing_everywhere_raises_when_required(self, pc_and_warn):
        make, _ = pc_and_warn
        with pytest.raises(ParserError, match="C9H9"):
            make(_net()).get_xsec(_rxn("C9H9._PHOTON__C9H9+.e-"), required=True)

    def test_dissociation_ignores_pi_database(self, pc_and_warn):
        make, warn = pc_and_warn
        x = make(_net("verner")).get_xsec(_rxn("CO._PHOTON__C.O"))
        assert x["database"] == "leiden"
        assert x["_equations"]["decay_type"] == "dissociation"
        warn.assert_not_called()

    def test_override_on_dissociation_raises(self, pc_and_warn):
        make, _ = pc_and_warn
        with pytest.raises(ParserError, match="photoionization"):
            make(_net()).get_xsec(_rxn("CO._PHOTON__C.O", "verner"))


class _FakeReaction:
    """Hashable stand-in for the Reaction attributes Radiation touches."""

    def __init__(self, xsecs):
        self.xsecs_dict = xsecs
        self.dRad = sp.Float(0.0)
        self._metadata = {}
        self.rad_groups = []
        self.rate = None
        self.rad_xsecs = None


class TestVernerIntegration:
    def _rad(self, bands=BANDS):
        with patch("jaff.physics.photo_reactions._radiation.BackgroundField"):
            return Radiation(
                None,
                RadiationProps(bands=list(bands), profile_index=[2.0, 1.0, 0.5], c=1.0),
            )

    @staticmethod
    def _reaction(xsecs):
        return _FakeReaction(xsecs)

    @pytest.mark.parametrize("bands", [BANDS, [6.0, 13.6, 100.0, "inf"]])
    def test_symbolic_matches_dense_tabulation(self, pc_and_warn, bands):
        make, _ = pc_and_warn
        sym = make(_net("verner")).get_xsec(_rxn("C._PHOTON__C+.e-"))
        E = sp.Symbol("E")
        # The fit is zero above E_max, so a dense grid up to E_max suffices.
        e_max = float(get_bounds(sym["photodecay_expr"], E)[-1])
        grid = np.geomspace(11.0, e_max, 400_000)
        tab = {
            **sym,
            "database": "verner",
            "photon_energy": grid,
            "photodecay": np.asarray(
                sp.lambdify(E, sym["photodecay_expr"], "numpy")(grid), dtype=float
            ),
            "photodecay_expr": None,
        }
        rad = self._rad(bands)
        r_sym, r_tab = self._reaction(sym), self._reaction(tab)
        rad.set_reaction_rate_coefficient(r_sym)
        rad.set_reaction_rate_coefficient(r_tab)
        for g in rad.groups:
            a, b = float(g.props[r_sym]["xsec"]), float(g.props[r_tab]["xsec"])
            assert a == pytest.approx(b, rel=1e-4, abs=1e-30)
        assert float(r_sym.rad_xsecs) == pytest.approx(float(r_tab.rad_xsecs), rel=1e-4)


PHOTO_JET = (
    "C -> C+ + E  []  PHOTO, 11.3\n"
    "Ca+++++ -> Ca++++++ + E  []  PHOTO, 100\n"
    "CO -> C + O  []  PHOTO, 11\n"
)


def _by_key(net):
    return {r.serialized: r for r in net.reactions}


class TestNetworkWiring:
    def test_global_choice_reaches_reactions(self, make_network):
        net = make_network(
            PHOTO_JET,
            name="t.jet",
            radiation_props=RadiationProps(bands=BANDS, pi_database="verner"),
        )
        r = _by_key(net)
        assert r["C._PHOTON__C+.e-"].xsecs_dict["database"] == "verner"
        assert r["CO._PHOTON__C.O"].xsecs_dict["database"] == "leiden"

    def test_jaff_toml_per_reaction_override(self, make_network, tmp_path):
        cfg = tmp_path / "jaff.toml"
        cfg.write_text('[network.reactions."C._PHOTON__C+.e-"]\npi_database = "Leiden"\n')
        net = make_network(
            PHOTO_JET,
            name="t.jet",
            config=str(cfg),
            radiation_props=RadiationProps(bands=BANDS),
        )
        rxn = _by_key(net)["C._PHOTON__C+.e-"]
        assert rxn.pi_database == "leiden"
        assert rxn.xsecs_dict["database"] == "leiden"

    def test_override_on_non_photo_reaction_raises(self, make_network, tmp_path):
        cfg = tmp_path / "jaff.toml"
        cfg.write_text('[network.reactions."H.H__H2"]\npi_database = "verner"\n')
        with pytest.raises(ParserError, match="photo"):
            make_network("H + H -> H2  []  1e-10\n", name="t.jet", config=str(cfg))

    def test_invalid_override_value_raises(self, make_network, tmp_path):
        cfg = tmp_path / "jaff.toml"
        cfg.write_text('[network.reactions."C._PHOTON__C+.e-"]\npi_database = "x"\n')
        with pytest.raises(ParserError, match="pi_database"):
            make_network(PHOTO_JET, name="t.jet", config=str(cfg))

    def test_empty_override_value_raises(self, make_network, tmp_path):
        cfg = tmp_path / "jaff.toml"
        cfg.write_text('[network.reactions."C._PHOTON__C+.e-"]\npi_database = ""\n')
        with pytest.raises(ParserError, match="pi_database"):
            make_network(PHOTO_JET, name="t.jet", config=str(cfg))

    def test_override_keyed_by_written_reaction_under_proxy(self, make_network, tmp_path):
        cfg = tmp_path / "jaff.toml"
        cfg.write_text(
            '[network.reactions."Cx._PHOTON__C+.e-"]\npi_database = "verner"\n'
        )
        net = make_network(
            "Cx -> C+ + E  []  PHOTO, 11.3\n",
            name="p.jet",
            config=str(cfg),
            use_proxy_photoreaction=True,
            radiation_props=RadiationProps(bands=BANDS),
        )
        rxn = _by_key(net)["Cx._PHOTON__C+.e-"]
        assert rxn.xsecs_dict["database"] == "verner"

    def test_missing_xsec_raises_only_with_radiation(self, make_network):
        line = "C9H9 -> C9H9+ + E  []  PHOTO, 9\n"
        net = make_network(line, name="a.jet")  # no radiation: fine
        assert _by_key(net)["C9H9._PHOTON__C9H9+.e-"].xsecs_dict is None
        with pytest.raises(ParserError, match="C9H9"):
            make_network(line, name="b.jet", radiation_props=RadiationProps(bands=BANDS))


class TestFallbackSummary:
    def test_one_warning_lists_all_fallbacks(self, pc_and_warn):
        make, warn = pc_and_warn
        pc = make(_net())
        for key in (
            "CO._PHOTON__CO+.e-",
            "Ca+++++._PHOTON__Ca++++++.e-",
            "C._PHOTON__C+.e-",
        ):
            pc.get_xsec(_rxn(key))
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
        pc.get_xsec(_rxn("C._PHOTON__C+.e-"))
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
            make_network(lines, name="t.jet", radiation_props=RadiationProps(bands=BANDS))
        warn.assert_called_once()
        msg = warn.call_args.args[0]
        assert "CO._PHOTON__CO+.e-" in msg and "Ca+++++._PHOTON__Ca++++++.e-" in msg


REPO = Path(__file__).resolve().parent.parent
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
            errors=None,
            network_config=None,
            outdir=str(tmp_path / "out"),
            indir=None,
            files=None,
            template="microphysics",
            lang="cxx",
        )
    )


class TestJaffgenWiring:
    RAD = "[network.radiation]\nbands = [6, 13.6, 100]\nuse_proxy_photoreaction = true\n"

    def test_global_pi_database(self, tmp_path):
        gen = _jaffgen(tmp_path, self.RAD + 'pi_database = "verner"\n')
        assert gen.net.radiation.pi_database == "verner"
        (rxn,) = [r for r in gen.net.reactions if r.type == "photo"]
        assert rxn.xsecs_dict["database"] == "verner"

    def test_per_reaction_pi_database(self, tmp_path):
        gen = _jaffgen(
            tmp_path,
            self.RAD
            + '\n[network.reactions."H._PHOTON__H+.e-"]\npi_database = "leiden"\n',
        )
        (rxn,) = [r for r in gen.net.reactions if r.type == "photo"]
        assert rxn.pi_database == "leiden"
        assert rxn.xsecs_dict["database"] == "leiden"


class TestJaffRoundTrip:
    def test_override_survives_round_trip(self, make_network, tmp_path):
        cfg = tmp_path / "jaff.toml"
        cfg.write_text('[network.reactions."C._PHOTON__C+.e-"]\npi_database = "verner"\n')
        props = RadiationProps(bands=BANDS)
        net = make_network(
            PHOTO_JET, name="t.jet", config=str(cfg), radiation_props=props
        )
        out = tmp_path / "t.jaff"
        net.to_jaff(out)
        back = Network(str(out), radiation_props=props)
        rxn = _by_key(back)["C._PHOTON__C+.e-"]
        assert rxn.pi_database == "verner"
        assert rxn.xsecs_dict["database"] == "verner"
        assert _by_key(back)["Ca+++++._PHOTON__Ca++++++.e-"].pi_database is None


class TestPlotSampling:
    def test_expression_sampled_on_threshold_range(self, pc_and_warn):
        make, _ = pc_and_warn
        x = make(_net("verner")).get_xsec(_rxn("C._PHOTON__C+.e-"))
        s = _sampled_xsecs(x)
        E = s["photon_energy"]
        assert E[0] == pytest.approx(11.26, rel=1e-3)
        assert np.all(np.diff(E) > 0) and np.all(s["photodecay"] > 0)

    def test_tabulated_passthrough(self, pc_and_warn):
        make, _ = pc_and_warn
        x = make(_net("norad")).get_xsec(_rxn("C._PHOTON__C+.e-"))
        assert _sampled_xsecs(x) is x

    def test_plotter_plot_xsec_accepts_verner(self, pc_and_warn):
        import matplotlib

        matplotlib.use("Agg")
        from jaff.plotting.plotter import Plotter

        make, _ = pc_and_warn
        x = make(_net("verner")).get_xsec(_rxn("C._PHOTON__C+.e-"))
        fig, _ax = Plotter().plot_xsec(x, show=False)
        assert fig is not None
