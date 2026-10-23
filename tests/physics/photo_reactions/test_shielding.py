# ABOUTME: Photo-reaction line shielding: the configured shielding model multiplies
# ABOUTME: the band rates; each shielding plugin validates its own config keys

from unittest.mock import patch

import pytest
from pydantic import ValidationError
from sympy import Float

from jaff.physics.photo_reactions.shielding.config import ShieldingParams
from jaff.errors import ParserError
from jaff.physics.photo_reactions.config import RadiationConfig
from jaff.physics.photo_reactions.shielding import ShieldingFunction
from jaff.physics.photo_reactions.shielding.global_.leiden import Leiden
from jaff.physics.photo_reactions.shielding.H2__PHOTON__H_H.db1996 import DB1996
from jaff.physics.photo_reactions.shielding.H2__PHOTON__H_H.hg2015 import HG2015

H2_PHOTO = "H2 -> H + H  []  PHOTO, 11.2\n"
H2_KEY = "H2._PHOTON__H.H"


@pytest.fixture
def shielded_network(make_network):
    """H2 photodissociation shielded by *shielding*; returns its band rates."""

    def _make(shielding):
        own = {"photo": {"shielding": shielding}}
        config = {"network": {"reactions": {"reaction": {H2_KEY: own}}}}
        with patch("jaff.physics.photo_reactions._radiation.BackgroundField"):
            net = make_network(
                H2_PHOTO,
                name="h2.jet",
                funcfile=False,
                radiation=RadiationConfig(bands=[6.0, 13.6]),
                config=config,
            )

        rxn = net.reactions[H2_KEY, "photo"]
        assert rxn.photo is not None and rxn.photo.bands is not None

        return [band.k for band in rxn.photo.bands]

    return _make


def test_value_key_is_rejected(shielded_network):
    # "value" was once a cached factor that, set by hand, bypassed the model.
    with pytest.raises(ParserError, match="value"):
        shielded_network({"type": "hg2015", "value": 0.5})


class TestConfigModel:
    def test_base_default_accepts_no_keys(self):
        class Plain(ShieldingFunction):
            name = "plain"

            def get_shielding(self, reaction, network):
                return Float(1.0)

        model = Plain.config_model()

        assert model.model_validate({}) == model()
        with pytest.raises(ValidationError, match="x"):
            model.model_validate({"x": 1})

    @pytest.mark.parametrize("cls", [Leiden, HG2015, DB1996])
    def test_plugins_expose_shielding_params(self, cls):
        assert issubclass(cls.config_model(), ShieldingParams)

    def test_plugins_do_not_share_a_model(self):
        models = [cls.config_model() for cls in (Leiden, HG2015, DB1996)]

        assert len(set(models)) == len(models)


class TestLeidenParams:
    def test_requires_shielded_by(self):
        with pytest.raises(ValidationError, match="shielded_by"):
            Leiden.config_model().model_validate({})

    def test_rejects_unknown_species(self):
        with pytest.raises(ValidationError, match="shielded_by"):
            Leiden.config_model().model_validate({"shielded_by": ["O2"]})

    def test_radiation_defaults_to_isrf(self):
        params = Leiden.config_model().model_validate({"shielded_by": ["self"]})

        assert params.radiation == "isrf"

    def test_radiation_is_lower_cased(self):
        params = Leiden.config_model().model_validate(
            {"shielded_by": ["self", "H2"], "radiation": "Draine"}
        )

        assert params.radiation == "draine"


@pytest.mark.parametrize("cls", [HG2015, DB1996])
class TestH2Params:
    def test_floor_defaults(self, cls):
        params = cls.config_model().model_validate({})

        assert params.min_ncol == pytest.approx(1e-50)
        assert params.min_vdisp == pytest.approx(1e-50)

    def test_accepts_int(self, cls):
        params = cls.config_model().model_validate({"min_ncol": 1})

        assert params.min_ncol == pytest.approx(1.0)

    def test_rejects_string(self, cls):
        with pytest.raises(ValidationError, match="min_ncol"):
            cls.config_model().model_validate({"min_ncol": "1e-35"})

    def test_rejects_unknown_key(self, cls):
        with pytest.raises(ValidationError, match="shielded_by"):
            cls.config_model().model_validate({"shielded_by": ["self"]})


class TestNetworkShielding:
    def test_invalid_setting_is_parser_error(self, shielded_network):
        with pytest.raises(ParserError, match="min_ncol"):
            shielded_network({"type": "hg2015", "min_ncol": "abc"})

    def test_unknown_setting_is_parser_error(self, shielded_network):
        with pytest.raises(ParserError, match="foo"):
            shielded_network({"type": "hg2015", "foo": 1.0})


CO_PHOTO = "CO -> C + O  []  PHOTO, 11.1\n"
CO_KEY = "CO._PHOTON__C.O"
LEIDEN = {"type": "leiden", "shielded_by": ["self", "H2"]}


@pytest.fixture
def leiden_network(make_network):
    """CO photodissociation with Leiden shielding, built without jaffgen."""
    config = {
        "network": {"reactions": {"reaction": {CO_KEY: {"photo": {"shielding": LEIDEN}}}}}
    }
    with patch("jaff.physics.photo_reactions._radiation.BackgroundField"):
        return make_network(
            CO_PHOTO,
            name="co.jet",
            funcfile=False,
            radiation=RadiationConfig(bands=[6.0, 13.6]),
            config=config,
        )


class TestLeidenTables:
    def test_table_is_kept_on_the_reaction(self, leiden_network):
        tables = leiden_network.reactions[CO_KEY, "photo"].photo.shielding_tables

        (name,) = tables
        assert name == f"shielding_{CO_KEY}.hdf5"
        assert {"N", "self", "H2"} <= set(tables[name])

    def test_network_writes_the_tables(self, leiden_network, tmp_path):
        written = leiden_network.write_shielding_tables(tmp_path)

        assert written == [tmp_path / f"shielding_{CO_KEY}.hdf5"]
        assert written[0].exists()

    def test_network_without_shielding_writes_nothing(self, make_network, tmp_path):
        net = make_network("H + H -> H2 [10,1000] 1e-10\n", name="n.dat", funcfile=False)

        assert net.write_shielding_tables(tmp_path) == []


def test_jaffgen_writes_shielding_tables_to_output(tmp_path):
    from typer.testing import CliRunner

    from jaff.cli.jaffgen._engine import app

    (tmp_path / "co.jet").write_text(CO_PHOTO)
    (tmp_path / "rates.py").write_text(
        "$$JAFF REPEAT idx, rate IN rates\nk[$idx$] = $rate$\n$$JAFF END\n"
    )
    cfg = tmp_path / "jaffgen.toml"
    cfg.write_text(
        '[jaffgen]\nnetwork_file = "co.jet"\ninput_files = ["rates.py"]\n'
        "[network]\nfuncfile = false\n"
        "[network.radiation]\nbands = [6.0, 13.6]\n"
        f'[network.reactions.reaction."{CO_KEY}".photo.shielding]\n'
        'shielded_by = ["self"]\n'
    )

    out = tmp_path / "out"
    with patch("jaff.physics.photo_reactions._radiation.BackgroundField"):
        result = CliRunner().invoke(app, ["--config", str(cfg), "--outdir", str(out)])

    assert result.exit_code == 0, result.output
    assert (out / f"shielding_{CO_KEY}.hdf5").exists()
