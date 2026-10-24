# ABOUTME: Photo-reaction line shielding: the configured shielding model multiplies
# ABOUTME: the band rates; each shielding plugin validates its own config keys

from unittest.mock import patch

import pytest
from pydantic import ValidationError
from sympy import Float

from jaff.physics.photo_reactions.shielding.config import ShieldingParams
from jaff.errors import ParserError
from jaff.physics.photo_reactions.config import RadiationConfig
from jaff.physics.photo_reactions.shielding import ShieldingFunction, ShieldingResult
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

            def get_shielding(self, reaction):
                return ShieldingResult(Float(1.0))

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


class TestPhotoReactionShielding:
    def test_factor_is_on_the_photo_reaction(self, leiden_network):
        photo = leiden_network.reactions[CO_KEY, "photo"].photo

        from sympy.core.function import AppliedUndef

        calls = photo.shielding_factor.atoms(AppliedUndef)
        assert len(calls) == 2  # one interpolation per shielding species
        assert all(band.k.has(*calls) for band in photo.bands)

    def test_unshielded_reaction_has_no_factor(self, make_network):
        net = make_network(
            CO_PHOTO, name="co.jet", funcfile=False, radiation={"bands": [6.0, 13.6]}
        )

        assert net.reactions[CO_KEY, "photo"].photo.shielding_factor is None

    def test_plugin_returns_its_factor_and_tables(self, leiden_network):
        rxn = leiden_network.reactions[CO_KEY, "photo"]

        result = Leiden().get_shielding(rxn)

        assert isinstance(result, ShieldingResult)
        assert result.factor == rxn.photo.shielding_factor
        assert set(result.tables) == {f"shielding_{CO_KEY}.hdf5"}

    def test_shielding_without_radiation_raises(self, make_network):
        config = {
            "network": {
                "reactions": {"reaction": {CO_KEY: {"photo": {"shielding": LEIDEN}}}}
            }
        }

        with pytest.raises(ParserError, match="radiation"):
            make_network(CO_PHOTO, name="co.jet", funcfile=False, config=config)

    def test_photochemistry_has_no_shielding(self):
        from jaff.physics import Photochemistry

        assert not hasattr(Photochemistry, "shielding")


class TestRegistry:
    def test_public_names(self):
        from jaff.physics.photo_reactions import shielding

        assert shielding.shielding_function("leiden", CO_KEY).name == "leiden"
        assert ("leiden", None) in shielding.shielding_functions()
        assert callable(shielding.register)

    def test_private_names_are_gone(self):
        from jaff.physics.photo_reactions import shielding

        for name in ("_register", "_get_shielding_function", "_discover_shielding"):
            assert not hasattr(shielding, name)


class TestTypedParams:
    def test_params_are_typed_on_load(self):
        from jaff.physics.photo_reactions.shielding.config import ShieldingConfig
        from jaff.physics.photo_reactions.shielding.global_.leiden import LeidenParams

        config = ShieldingConfig.model_validate(
            {"type": "leiden", "shielded_by": ["self"]}
        )

        assert isinstance(config.params, LeidenParams)
        assert config.params.radiation == "isrf"

    def test_dump_is_the_table_as_written(self):
        from jaff.physics.photo_reactions.shielding.config import ShieldingConfig

        table = {"type": "leiden", "shielded_by": ["self"]}

        assert ShieldingConfig.model_validate(table).model_dump() == table

    def test_unknown_type_is_rejected_on_load(self):
        from jaff.physics.photo_reactions.shielding.config import ShieldingConfig

        with pytest.raises(ValidationError, match="nope"):
            ShieldingConfig.model_validate({"type": "nope"})

    def test_plugin_reads_the_stored_settings(self, leiden_network):
        rxn = leiden_network.reactions[CO_KEY, "photo"]

        assert Leiden().params(rxn) is rxn.config.photo.shielding.params

    def test_one_name_has_one_settings_model(self):
        from jaff.errors import RegistrationError
        from jaff.physics.photo_reactions.shielding import register

        class Other(ShieldingFunction):
            name = "leiden"
            reaction = "CO._PHOTON__C.O"

            def get_shielding(self, reaction):
                return ShieldingResult(Float(1.0))

        with pytest.raises(RegistrationError, match="leiden"):
            register(Other)
