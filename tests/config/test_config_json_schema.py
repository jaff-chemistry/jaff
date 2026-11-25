# ABOUTME: JSON Schemas of jaff.toml / jaffgen.toml (for editor linting): real TOML text
# ABOUTME: validates against them, and the committed schema files match the models

import json
import subprocess
import sys
import tomllib
from pathlib import Path

import pytest
from jsonschema import Draft202012Validator

from jaff.config.base import attribute_descriptions
from jaff.config.json_schema import SCHEMA_DIR, schema_files


@pytest.fixture(scope="module")
def schemas():
    return {name: Draft202012Validator(schema) for name, schema in schema_files().items()}


def _errors(validator, text):
    return list(validator.iter_errors(tomllib.loads(text)))


VALID_JAFF_TOML = [
    "",
    '[network]\nlabel = "gow"\nstrict = true\nfuncfile = "aux.jfunc"\n',
    "[network]\nfuncfile = false\n",
    '[network.reactions]\nt_cutoff = "extrapolate"\nduplicate_policy = "error"\n',
    '[network.reactions.photo]\npi_database = "verner"\nuse_proxy_photoreaction = false\n',
    '[network.reactions.photo.shielding]\nshielded_by = ["self", "H2"]\nradiation = "draine"\n',
    '[network.reactions.reaction."H2._PHOTON__H.H".photo.shielding]\n'
    'type = "hg2015"\nmin_ncol = 1e-35\n',
    '[network.reactions.reaction."CO._PHOTON__C.O"]\nt_cutoff = "clip"\n'
    '[network.reactions.reaction."CO._PHOTON__C.O".photo]\n'
    'shielding = { type = "leiden", shielded_by = ["self"] }\n',
    '[network.radiation]\nbands = [6, 13.6, "inf"]\nprofile_index = [0, 1.5]\nc = "c_hat"\n',
    "[network.dust]\n",
    "[network.dust]\nrv = 4.0\n[network.dust.pe]\nthreshold_low = 5.0\n",
    "[network.thermodynamics.eos]\n",
    "[network.thermodynamics.eos]\ngamma = 1.4\n",
    '[network.thermodynamics.eos]\ntype = "multi_gamma"\ndefault_gamma = 1.4\n'
    "gamma_map = { H2 = 1.4 }\n",
    '[network.thermodynamics.eos]\ntype = "fermi_degenerate"\n',
]

INVALID_JAFF_TOML = {
    "old-rsl": '[network.radiation]\nbands = [6.0, 13.6]\nrsl = "c_hat"\n',
    "old-T_cutoff": '[network.reactions]\nT_cutoff = "clip"\n',
    "old-rates": '[network.rates]\nt_cutoff = "clip"\n',
    "old-eos-place": '[network.eos]\ntype = "ideal"\n',
    "bad-choice": '[network.reactions]\nt_cutoff = "clamp"\n',
    "string-bool": '[network]\nstrict = "yes"\n',
    "bad-band": '[network.radiation]\nbands = [6.0, "big"]\n',
    "missing-bands": '[network.radiation]\nmode = "u"\n',
    "gamma-too-small": "[network.thermodynamics.eos]\ngamma = 1.0\n",
    "eos-unknown-type": '[network.thermodynamics.eos]\ntype = "polytrope"\n',
    "eos-foreign-key": '[network.thermodynamics.eos]\ntype = "fermi_degenerate"\ngamma = 1.4\n',
    "shielding-typo": '[network.reactions.reaction."H2._PHOTON__H.H".photo.shielding]\n'
    'type = "hg2015"\nmin_ncl = 1e-35\n',
    "shielding-unknown-type": '[network.reactions.photo.shielding]\ntype = "nope"\n',
    "leiden-missing-shielded_by": '[network.reactions.photo.shielding]\ntype = "leiden"\n',
    "jaffgen-in-jaff-toml": '[jaffgen]\ntemplate = "microphysics"\n',
}


@pytest.mark.parametrize("text", VALID_JAFF_TOML)
def test_valid_jaff_toml(schemas, text):
    assert _errors(schemas["jaff"], text) == []


@pytest.mark.parametrize("text", INVALID_JAFF_TOML.values(), ids=INVALID_JAFF_TOML.keys())
def test_invalid_jaff_toml(schemas, text):
    assert _errors(schemas["jaff"], text) != []


def test_jaffgen_toml_holds_jaffgen_and_network(schemas):
    text = (
        '[jaffgen]\ntemplate = "microphysics"\noutput_dir = "out"\necho = false\n'
        'input_files = ["a.cpp"]\n'
        '[[jaffgen.table]]\n[jaffgen.table.source]\npath = "default"\n'
        '[jaffgen.table.target."/co".attrs]\ntmax = "/co.max"\n'
        '[network.radiation]\nbands = [13.6, "inf"]\n'
    )

    assert _errors(schemas["jaffgen"], text) == []


def test_jaffgen_toml_rejects_top_level_table(schemas):
    assert _errors(schemas["jaffgen"], '[[table]]\n[table.source]\npath = "x"\n') != []


def test_speed_of_light_default_is_in_schema():
    radiation = schema_files()["jaff"]["$defs"]["RadiationConfig"]

    assert radiation["properties"]["c"]["default"] == pytest.approx(2.99792458e10)


@pytest.mark.parametrize("name", ["jaff", "jaffgen"])
def test_committed_schema_matches_models(name):
    committed = json.loads(
        (SCHEMA_DIR / f"{name}.schema.json").read_text(encoding="utf-8")
    )

    assert committed == schema_files()[name], (
        "Run `python -m jaff.config.json_schema` to regenerate the schemas"
    )


@pytest.fixture(scope="module")
def defs():
    return schema_files()["jaff"]["$defs"]


class TestDescriptions:
    """Keys carry their numpydoc Attributes entry as description (hover/completion)."""

    def test_key_description_from_attributes(self, defs):
        description = defs["RadiationConfig"]["properties"]["bands"]["description"]

        assert "band edges" in description

    def test_table_description_is_the_summary(self, defs):
        description = defs["RadiationConfig"]["description"]

        assert "Attributes" not in description
        assert description.startswith("Radiation field")

    def test_reference_to_a_table_is_described(self, defs):
        assert (
            "Radiation field"
            in defs["NetworkConfig"]["properties"]["radiation"]["description"]
        )

    def test_rst_markup_becomes_markdown(self, defs):
        description = defs["ReactionsConfig"]["properties"]["reaction"]["description"]

        assert "``" not in description
        assert ":meth:" not in description

    def test_shielding_plugin_keys_are_described(self, defs):
        branches = defs["ShieldingConfig"]["anyOf"]
        leiden = next(b for b in branches if b["properties"]["type"]["const"] == "leiden")

        assert "Shielding species" in leiden["properties"]["shielded_by"]["description"]

    def test_eos_still_requires_type(self, defs):
        assert "type" in defs["MultiGammaEosConfig"]["required"]


class TestAttributeDescriptions:
    def test_parses_multi_line_entries(self):
        doc = """Summary line.

        More text.

        Attributes
        ----------
        alpha : float
            First line
            continues here.
        beta : str or None
            Uses ``code`` and :class:`~jaff.Thing`.

        Notes
        -----
        Not an attribute.
        """

        summary, fields = attribute_descriptions(doc)

        assert summary == "Summary line.\n\nMore text."
        assert fields == {
            "alpha": "First line continues here.",
            "beta": "Uses `code` and `Thing`.",
        }

    def test_without_attributes_section(self):
        assert attribute_descriptions("Just a summary.") == ("Just a summary.", {})


class TestNoNull:
    """TOML has no null: optional keys are simply absent, never typed null."""

    def test_no_null_anywhere(self):
        for name, schema in schema_files().items():
            assert '"null"' not in json.dumps(schema), name

    def test_optional_table_is_a_plain_reference(self, defs):
        assert defs["NetworkConfig"]["properties"]["radiation"]["$ref"] == (
            "#/$defs/RadiationConfig"
        )

    def test_optional_value_keeps_its_type(self, defs):
        assert defs["NetworkConfig"]["properties"]["label"]["type"] == "string"


class TestSchemaPath:
    @pytest.mark.parametrize("name", ["jaff", "jaffgen"])
    def test_prints_installed_schema_path(self, name):
        result = subprocess.run(
            [sys.executable, "-m", "jaff.config.json_schema", "--path", name],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0, result.stderr
        assert Path(result.stdout.strip()) == SCHEMA_DIR / f"{name}.schema.json"

    def test_unknown_name_fails(self):
        result = subprocess.run(
            [sys.executable, "-m", "jaff.config.json_schema", "--path", "nope"],
            capture_output=True,
            text=True,
        )

        assert result.returncode != 0
