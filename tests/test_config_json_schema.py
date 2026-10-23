# ABOUTME: JSON Schemas of jaff.toml / jaffgen.toml (for editor linting): real TOML text
# ABOUTME: validates against them, and the committed schema files match the models

import json
import tomllib

import pytest
from jsonschema import Draft202012Validator

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
    committed = json.loads((SCHEMA_DIR / f"{name}.schema.json").read_text())

    assert committed == schema_files()[name], (
        "Run `python -m jaff.config.json_schema` to regenerate the schemas"
    )
