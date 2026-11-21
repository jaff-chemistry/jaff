# ABOUTME: Radiation plugins: [network.radiation] emissivity / absorption_coefficient
# ABOUTME: lists, summed with the .jfunc global functions, then integrated / averaged

import pytest
from pydantic import ValidationError

from tests.radiation_plugins import install
from jaff.errors import ParserError
from jaff.physics.photo_reactions.config import RadiationConfig

BANDS = [6.0, 13.6, 100.0]
NET = "H + H -> H2  []  1e-17\n"
U_MODE = {"mode": "u", "profile_index": 0.5}
FLAT_PHOTONS = {"profile_index": 2}


@pytest.fixture(autouse=True)
def test_plugins(monkeypatch):
    """Only the test plugins are registered, for each test."""
    install(monkeypatch)


@pytest.fixture
def radiation(make_network, tmp_path):
    """The radiation field of a network with *jfunc* and the radiation table *table*."""

    def _make(table, jfunc=""):
        (tmp_path / "r.jfunc").write_text(jfunc)
        return make_network(NET, name="r.jet", radiation=table).radiation

    return _make


def _emission(value):
    return {"type": "flat_emission_test", "value": value}


def _absorption(slope):
    return {"type": "linear_absorption_test", "slope": slope}


def test_without_entries_or_jfunc_both_are_none(radiation):
    rad = radiation({"bands": BANDS})

    assert rad.emissivity is None
    assert rad.absorption_coefficient is None


def test_an_emissivity_plugin_is_integrated_per_group(radiation):
    rad = radiation({"bands": BANDS, "emissivity": [_emission(1e-20)], **U_MODE})

    assert float(rad.emissivity) == pytest.approx(1e-20)
    assert [float(g.emissivity) for g in rad.groups] == pytest.approx(
        [1e-20 * (13.6 - 6.0), 1e-20 * (100.0 - 13.6)]
    )


def test_emissivity_plugins_and_the_jfunc_are_summed(radiation):
    jfunc = "@function globalEmissivity()\n@return 3e-20\n"
    entries = [_emission(1e-20), _emission(2e-20)]
    rad = radiation({"bands": BANDS, "emissivity": entries, **U_MODE}, jfunc)

    assert float(rad.emissivity) == pytest.approx(6e-20)
    assert float(rad.total_emissivity) == pytest.approx(6e-20 * (100.0 - 6.0))


def test_an_absorption_plugin_is_averaged_per_group(radiation):
    table = {"bands": BANDS, "absorption_coefficient": [_absorption(1e-3)]}
    rad = radiation({**table, **FLAT_PHOTONS})

    assert [float(g.absorption_coefficient) for g in rad.groups] == pytest.approx(
        [1e-3 * (6.0 + 13.6) / 2, 1e-3 * (13.6 + 100.0) / 2]
    )


def test_absorption_plugins_and_the_jfunc_are_summed(radiation):
    jfunc = "@function globalAbsorptionCoefficient()\n@return 2e-3\n"
    table = {"bands": BANDS, "absorption_coefficient": [_absorption(1e-3)]}
    rad = radiation({**table, **FLAT_PHOTONS}, jfunc)

    assert float(rad.absorption_coefficient.subs(rad.ph_E, 10.0)) == pytest.approx(
        1e-2 + 2e-3
    )


@pytest.mark.parametrize(
    ("key", "entry"),
    [("emissivity", _emission(1e-20)), ("absorption_coefficient", _absorption(1e-3))],
)
def test_an_exact_copy_is_rejected(key, entry):
    with pytest.raises(ValidationError, match=rf"{key}\.1: repeats entry 0"):
        RadiationConfig.model_validate({"bands": BANDS, key: [entry, dict(entry)]})


@pytest.mark.parametrize("key", ["emissivity", "absorption_coefficient"])
def test_unknown_type_is_rejected(key):
    with pytest.raises((ParserError, ValidationError), match="missing_model"):
        RadiationConfig.model_validate({"bands": BANDS, key: [{"type": "missing_model"}]})


def test_a_missing_required_setting_is_rejected():
    entry = {"type": "flat_emission_test"}

    with pytest.raises(ValidationError, match="value"):
        RadiationConfig.model_validate({"bands": BANDS, "emissivity": [entry]})


def test_dumps_as_flat_tables():
    config = RadiationConfig.model_validate(
        {"bands": BANDS, "emissivity": [_emission(1e-20)]}
    )

    assert config.model_dump()["emissivity"] == [_emission(1e-20)]
