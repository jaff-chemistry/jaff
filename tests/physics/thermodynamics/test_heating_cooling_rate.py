# ABOUTME: Global heating/cooling-rate plugins: [network.thermodynamics] heating_cooling_rate,
# ABOUTME: a list summed with the .jfunc heatingCoolingRate into dEdt_extra; exact copies rejected

import pytest
from pydantic import ValidationError
from sympy import Float

from tests.helpers.heating_cooling_rate_plugins import install
from jaff import Network
from jaff.errors import ParserError
from jaff.physics.thermodynamics.config import ThermodynamicsConfig

LINES = "@format:idx, R, R, P, rate\n1, H, D, He, n_H\n"


@pytest.fixture(autouse=True)
def test_plugins(monkeypatch):
    """Only the test plugins are registered, for each test."""
    install(monkeypatch)


def _network(tmp_path, entries=None, jfunc=None):
    path = tmp_path / "net.dat"
    path.write_text(LINES, encoding="utf-8")
    if jfunc is not None:
        (tmp_path / "net.jfunc").write_text(jfunc, encoding="utf-8")
    thermodynamics = {"heating_cooling_rate": entries} if entries is not None else {}

    return Network(str(path), thermodynamics=thermodynamics)


def _extra(net):
    return net.thermodynamics.dEdt_extra.volumetric


def test_without_plugins_or_jfunc_the_extra_rate_is_zero(tmp_path):
    assert float(_extra(_network(tmp_path))) == pytest.approx(0.0)


def test_a_configured_plugin_sets_the_extra_rate(tmp_path):
    net = _network(tmp_path, [{"type": "constant_rate_test", "value": 2e-25}])

    assert float(_extra(net)) == pytest.approx(2e-25)


def test_plugins_are_summed(tmp_path):
    net = _network(
        tmp_path,
        [{"type": "constant_rate_test", "value": 2e-25}, {"type": "cooling_test"}],
    )
    tgas = net.symbols.tgas

    assert float(_extra(net).subs(tgas, 100.0)) == pytest.approx(2e-25 - 1e-28)


def test_the_jfunc_rate_is_added_to_the_plugins(tmp_path):
    net = _network(
        tmp_path,
        [{"type": "constant_rate_test", "value": 2e-25}],
        jfunc="@function heatingCoolingRate()\n    @return 3e-25\n",
    )

    assert float(_extra(net)) == pytest.approx(5e-25)


def test_same_type_with_other_settings_counts_each(tmp_path):
    net = _network(
        tmp_path,
        [
            {"type": "constant_rate_test", "value": 2e-25},
            {"type": "constant_rate_test", "value": 1e-25},
        ],
    )

    assert float(_extra(net)) == pytest.approx(3e-25)


def test_an_exact_copy_is_rejected():
    entry = {"type": "constant_rate_test", "value": 2e-25}

    with pytest.raises(ValidationError, match="heating_cooling_rate.1"):
        ThermodynamicsConfig.model_validate(
            {"heating_cooling_rate": [entry, dict(entry)]}
        )


def test_the_rate_enters_the_total_and_temperature_rates(tmp_path):
    net = _network(tmp_path, [{"type": "constant_rate_test", "value": 2e-25}])

    assert net.thermodynamics.dEdt_tot().volumetric.has(Float(2e-25))
    assert net.thermodynamics.dTdt_tot().has(Float(2e-25))


def test_unknown_type_is_rejected(tmp_path):
    with pytest.raises((ParserError, ValidationError), match="missing_model"):
        _network(tmp_path, [{"type": "missing_model"}])


def test_a_missing_required_setting_is_rejected(tmp_path):
    with pytest.raises((ParserError, ValidationError), match="value"):
        _network(tmp_path, [{"type": "constant_rate_test"}])


def test_type_is_required():
    with pytest.raises(ValidationError, match="type is required"):
        ThermodynamicsConfig.model_validate({"heating_cooling_rate": [{"value": 1.0}]})


def test_dumps_as_flat_tables():
    entries = [{"type": "constant_rate_test", "value": 2e-25}, {"type": "cooling_test"}]
    config = ThermodynamicsConfig.model_validate({"heating_cooling_rate": entries})

    assert config.model_dump()["heating_cooling_rate"] == entries
