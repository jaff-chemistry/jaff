# ABOUTME: Per-reaction energy_change plugins: selected by [...reaction."<s>".thermodynamics]
# ABOUTME: energy_change; the .jfunc deltaE<N> wins, then the plugin, the file's dE, then 0

import pytest
from pydantic import ValidationError
from sympy import Float

from tests.energy_change_plugins import install
from jaff import Network
from jaff.core.reaction.config import ReactionsConfig
from jaff.errors import ParserError

KEY = "D.H__He"


@pytest.fixture(autouse=True)
def test_plugins(monkeypatch):
    """Only the test plugins are registered, for each test."""
    install(monkeypatch)


LINES = "@format:idx, R, R, P, rate\n1, H, D, He, n_H\n"


def _network(tmp_path, energy_change=None, jfunc=None, reactions=None):
    path = tmp_path / "net.dat"
    path.write_text(LINES, encoding="utf-8")
    if jfunc is not None:
        (tmp_path / "net.jfunc").write_text(jfunc, encoding="utf-8")
    if reactions is None:
        own = (
            {"thermodynamics": {"energy_change": energy_change}} if energy_change else {}
        )
        reactions = {"reaction": {KEY: own}} if own else {}

    return Network(str(path), reactions=reactions)


def _reaction(net):
    (reaction,) = list(net.reactions)
    return reaction


def test_configured_plugin_sets_the_energy_change(tmp_path):
    net = _network(tmp_path, {"type": "constant_test", "value": 2e-12})

    assert _reaction(net).energy.per_event == Float(2e-12)


def test_the_energy_change_enters_the_chemical_energy_rate(tmp_path):
    net = _network(tmp_path, {"type": "constant_test", "value": 2e-12})

    assert net.thermodynamics.dEdt_chemical().volumetric.has(Float(2e-12))


def test_a_plugin_bound_to_the_reaction_applies(tmp_path):
    assert _reaction(_network(tmp_path, {"type": "bound_test"})).energy.per_event == Float(-1e-12)


def test_the_jfunc_deltaE_wins_over_the_plugin(tmp_path):
    net = _network(
        tmp_path,
        {"type": "constant_test", "value": 2e-12},
        jfunc="@function deltaE0()\n    @return 3e-12\n",
    )

    assert float(_reaction(net).energy.per_event) == pytest.approx(3e-12)


def test_without_a_plugin_or_jfunc_the_energy_change_is_zero(tmp_path):
    assert float(_reaction(_network(tmp_path)).energy.per_event) == pytest.approx(0.0)


def test_unknown_type_is_rejected(tmp_path):
    with pytest.raises((ParserError, ValidationError), match="missing_model"):
        _network(tmp_path, {"type": "missing_model"})


def test_a_missing_required_setting_is_rejected(tmp_path):
    with pytest.raises((ParserError, ValidationError), match="value"):
        _network(tmp_path, {"type": "constant_test"})


def test_a_plugin_bound_to_another_reaction_is_rejected(tmp_path):
    other = {
        "reaction": {
            "H.H__H2": {"thermodynamics": {"energy_change": {"type": "bound_test"}}}
        }
    }

    with pytest.raises(ValidationError, match="bound_test"):
        ReactionsConfig.model_validate(other)


def test_energy_change_is_per_reaction_only():
    with pytest.raises(ValidationError, match="thermodynamics"):
        ReactionsConfig.model_validate(
            {"thermodynamics": {"energy_change": {"type": "constant_test", "value": 1.0}}}
        )


def test_type_is_required():
    entry = {"reaction": {KEY: {"thermodynamics": {"energy_change": {"value": 1.0}}}}}

    with pytest.raises(ValidationError, match="type is required"):
        ReactionsConfig.model_validate(entry)


def test_model_resolver_returns_none_without_a_plugin(tmp_path):
    assert _reaction(_network(tmp_path)).model("energy_change") is None


def test_model_resolver_returns_the_plugin_result(tmp_path):
    net = _network(tmp_path, {"type": "constant_test", "value": 5e-13})

    assert _reaction(net).model("energy_change") == Float(5e-13)


def test_model_resolver_rejects_an_unknown_kind(tmp_path):
    with pytest.raises(ValueError, match="kind"):
        _reaction(_network(tmp_path)).model("nonsense")
