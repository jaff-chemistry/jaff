# ABOUTME: Tests for common.registry.Registry: plugins keyed by (name, reaction), the
# ABOUTME: reaction's own beating a global one, duplicate rejection and one-time discovery

import pytest

import jaff.common.registry as registry_module
from jaff.common.registry import Registry
from jaff.errors import RegistrationError


def _plugin(name: str, reaction: str | None = None) -> type:
    return type("Plugin", (), {"name": name, "reaction": reaction})


@pytest.fixture
def registry() -> Registry:
    return Registry("test plugin", packages=[])


def test_registered_plugin_is_found_by_name(registry):
    registry.register(_plugin("model"))

    assert registry.find("model", None).name == "model"


def test_names_are_case_insensitive(registry):
    registry.register(_plugin("Model"))

    assert registry.find("MODEL", None) is not None


def test_a_reaction_plugin_beats_a_global_one(registry):
    registry.register(_plugin("model"))
    registry.register(_plugin("model", "H2._PHOTON__H.H"))

    assert registry.find("model", "H2._PHOTON__H.H").reaction == "H2._PHOTON__H.H"


def test_a_global_plugin_serves_any_reaction(registry):
    registry.register(_plugin("model"))

    assert registry.find("model", "CO._PHOTON__C.O").reaction is None


def test_a_reaction_plugin_serves_only_its_reaction(registry):
    registry.register(_plugin("model", "H2._PHOTON__H.H"))

    assert registry.find("model", "CO._PHOTON__C.O") is None
    assert registry.find("model", None) is None


def test_unknown_name_is_none(registry):
    assert registry.find("missing", None) is None


def test_register_returns_the_class(registry):
    cls = _plugin("model")

    assert registry.register(cls) is cls


def test_duplicate_key_is_rejected(registry):
    registry.register(_plugin("model", "H2._PHOTON__H.H"))

    with pytest.raises(RegistrationError, match="test plugin: MODEL, H2._PHOTON__H.H"):
        registry.register(_plugin("MODEL", "H2._PHOTON__H.H"))


def test_items_is_a_copy(registry):
    registry.register(_plugin("model"))
    items = registry.items()
    items.clear()

    assert registry.find("model", None) is not None
    assert list(registry.items()) == [("model", None)]


def test_packages_are_discovered_once_on_first_use(monkeypatch):
    imported = []
    monkeypatch.setattr(registry_module, "import_subpackages", imported.append)
    registry = Registry("test plugin", packages=["pkg.a", "pkg.b"])

    assert imported == []
    registry.find("model", None)
    registry.items()

    assert imported == ["pkg.a", "pkg.b"]


def test_registered_does_not_discover(monkeypatch):
    imported = []
    monkeypatch.setattr(registry_module, "import_subpackages", imported.append)
    registry = Registry("test plugin", packages=["pkg.a"])
    registry.register(_plugin("model"))

    assert list(registry.registered()) == [("model", None)]
    assert imported == []
