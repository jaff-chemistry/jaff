# ABOUTME: Tests for common.plugins.load_plugins: the repo-root plugins folder imported as
# ABOUTME: jaff.plugins, every public subpackage at any depth, once; missing folder is fine

import sys
from pathlib import Path

import pytest

import jaff.common.plugins as plugins_module
from jaff.common.plugins import load_plugins
from jaff.common.registry import Registry


def _package(path: Path, body: str = "") -> None:
    """A package at *path*, with an ``__init__.py`` in every folder up to ``plugins``."""
    path.mkdir(parents=True, exist_ok=True)
    folder = path
    while folder.name != "plugins" and not (folder / "__init__.py").exists():
        (folder / "__init__.py").write_text("", encoding="utf-8")
        folder = folder.parent
    (path / "__init__.py").write_text(body, encoding="utf-8")


def _plugin_modules() -> list[str]:
    return [
        m for m in sys.modules if m == "jaff.plugins" or m.startswith("jaff.plugins.")
    ]


@pytest.fixture
def fresh(monkeypatch):
    """No plugins loaded yet; the plugin modules the test imports are dropped after it."""
    monkeypatch.setattr(plugins_module, "_LOADED", False)
    for name in _plugin_modules():
        monkeypatch.delitem(sys.modules, name)
    yield
    for name in _plugin_modules():
        sys.modules.pop(name)


@pytest.fixture
def tree(tmp_path, monkeypatch, fresh):
    """An empty plugins folder in place of the repo's."""
    root = tmp_path / "plugins"
    _package(root)
    monkeypatch.setattr(plugins_module, "PLUGINS_DIR", root)

    return root


def test_nested_packages_are_imported_as_jaff_plugins(tree):
    _package(tree / "reactions" / "H2__PHOTON__H_H" / "shielding", "LOADED = True\n")
    load_plugins()

    module = sys.modules["jaff.plugins.reactions.H2__PHOTON__H_H.shielding"]
    assert module.LOADED


def test_private_folders_are_not_imported(tree):
    _package(tree / "reactions" / "_utils", "raise RuntimeError('imported')\n")

    load_plugins()

    assert "jaff.plugins.reactions._utils" not in sys.modules


def test_relative_imports_work_inside_a_plugin(tree):
    _package(tree / "physics" / "_helpers", "VALUE = 3\n")
    _package(
        tree / "physics" / "model", "from .._helpers import VALUE\nDOUBLE = 2 * VALUE\n"
    )

    load_plugins()

    assert sys.modules["jaff.plugins.physics.model"].DOUBLE == 6


def test_plugins_load_once(tree, monkeypatch):
    monkeypatch.setattr(plugins_module, "IMPORTS", [], raising=False)
    _package(
        tree / "physics", "from jaff.common import plugins\nplugins.IMPORTS.append(1)\n"
    )

    load_plugins()
    load_plugins()

    assert plugins_module.IMPORTS == [1]


def test_a_missing_folder_loads_nothing(tmp_path, monkeypatch, fresh):
    monkeypatch.setattr(plugins_module, "PLUGINS_DIR", tmp_path / "absent")

    load_plugins()

    assert "jaff.plugins" not in sys.modules


def test_a_plugin_registers_into_a_registry_on_first_lookup(tree, monkeypatch):
    models = Registry("model", packages=[])
    monkeypatch.setattr(plugins_module, "MODELS", models, raising=False)
    _package(
        tree / "physics" / "simple",
        "from jaff.common import plugins\n"
        "@plugins.MODELS.register\n"
        "class Simple:\n"
        "    name, reaction = 'simple', None\n",
    )

    assert models.find("simple", None) is not None


def test_the_repo_plugins_folder_is_used():
    from jaff.config import PLUGINS_DIR

    assert PLUGINS_DIR.name == "plugins"
    assert (PLUGINS_DIR / "__init__.py").is_file()
