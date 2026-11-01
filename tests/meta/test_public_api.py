# ABOUTME: Every jaff package's __all__ is a list of unique strings naming its
# ABOUTME: attributes, as the generated API reference (jaff._utils.api_docs) requires

import importlib
import pkgutil

import pytest

import jaff


def _packages():
    names = ["jaff"]
    for info in pkgutil.walk_packages(jaff.__path__, "jaff."):
        if info.ispkg and ".templates" not in info.name:
            names.append(info.name)

    return names


@pytest.mark.parametrize("name", _packages())
def test_all_names_its_attributes(name):
    module = importlib.import_module(name)
    exported = getattr(module, "__all__", None)
    if exported is None:
        pytest.skip(f"{name} has no __all__")

    assert all(isinstance(item, str) for item in exported), name
    assert len(exported) == len(set(exported)), f"{name}: duplicate entries"
    assert [item for item in exported if not hasattr(module, item)] == []
