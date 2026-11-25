# ABOUTME: Tests for config.plugin_table: a flat {type = ..., settings...} table selecting a
# ABOUTME: registered plugin, its settings validated by the plugin's model; dump and schema

from typing import ClassVar

import pytest
from pydantic import ValidationError

from jaff.common.registry import Registry
from jaff.config.plugin_table import PluginParams, PluginTable, settings_model
from jaff.errors import RegistrationError

MODELS: Registry = Registry("test model", packages=[])


class ScaleParams(PluginParams):
    """Settings of the scale model.

    Attributes
    ----------
    factor : float
        Scale factor.
    """

    factor: float = 1.0


@MODELS.register
class Scale:
    name, reaction = "scale", None

    @staticmethod
    def config_model() -> type[PluginParams]:
        return ScaleParams


@MODELS.register
class Plain:
    name, reaction = "plain", "H2._PHOTON__H.H"

    @staticmethod
    def config_model() -> type[PluginParams]:
        return PluginParams


class ModelTable(PluginTable):
    """A test table (``type`` required)."""

    registry: ClassVar = staticmethod(lambda: MODELS)
    title: ClassVar[str] = "Model"


class DefaultTable(PluginTable):
    """A test table defaulting to ``scale``."""

    registry: ClassVar = staticmethod(lambda: MODELS)
    title: ClassVar[str] = "Model"
    default_type: ClassVar[str | None] = "scale"


def test_flat_table_selects_the_plugin_and_validates_its_settings():
    table = ModelTable.model_validate({"type": "Scale", "factor": 2.5})

    assert table.type == "scale"
    assert isinstance(table.params, ScaleParams)
    assert table.params.factor == pytest.approx(2.5)


def test_unknown_setting_is_rejected():
    with pytest.raises(ValidationError, match="nope"):
        ModelTable.model_validate({"type": "scale", "nope": 1})


def test_unknown_type_is_rejected():
    with pytest.raises(ValidationError, match="missing"):
        ModelTable.model_validate({"type": "missing"})


def test_type_is_required_without_a_default():
    with pytest.raises(ValidationError, match="type"):
        ModelTable.model_validate({"factor": 2.0})


def test_default_type_applies_when_omitted():
    assert DefaultTable.model_validate({"factor": 3.0}).type == "scale"


def test_dumps_flat():
    table = ModelTable.model_validate({"type": "scale", "factor": 2.5})

    assert table.model_dump() == {"type": "scale", "factor": 2.5}


def test_schema_has_one_alternative_per_plugin():
    schema = ModelTable.model_json_schema()
    names = [branch["properties"]["type"]["const"] for branch in schema["anyOf"]]

    assert names == ["plain", "scale"]
    assert "factor" in schema["anyOf"][1]["properties"]
    assert schema["anyOf"][1]["required"] == ["type"]
    assert schema["title"] == "Model"


def test_settings_model_by_name():
    assert settings_model(MODELS, "SCALE") is ScaleParams


def test_plugins_of_one_name_must_share_a_settings_model():
    registry = Registry("test model", packages=[])
    registry.register(
        type(
            "A",
            (),
            {
                "name": "x",
                "reaction": None,
                "config_model": staticmethod(lambda: ScaleParams),
            },
        )
    )
    other = type(
        "B",
        (),
        {
            "name": "x",
            "reaction": "H2._PHOTON__H.H",
            "config_model": staticmethod(lambda: PluginParams),
        },
    )

    from jaff.config.plugin_table import check_shared_settings

    with pytest.raises(RegistrationError, match="config_model"):
        check_shared_settings(registry, other)


def test_schema_without_plugins_accepts_any_typed_table():
    empty = Registry("test model", packages=[])

    class EmptyTable(PluginTable):
        registry: ClassVar = staticmethod(lambda: empty)
        title: ClassVar[str] = "Empty"

    (branch,) = EmptyTable.model_json_schema()["anyOf"]

    assert branch["properties"]["type"] == {"type": "string"}
    assert branch["required"] == ["type"]
    assert branch["additionalProperties"] is True


def test_an_exact_copy_in_a_list_is_rejected():
    from jaff.config.plugin_table import check_no_copies

    entries = [ModelTable.model_validate({"type": "scale", "factor": 2.0})] * 2

    with pytest.raises(ValueError, match=r"models\.1: repeats entry 0"):
        check_no_copies(entries, "models")


def test_same_type_with_other_settings_is_not_a_copy():
    from jaff.config.plugin_table import check_no_copies

    entries = [
        ModelTable.model_validate({"type": "scale", "factor": 2.0}),
        ModelTable.model_validate({"type": "scale", "factor": 3.0}),
    ]

    check_no_copies(entries, "models")
