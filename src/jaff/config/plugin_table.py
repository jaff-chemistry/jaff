# ABOUTME: PluginTable: a flat config table {type = "<plugin>", <its settings>...} selecting a
# ABOUTME: registered plugin, its settings validated by the plugin's model (PluginParams)
"""Config tables that select a plugin by ``type`` (e.g. ``shielding``, ``energy_change``)."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable, ClassVar

from pydantic import GetJsonSchemaHandler, model_serializer, model_validator
from pydantic.json_schema import JsonSchemaValue
from pydantic_core import CoreSchema

from ..errors import ParserError, RegistrationError
from .base import ConfigModel, LowerStr

if TYPE_CHECKING:
    from ..common.registry import Registry

__all__ = ["PluginParams", "PluginTable", "check_shared_settings", "settings_model"]


class PluginParams(ConfigModel):
    """Base of one plugin type's settings (every key of its table but ``type``).

    A plugin returns its subclass from ``config_model()``; this base has no keys, so a
    plugin without settings accepts none.
    """


def settings_model(
    registry: Registry[Any], name: str, label: str | None = None
) -> type[PluginParams]:
    """The settings model (``config_model()``) of the plugins named *name*.

    Every plugin of one name shares it (see :func:`check_shared_settings`).

    Parameters
    ----------
    registry : Registry
        The plugins' registry.
    name : str
        The ``type`` value (case-insensitive).
    label : str, optional
        What the plugins are, for the message; the registry's kind by default.

    Returns
    -------
    type[PluginParams]

    Raises
    ------
    ParserError
        If no plugin has that name.
    """
    for (plugin_name, _), plugin in registry.items().items():
        if plugin_name == name.lower():
            return plugin.config_model()

    raise ParserError(f"Unknown {label or registry.kind} type {name!r}")


def check_shared_settings(registry: Registry[Any], cls: Any) -> None:
    """Reject *cls* if a registered plugin of its name has another ``config_model()``.

    A table's settings model is chosen by ``type`` alone (the reaction is not known
    when the table is read), so every plugin of one name must share it.  Uses the
    plugins registered so far, without importing the plugin packages (called while
    they register).

    Raises
    ------
    RegistrationError
        If one of the same name has another ``config_model()``.
    """
    for (name, _), plugin in registry.registered().items():
        if name == cls.name.lower() and plugin.config_model() is not cls.config_model():
            raise RegistrationError(
                f"{registry.kind.capitalize()}s named {cls.name} must share one "
                "config_model()"
            )


class PluginTable(ConfigModel):
    """A config table selecting a plugin: ``{ type = "<name>", <settings>... }``.

    Written flat in TOML: ``type`` names the plugin and every other key is one of its
    settings, validated by the plugin's ``config_model()`` (a :class:`PluginParams`
    subclass) when the table is read.  A subclass sets :attr:`registry` (a callable
    returning the plugins' registry), :attr:`label`, :attr:`title`, :attr:`summary` and
    :attr:`default_type`; without a default, ``type`` is required.  The table is
    atomic: a more specific one replaces a more general one whole.

    Attributes
    ----------
    type : str
        The plugin (case-insensitive).
    params : PluginParams
        The plugin's settings, validated by its ``config_model()``.
    """

    _atomic: ClassVar[bool] = True

    registry: ClassVar[Callable[[], Registry[Any]]]
    """Returns the plugins' registry (a callable: the registry imports the config)."""

    label: ClassVar[str] = "plugin"
    """What the plugins are, in messages and schema titles (e.g. ``"shielding"``)."""

    title: ClassVar[str] = "Plugin"
    """The editor schema's title."""

    summary: ClassVar[str] = "A plugin"
    """The start of the editor schema's description."""

    default_type: ClassVar[str | None] = None
    """The plugin used when ``type`` is omitted; ``None`` makes ``type`` required."""

    type: LowerStr
    params: PluginParams = PluginParams()

    @model_validator(mode="before")
    @classmethod
    def _split_type(cls, data: Any) -> Any:
        """Split a flat table into ``type`` and the plugin's validated ``params``."""
        if not isinstance(data, dict):
            return data

        if "type" not in data and cls.default_type is None:
            raise ValueError(f"{cls.label}: type is required")

        name = str(data.get("type", cls.default_type)).lower()
        settings = {k: v for k, v in data.items() if k != "type"}
        try:
            model = settings_model(cls.registry(), name, cls.label)
        except ParserError as e:
            raise ValueError(str(e)) from e

        # A ValidationError raised here keeps its key paths (``<table>.<key>``).
        return {"type": name, "params": model.model_validate(settings)}

    @model_serializer
    def _flat(self) -> dict[str, Any]:
        """Dump as the flat TOML table it is written as."""
        return {"type": self.type, **self.params.model_dump(exclude_unset=True)}

    @classmethod
    def __get_pydantic_json_schema__(
        cls, core_schema: CoreSchema, handler: GetJsonSchemaHandler
    ) -> JsonSchemaValue:
        """Describe the flat TOML table: one alternative per registered plugin.

        Each alternative fixes ``type`` to the plugin's name (required, but for the
        default type) and holds the plugin's ``config_model()`` keys.  Whether a
        reaction-bound plugin fits the reaction is checked on load.  Without any
        registered plugin, the one alternative is any table with a ``type``.
        """
        default = cls.default_type
        plugins = sorted(
            cls.registry().items().items(),
            key=lambda item: (item[0][0], item[0][1] or ""),
        )
        branches: list[JsonSchemaValue] = []
        for (name, reaction), plugin in plugins:
            params = plugin.config_model().model_json_schema()
            if "$defs" in params:
                raise NotImplementedError(
                    f"{cls.label} {name}: nested models in its settings are not supported"
                )

            only = f" (only for {reaction})" if reaction is not None else ""
            branches.append(
                {
                    "type": "object",
                    "title": f"{name} {cls.label}{only}",
                    "properties": {
                        "type": {"const": name},
                        **params.get("properties", {}),
                    },
                    "required": [
                        *params.get("required", []),
                        *(["type"] if name != default else []),
                    ],
                    "additionalProperties": False,
                }
            )

        if not branches:
            # No plugin registered: any typed table; jaff checks it on load.
            branches.append(
                {
                    "type": "object",
                    "title": f"{cls.label} (no built-in plugins)",
                    "properties": {"type": {"type": "string"}},
                    "required": ["type"],
                    "additionalProperties": True,
                }
            )

        choice = f" (default {default!r})" if default is not None else ""

        return {
            "title": cls.title,
            "description": f"{cls.summary}; `type` selects the model{choice} and its keys.",
            "anyOf": branches,
        }
