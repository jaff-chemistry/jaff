"""Config of the photo-reaction line shielding (``shielding`` tables).

A ``shielding`` table names a shielding plugin by ``type``; its other keys are
that plugin's settings, checked by the plugin's ``config_model()`` (a
:class:`ShieldingParams` subclass) once the reaction it applies to is known.
"""

from __future__ import annotations

from typing import Any, ClassVar, Dict, List

from pydantic import (
    GetJsonSchemaHandler,
    ValidationError,
    model_serializer,
    model_validator,
)
from pydantic.json_schema import JsonSchemaValue
from pydantic_core import CoreSchema

from ....config.base import (
    ConfigModel,
    dotted,
    format_validation_errors,
    KeyPath,
    LowerStr,
)
from ....errors import ParserError


class ShieldingParams(ConfigModel):
    """Base of the settings of one shielding type (every key but ``type``).

    Each shielding plugin exposes its subclass through
    :meth:`~jaff.physics.photo_reactions.shielding.ShieldingFunction.config_model`;
    this base has no keys, so a plugin without settings accepts none.
    """


class ShieldingConfig(ConfigModel):
    """Shielding of a photo reaction (a ``shielding`` table).

    Written flat in TOML (``{ type = "leiden", shielded_by = ["self"] }``):
    ``type`` names the shielding plugin and every other key is one of its
    settings, checked by the plugin's ``config_model()`` once the reaction it
    applies to is known (see :class:`~jaff.core.reaction.config.ReactionsConfig`).

    The table is atomic: a per-reaction one replaces the global one whole.

    Attributes
    ----------
    type : str
        Shielding model (case-insensitive).
    params : dict
        The plugin's settings, as written.
    """

    _atomic: ClassVar[bool] = True

    type: LowerStr = "leiden"
    params: Dict[str, Any] = {}

    @model_validator(mode="before")
    @classmethod
    def _split_type(cls, data: Any) -> Any:
        """Split a flat table into ``type`` and the plugin ``params``."""
        if not isinstance(data, dict):
            return data

        params = {k: v for k, v in data.items() if k != "type"}

        return {"type": data.get("type", "leiden"), "params": params}

    @model_serializer
    def _flat(self) -> Dict[str, Any]:
        """Dump as the flat TOML table it is written as."""
        return {"type": self.type, **self.params}

    @classmethod
    def __get_pydantic_json_schema__(
        cls, core_schema: CoreSchema, handler: GetJsonSchemaHandler
    ) -> JsonSchemaValue:
        """Describe the flat TOML table: one alternative per shielding plugin.

        Each alternative fixes ``type`` to the plugin's name (required, but for
        the default type) and holds the plugin's ``config_model()`` keys.
        Whether a reaction-bound plugin fits the reaction is checked on load.
        """
        # Imported here: the plugin registry imports this module.
        from . import _REGISTER, _discover_shielding

        _discover_shielding()
        default = cls.model_fields["type"].default
        branches: List[JsonSchemaValue] = []
        for (name, reaction), plugin in sorted(
            _REGISTER.items(), key=lambda item: (item[0][0], item[0][1] or "")
        ):
            params = plugin.config_model().model_json_schema()
            if "$defs" in params:
                raise NotImplementedError(
                    f"shielding {name}: nested models in its settings are not supported"
                )

            only = f" (only for {reaction})" if reaction is not None else ""
            branches.append(
                {
                    "type": "object",
                    "title": f"{name} shielding{only}",
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

        return {
            "title": "Shielding",
            "description": "Line shielding of a photo reaction; `type` selects the "
            f"model (default {default!r}) and its keys.",
            "anyOf": branches,
        }


def check_shielding(shielding: ShieldingConfig | None, reaction: str | None) -> None:
    """Check *shielding* against its plugin, for *reaction* (``None``: global).

    Raises
    ------
    ValueError
        If no plugin of that type applies to *reaction*, or the plugin rejects
        a setting; one ``"<key path>: <message>"`` line per problem, the path
        relative to ``[network.reactions]``.
    """
    if shielding is None:
        return

    # Imported here: the plugin registry imports this module.
    from . import _get_shielding_function

    # Key path of the table within [network.reactions].
    where: KeyPath = (
        ("photo", "shielding")
        if reaction is None
        else ("reaction", reaction, "photo", "shielding")
    )
    try:
        plugin = _get_shielding_function(shielding.type, reaction)
        plugin.config_model().model_validate(shielding.params)
    except ParserError as e:
        raise ValueError(f"{dotted(where)}: {e}") from e
    except ValidationError as e:
        raise ValueError("\n".join(format_validation_errors(e, where))) from e
