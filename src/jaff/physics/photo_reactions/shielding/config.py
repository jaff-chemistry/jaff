"""Config of the photo-reaction line shielding (``shielding`` tables).

A ``shielding`` table names a shielding plugin by ``type``; its other keys are
that plugin's settings, checked by the plugin's ``config_model()`` (a
:class:`ShieldingParams` subclass) once the reaction it applies to is known.
"""

from __future__ import annotations

from typing import Any, ClassVar, Dict

from pydantic import ValidationError, model_validator

from ....config.base import ConfigModel, LowerStr
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
        Shielding plugin name, lower-cased.
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


def check_shielding(shielding: ShieldingConfig | None, reaction: str | None) -> None:
    """Check *shielding* against its plugin, for *reaction* (``None``: global).

    Raises
    ------
    ValueError
        If no plugin of that type applies to *reaction*, or the plugin rejects
        a setting.
    """
    if shielding is None:
        return

    # Imported here: the plugin registry imports this module.
    from . import _get_shielding_function

    where = (
        "reactions.photo.shielding"
        if reaction is None
        else f'reactions.reaction."{reaction}".photo.shielding'
    )
    try:
        plugin = _get_shielding_function(shielding.type, reaction)
        plugin.config_model().model_validate(shielding.params)
    except (ParserError, ValidationError) as e:
        raise ValueError(f"{where}: {e}") from e
