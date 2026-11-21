# ABOUTME: Config models of the photo-reaction line shielding: the flat
# ABOUTME: shielding table (type + plugin settings) and the plugins' settings base
"""Config of the photo-reaction line shielding (``shielding`` tables).

A ``shielding`` table names a shielding plugin by ``type``; its other keys are
that plugin's settings, validated by the type's ``config_model()`` (a
:class:`ShieldingParams` subclass) when the table is read.  Whether the type
applies to the reaction is checked once the reaction is known
(:func:`check_shielding`).
"""

from __future__ import annotations

from typing import ClassVar

from ....config.base import KeyPath, LowerStr, dotted
from ....config.plugin_table import PluginParams, PluginTable
from ....errors import ParserError

__all__ = [
    "ShieldingConfig",
    "ShieldingParams",
    "check_shielding",
]


class ShieldingParams(PluginParams):
    """Base of the settings of one shielding type (every key but ``type``).

    Each shielding plugin exposes its subclass through
    :meth:`~jaff.physics.photo_reactions.shielding.ShieldingFunction.config_model`;
    this base has no keys, so a plugin without settings accepts none.
    """


class ShieldingConfig(PluginTable):
    """Shielding of a photo reaction (a ``shielding`` table).

    Written flat in TOML (``{ type = "leiden", shielded_by = ["self"] }``):
    ``type`` names the shielding plugin and every other key is one of its
    settings, validated by the type's ``config_model()``; whether the type
    applies to the reaction is checked once it is known (see
    :class:`~jaff.core.reaction.config.ReactionsConfig`).

    The table is atomic: a per-reaction one replaces the global one whole.

    Attributes
    ----------
    type : str
        Shielding model (case-insensitive): ``"leiden"`` (the default; any
        photo reaction with a Leiden table), or ``"db1996"`` / ``"hg2015"``
        (only for ``"H2._PHOTON__H.H"``).
    params : ShieldingParams
        The plugin's settings, validated by its ``config_model()``.
    """

    registry: ClassVar = staticmethod(lambda: _registry())
    label: ClassVar[str] = "shielding"
    title: ClassVar[str] = "Shielding"
    summary: ClassVar[str] = "Line shielding of a photo reaction"
    default_type: ClassVar[str | None] = "leiden"

    type: LowerStr = "leiden"
    params: ShieldingParams = ShieldingParams()


def _registry():
    """The shielding plugins' registry (imported here: it imports this module)."""
    from . import _REGISTRY

    return _REGISTRY


def check_shielding(shielding: ShieldingConfig | None, reaction: str | None) -> None:
    """Check that a plugin of *shielding*'s type applies to *reaction* (``None``: global).

    Its settings were validated when the table was read.

    Raises
    ------
    ValueError
        If no plugin of that type applies to *reaction*, as
        ``"<key path>: <message>"`` with the path relative to
        ``[network.reactions]``.
    """
    if shielding is None:
        return

    # Imported here: the plugin registry imports this module.
    from . import shielding_function

    # Key path of the table within [network.reactions].
    where: KeyPath = (
        ("photo", "shielding")
        if reaction is None
        else ("reaction", reaction, "photo", "shielding")
    )
    try:
        shielding_function(shielding.type, reaction)
    except ParserError as e:
        raise ValueError(f"{dotted(where)}: {e}") from e
