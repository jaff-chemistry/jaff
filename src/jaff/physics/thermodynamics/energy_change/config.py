# ABOUTME: Config of the per-reaction energy-change plugins: the energy_change table (type +
# ABOUTME: settings) in a reaction's thermodynamics table, and its settings base
"""Config of the per-reaction energy change (``energy_change`` tables)."""

from __future__ import annotations

from typing import ClassVar

from ....config.base import ConfigModel, LowerStr, dotted
from ....config.plugin_table import PluginParams, PluginTable
from ....errors import ParserError

__all__ = [
    "EnergyChangeConfig",
    "EnergyChangeParams",
    "ReactionThermodynamicsConfig",
    "check_energy_change",
]


class EnergyChangeParams(PluginParams):
    """Base of the settings of one energy-change type (every key but ``type``).

    Each plugin exposes its subclass through
    :meth:`~jaff.physics.thermodynamics.energy_change.EnergyChangeFunction.config_model`;
    this base has no keys, so a plugin without settings accepts none.
    """


class EnergyChangeConfig(PluginTable):
    """Energy change of one reaction (an ``energy_change`` table).

    Written flat in TOML (``{ type = "<plugin>", <its settings>... }``): ``type`` names
    the energy-change plugin and every other key is one of its settings.  The energy
    one event of the reaction adds to the gas (erg; negative cools); a ``.jfunc``
    ``deltaE<N>`` takes preference over it.

    Attributes
    ----------
    type : str
        Energy-change model (case-insensitive); required.
    params : EnergyChangeParams
        The plugin's settings, validated by its ``config_model()``.
    """

    registry: ClassVar = staticmethod(lambda: _registry())
    label: ClassVar[str] = "energy change"
    title: ClassVar[str] = "Energy change"
    summary: ClassVar[str] = (
        "Energy one event of the reaction adds to the gas (erg; negative cools)"
    )

    type: LowerStr
    params: EnergyChangeParams = EnergyChangeParams()


class ReactionThermodynamicsConfig(ConfigModel):
    """Thermodynamic settings of one reaction (its ``thermodynamics`` table).

    Attributes
    ----------
    energy_change : EnergyChangeConfig or None
        The energy-change plugin of the reaction (``{ type = "<plugin>", ... }``): the
        energy one event adds to the gas (erg; negative cools).  A ``.jfunc``
        ``deltaE<N>`` takes preference over it.
    """

    energy_change: EnergyChangeConfig | None = None


def _registry():
    """The energy-change plugins' registry (imported here: it imports this module)."""
    from . import _REGISTRY

    return _REGISTRY


def check_energy_change(table: EnergyChangeConfig | None, reaction: str) -> None:
    """Check that a plugin of *table*'s type applies to the serialized *reaction*.

    Its settings were validated when the table was read.

    Raises
    ------
    ValueError
        If no plugin of that type applies to *reaction*, as ``"<key path>: <message>"``
        with the path relative to ``[network.reactions]``.
    """
    if table is None:
        return

    from . import energy_change_function

    where = ("reaction", reaction, "thermodynamics", "energy_change")
    try:
        energy_change_function(table.type, reaction)
    except ParserError as e:
        raise ValueError(f"{dotted(where)}: {e}") from e
