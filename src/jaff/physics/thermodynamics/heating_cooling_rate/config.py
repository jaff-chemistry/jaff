# ABOUTME: Config of the global heating/cooling-rate plugins: one heating_cooling_rate table
# ABOUTME: (type + settings) per entry of [network.thermodynamics] heating_cooling_rate
"""Config of the global heating/cooling rates (``heating_cooling_rate`` tables)."""

from __future__ import annotations

from typing import ClassVar

from ....config.base import LowerStr
from ....config.plugin_table import PluginParams, PluginTable

__all__ = [
    "HeatingCoolingRateConfig",
    "HeatingCoolingRateParams",
    "check_heating_cooling_rates",
]


class HeatingCoolingRateParams(PluginParams):
    """Base of the settings of one heating/cooling-rate type (every key but ``type``).

    Each plugin exposes its subclass through
    :meth:`~jaff.physics.thermodynamics.heating_cooling_rate.HeatingCoolingRateFunction.config_model`;
    this base has no keys, so a plugin without settings accepts none.
    """


class HeatingCoolingRateConfig(PluginTable):
    """One global heating/cooling mechanism (an entry of ``heating_cooling_rate``).

    Written flat in TOML (``{ type = "<plugin>", <its settings>... }``): ``type`` names
    the heating/cooling-rate plugin and every other key is one of its settings.  Its
    rate (erg cm⁻³ s⁻¹; negative cools) is added to the ``.jfunc``
    ``heatingCoolingRate`` and to the other entries.

    Attributes
    ----------
    type : str
        Heating/cooling-rate model (case-insensitive); required.
    params : HeatingCoolingRateParams
        The plugin's settings, validated by its ``config_model()``.
    """

    registry: ClassVar = staticmethod(lambda: _registry())
    label: ClassVar[str] = "heating/cooling rate"
    title: ClassVar[str] = "Heating/cooling rate"
    summary: ClassVar[str] = (
        "A global heating/cooling mechanism (erg cm⁻³ s⁻¹; negative cools)"
    )

    type: LowerStr
    params: HeatingCoolingRateParams = HeatingCoolingRateParams()


def _registry():
    """The heating/cooling-rate plugins' registry (imported here: it imports this module)."""
    from . import _REGISTRY

    return _REGISTRY


def check_heating_cooling_rates(entries: list[HeatingCoolingRateConfig]) -> None:
    """Reject an entry that repeats an earlier one exactly (same type and settings).

    Every entry is summed, so an exact copy would count one mechanism twice; the same
    type with other settings (e.g. another species) is a separate term.

    Raises
    ------
    ValueError
        As ``"heating_cooling_rate.<index>: <message>"`` for the first copy.
    """
    seen: list[dict] = []
    for index, entry in enumerate(entries):
        dumped = entry.model_dump()
        if dumped in seen:
            raise ValueError(
                f"heating_cooling_rate.{index}: repeats entry {seen.index(dumped)} "
                f"({entry.type!r} with the same settings); each entry is summed"
            )
        seen.append(dumped)
