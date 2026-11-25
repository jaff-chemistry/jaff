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
