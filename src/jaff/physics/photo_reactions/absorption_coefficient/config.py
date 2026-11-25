# ABOUTME: Config of the absorption coefficient plugins: one table (type + settings) per
# ABOUTME: [network.radiation] absorption_coefficient entry, summed with the .jfunc one
"""Config of the absorption coefficient plugins (``absorption_coefficient`` entries)."""

from __future__ import annotations

from typing import ClassVar

from ....config.base import LowerStr
from ....config.plugin_table import PluginParams, PluginTable

__all__ = [
    "AbsorptionCoefficientConfig",
    "AbsorptionCoefficientParams",
]


class AbsorptionCoefficientParams(PluginParams):
    """Base of the settings of one absorption coefficient type (every key but ``type``).

    Each plugin exposes its subclass through
    ``AbsorptionCoefficientFunction.config_model()``
    (:mod:`jaff.physics.photo_reactions.absorption_coefficient`); this base has no keys, so a plugin without settings accepts none.
    """


class AbsorptionCoefficientConfig(PluginTable):
    """One absorber (an entry of ``[network.radiation] absorption_coefficient``).

    Written flat in TOML (``{ type = "<plugin>", <its settings>... }``): ``type`` names
    the absorption coefficient plugin and every other key is one of its settings.  Its
    spectrum (cm⁻¹, in ``ph_E``) is added to the ``.jfunc``
    ``globalAbsorptionCoefficient`` and to the other entries, then averaged over each
    radiation group.

    Attributes
    ----------
    type : str
        Absorption coefficient model (case-insensitive); required.
    params : AbsorptionCoefficientParams
        The plugin's settings, validated by its ``config_model()``.
    """

    registry: ClassVar = staticmethod(lambda: _registry())
    label: ClassVar[str] = "absorption coefficient"
    title: ClassVar[str] = "Absorption coefficient"
    summary: ClassVar[str] = "An absorption coefficient spectrum in ph_E (cm⁻¹)"

    type: LowerStr
    params: AbsorptionCoefficientParams = AbsorptionCoefficientParams()


def _registry():
    """The absorption coefficient plugins' registry (imported here: it imports us)."""
    from . import _REGISTRY

    return _REGISTRY
