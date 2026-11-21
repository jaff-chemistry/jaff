# ABOUTME: Config of the emissivity plugins: one table (type + settings) per entry of
# ABOUTME: [network.radiation] emissivity, summed with the .jfunc globalEmissivity
"""Config of the emissivity plugins (``[network.radiation] emissivity`` entries)."""

from __future__ import annotations

from typing import ClassVar

from ....config.base import LowerStr
from ....config.plugin_table import PluginParams, PluginTable

__all__ = [
    "EmissivityConfig",
    "EmissivityParams",
]


class EmissivityParams(PluginParams):
    """Base of the settings of one emissivity type (every key but ``type``).

    Each plugin exposes its subclass through
    :meth:`~jaff.physics.photo_reactions.emissivity.EmissivityFunction.config_model`;
    this base has no keys, so a plugin without settings accepts none.
    """


class EmissivityConfig(PluginTable):
    """One emissivity source (an entry of ``[network.radiation] emissivity``).

    Written flat in TOML (``{ type = "<plugin>", <its settings>... }``): ``type`` names
    the emissivity plugin and every other key is one of its settings.  Its spectrum
    (erg cm⁻³ s⁻¹ eV⁻¹, in ``ph_E``) is added to the ``.jfunc`` ``globalEmissivity``
    and to the other entries, then integrated over each radiation group.

    Attributes
    ----------
    type : str
        Emissivity model (case-insensitive); required.
    params : EmissivityParams
        The plugin's settings, validated by its ``config_model()``.
    """

    registry: ClassVar = staticmethod(lambda: _registry())
    label: ClassVar[str] = "emissivity"
    title: ClassVar[str] = "Emissivity"
    summary: ClassVar[str] = "A volume emission spectrum in ph_E (erg cm⁻³ s⁻¹ eV⁻¹)"

    type: LowerStr
    params: EmissivityParams = EmissivityParams()


def _registry():
    """The emissivity plugins' registry (imported here: it imports this module)."""
    from . import _REGISTRY

    return _REGISTRY
