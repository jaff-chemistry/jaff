"""Config of the thermodynamics (``[network.thermodynamics]``)."""

from __future__ import annotations


from ...config.base import ConfigModel
from .eos.config import EosConfig, IdealEosConfig

__all__ = [
    "ThermodynamicsConfig",
]


class ThermodynamicsConfig(ConfigModel):
    """Thermodynamics (``[network.thermodynamics]``).

    Attributes
    ----------
    eos : EosConfig
        Equation of state; its ``type`` selects the model and its keys.
    """

    eos: EosConfig = IdealEosConfig()
