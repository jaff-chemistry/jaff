# ABOUTME: Config model of the thermodynamics ([network.thermodynamics]): the
# ABOUTME: equation of state table and the global heating/cooling-rate plugins
"""Config of the thermodynamics (``[network.thermodynamics]``)."""

from __future__ import annotations


from pydantic import field_validator

from ...config.base import ConfigModel
from .eos.config import EosConfig, IdealEosConfig
from .heating_cooling_rate.config import (
    HeatingCoolingRateConfig,
    check_heating_cooling_rates,
)

__all__ = [
    "ThermodynamicsConfig",
]


class ThermodynamicsConfig(ConfigModel):
    """Thermodynamics (``[network.thermodynamics]``).

    Attributes
    ----------
    eos : EosConfig
        Equation of state; its ``type`` selects the model and its keys.
    heating_cooling_rate : list[HeatingCoolingRateConfig]
        Global heating/cooling mechanisms, each ``{ type = "<plugin>", ... }`` (erg
        cm⁻³ s⁻¹; negative cools).  Their rates are summed and added to the ``.jfunc``
        ``heatingCoolingRate``; an entry repeating another exactly (same type and
        settings) is rejected.  Empty by default.
    """

    eos: EosConfig = IdealEosConfig()
    heating_cooling_rate: list[HeatingCoolingRateConfig] = []

    @field_validator("heating_cooling_rate")
    @classmethod
    def _no_copies(
        cls, entries: list[HeatingCoolingRateConfig]
    ) -> list[HeatingCoolingRateConfig]:
        """Reject an entry repeating another exactly (it would be counted twice)."""
        check_heating_cooling_rates(entries)

        return entries
