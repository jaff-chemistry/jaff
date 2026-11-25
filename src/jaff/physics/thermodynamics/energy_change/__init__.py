# ABOUTME: Registry of the per-reaction energy-change plugins, keyed by (name, reaction);
# ABOUTME: they live in plugins/ (none in jaff's core), found on the first lookup
"""Per-reaction energy-change plugins (``thermodynamics.energy_change`` tables)."""

from __future__ import annotations

from jaff.common.registry import Registry
from jaff.config.plugin_table import check_shared_settings
from jaff.errors import ParserError

from ._base import EnergyChangeFunction
from .config import EnergyChangeConfig, EnergyChangeParams, ReactionThermodynamicsConfig

_REGISTRY: Registry[EnergyChangeFunction] = Registry(
    "energy change function", packages=[]
)
"""The energy-change plugins, keyed by ``(name, reaction)``; all in ``plugins/``."""


def register(cls: type[EnergyChangeFunction]) -> type[EnergyChangeFunction]:
    """Class decorator registering an energy-change plugin under ``(name, reaction)``.

    Parameters
    ----------
    cls : type[EnergyChangeFunction]
        The plugin class; one instance is registered.

    Returns
    -------
    type[EnergyChangeFunction]
        *cls*, unchanged.

    Raises
    ------
    RegistrationError
        If a plugin with the same ``(name, reaction)`` is already registered, or one
        with the same name has another ``config_model()``.
    """
    check_shared_settings(_REGISTRY, cls)

    return _REGISTRY.register(cls)


def energy_change_functions() -> dict[tuple[str, str | None], EnergyChangeFunction]:
    """Every registered energy-change plugin, keyed by ``(name, reaction)`` (a copy)."""
    return _REGISTRY.items()


def energy_change_function(name: str, reaction: str | None) -> EnergyChangeFunction:
    """The energy-change plugin of type *name* for the serialized *reaction*.

    A plugin bound to *reaction* wins over one for any reaction of the same name.

    Raises
    ------
    ParserError
        If no plugin of that type applies to *reaction*.
    """
    plugin = _REGISTRY.find(name, reaction)
    if plugin is None:
        raise ParserError(f"Unknown energy change type {name!r} for reaction {reaction}")

    return plugin


__all__ = [
    "EnergyChangeConfig",
    "EnergyChangeFunction",
    "EnergyChangeParams",
    "ReactionThermodynamicsConfig",
    "energy_change_function",
    "energy_change_functions",
    "register",
]
