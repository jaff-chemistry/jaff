# ABOUTME: Registry of the global heating/cooling-rate plugins, keyed by name; they live in
# ABOUTME: plugins/ (none in jaff's core), found on the first lookup
"""Global heating/cooling-rate plugins (``[network.thermodynamics] heating_cooling_rate``)."""

from __future__ import annotations

from jaff.common.registry import Registry
from jaff.config.plugin_table import check_shared_settings
from jaff.errors import ParserError, RegistrationError

from ._base import HeatingCoolingRateFunction
from .config import HeatingCoolingRateConfig, HeatingCoolingRateParams

_REGISTRY: Registry[HeatingCoolingRateFunction] = Registry(
    "heating/cooling rate function", packages=[]
)
"""The heating/cooling-rate plugins, keyed by ``(name, None)``; all in ``plugins/``."""


def register(
    cls: type[HeatingCoolingRateFunction],
) -> type[HeatingCoolingRateFunction]:
    """Class decorator registering a heating/cooling-rate plugin under its name.

    Parameters
    ----------
    cls : type[HeatingCoolingRateFunction]
        The plugin class; one instance is registered.

    Returns
    -------
    type[HeatingCoolingRateFunction]
        *cls*, unchanged.

    Raises
    ------
    RegistrationError
        If a plugin with the same name is already registered, or *cls* is bound to a
        reaction.
    """
    if getattr(cls, "reaction", None) is not None:
        raise RegistrationError(
            f"Heating/cooling-rate plugin {cls.name} cannot be bound to a reaction"
        )
    check_shared_settings(_REGISTRY, cls)

    return _REGISTRY.register(cls)


def heating_cooling_rate_functions() -> dict[
    tuple[str, str | None], HeatingCoolingRateFunction
]:
    """Every registered heating/cooling-rate plugin, keyed by ``(name, None)`` (a copy)."""
    return _REGISTRY.items()


def heating_cooling_rate_function(name: str) -> HeatingCoolingRateFunction:
    """The heating/cooling-rate plugin of type *name*.

    Raises
    ------
    ParserError
        If no plugin has that name.
    """
    plugin = _REGISTRY.find(name, None)
    if plugin is None:
        raise ParserError(f"Unknown heating/cooling rate type {name!r}")

    return plugin


__all__ = [
    "HeatingCoolingRateConfig",
    "HeatingCoolingRateFunction",
    "HeatingCoolingRateParams",
    "heating_cooling_rate_function",
    "heating_cooling_rate_functions",
    "register",
]
