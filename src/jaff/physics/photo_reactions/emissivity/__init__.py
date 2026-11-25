# ABOUTME: Registry of the emissivity plugins, keyed by name; they live in plugins/
# ABOUTME: (none in jaff's core), found on the first lookup
"""Emissivity plugins (``[network.radiation] emissivity``)."""

from __future__ import annotations

from jaff.common.registry import Registry
from jaff.config.plugin_table import check_shared_settings
from jaff.errors import ParserError, RegistrationError

from ._base import EmissivityFunction
from .config import EmissivityConfig, EmissivityParams

_REGISTRY: Registry[EmissivityFunction] = Registry("emissivity function", packages=[])
"""The emissivity plugins, keyed by ``(name, None)``; all in ``plugins/``."""


def register(cls: type[EmissivityFunction]) -> type[EmissivityFunction]:
    """Class decorator registering an emissivity plugin under its name.

    Parameters
    ----------
    cls : type[EmissivityFunction]
        The plugin class; one instance is registered.

    Returns
    -------
    type[EmissivityFunction]
        *cls*, unchanged.

    Raises
    ------
    RegistrationError
        If a plugin with the same name is already registered, or *cls* is bound to a
        reaction.
    """
    if getattr(cls, "reaction", None) is not None:
        raise RegistrationError(
            f"Emissivity plugin {cls.name} cannot be bound to a reaction"
        )
    check_shared_settings(_REGISTRY, cls)

    return _REGISTRY.register(cls)


def emissivity_functions() -> dict[tuple[str, str | None], EmissivityFunction]:
    """Every registered emissivity plugin, keyed by ``(name, None)`` (a copy)."""
    return _REGISTRY.items()


def emissivity_function(name: str) -> EmissivityFunction:
    """The emissivity plugin of type *name*.

    Raises
    ------
    ParserError
        If no plugin has that name.
    """
    plugin = _REGISTRY.find(name, None)
    if plugin is None:
        raise ParserError(f"Unknown emissivity type {name!r}")

    return plugin


__all__ = [
    "EmissivityConfig",
    "EmissivityFunction",
    "EmissivityParams",
    "emissivity_function",
    "emissivity_functions",
    "register",
]
