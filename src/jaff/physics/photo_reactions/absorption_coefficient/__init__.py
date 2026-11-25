# ABOUTME: Registry of the absorption coefficient plugins, keyed by name; they live in
# ABOUTME: plugins/ (none in jaff's core), found on the first lookup
"""Absorption coefficient plugins (``[network.radiation] absorption_coefficient``)."""

from __future__ import annotations

from jaff.common.registry import Registry
from jaff.config.plugin_table import check_shared_settings
from jaff.errors import ParserError, RegistrationError

from ._base import AbsorptionCoefficientFunction
from .config import AbsorptionCoefficientConfig, AbsorptionCoefficientParams

_REGISTRY: Registry[AbsorptionCoefficientFunction] = Registry(
    "absorption coefficient function", packages=[]
)
"""The absorption coefficient plugins, keyed by ``(name, None)``; all in ``plugins/``."""


def register(
    cls: type[AbsorptionCoefficientFunction],
) -> type[AbsorptionCoefficientFunction]:
    """Class decorator registering an absorption coefficient plugin under its name.

    Parameters
    ----------
    cls : type[AbsorptionCoefficientFunction]
        The plugin class; one instance is registered.

    Returns
    -------
    type[AbsorptionCoefficientFunction]
        *cls*, unchanged.

    Raises
    ------
    RegistrationError
        If a plugin with the same name is already registered, or *cls* is bound to a
        reaction.
    """
    if getattr(cls, "reaction", None) is not None:
        raise RegistrationError(
            f"Absorption coefficient plugin {cls.name} cannot be bound to a reaction"
        )
    check_shared_settings(_REGISTRY, cls)

    return _REGISTRY.register(cls)


def absorption_coefficient_functions() -> dict[
    tuple[str, str | None], AbsorptionCoefficientFunction
]:
    """Every registered absorption coefficient plugin, by ``(name, None)`` (a copy)."""
    return _REGISTRY.items()


def absorption_coefficient_function(name: str) -> AbsorptionCoefficientFunction:
    """The absorption coefficient plugin of type *name*.

    Raises
    ------
    ParserError
        If no plugin has that name.
    """
    plugin = _REGISTRY.find(name, None)
    if plugin is None:
        raise ParserError(f"Unknown absorption coefficient type {name!r}")

    return plugin


__all__ = [
    "AbsorptionCoefficientConfig",
    "AbsorptionCoefficientFunction",
    "AbsorptionCoefficientParams",
    "absorption_coefficient_function",
    "absorption_coefficient_functions",
    "register",
]
