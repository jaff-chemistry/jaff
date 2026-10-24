# ABOUTME: Registry of the line-shielding plugins, keyed by (name, reaction);
# ABOUTME: discovers the plugin subpackages on first lookup
from __future__ import annotations

from jaff.common._helper import import_subpackages
from jaff.errors import ParserError, RegistrationError

from ._base import ShieldingFunction, ShieldingResult
from .config import ShieldingParams

_REGISTER: dict[tuple[str, str | None], ShieldingFunction] = {}
_DISCOVERED: bool = False


def register(cls: type[ShieldingFunction]) -> type[ShieldingFunction]:
    """Class decorator registering a shielding plugin under ``(name, reaction)``.

    Parameters
    ----------
    cls : type[ShieldingFunction]
        The plugin class; one instance is registered.

    Returns
    -------
    type[ShieldingFunction]
        *cls*, unchanged.

    Raises
    ------
    RegistrationError
        If a plugin with the same ``(name, reaction)`` is already registered, or
        one with the same name has another ``config_model()``.
    """
    key = (cls.name, cls.reaction)
    if key in _REGISTER:
        raise RegistrationError(
            f"Key already exists for shielding function: {cls.name}, {cls.reaction}"
        )

    # A shielding type's settings model is chosen by its name alone (the
    # reaction is not known when a config table is read).
    for (name, _), plugin in _REGISTER.items():
        if name == cls.name and plugin.config_model() is not cls.config_model():
            raise RegistrationError(
                f"Shielding functions named {cls.name} must share one config_model()"
            )

    _REGISTER[key] = cls()

    return cls


def _discover() -> None:
    """Import the plugin subpackages once, so their ``register`` calls run."""
    global _DISCOVERED
    if not _DISCOVERED:
        import_subpackages(__name__)
        _DISCOVERED = True


def shielding_functions() -> dict[tuple[str, str | None], ShieldingFunction]:
    """Every registered shielding plugin, keyed by ``(name, reaction)``.

    ``reaction`` is ``None`` for a global plugin (any reaction).

    Returns
    -------
    dict[tuple[str, str or None], ShieldingFunction]
        A copy of the registry.
    """
    _discover()

    return dict(_REGISTER)


def settings_model(name: str) -> type[ShieldingParams]:
    """The settings model (``config_model()``) of the shielding type *name*.

    Every plugin of one name shares it (see :func:`register`).

    Parameters
    ----------
    name : str
        The ``shielding.type`` value (case-insensitive).

    Returns
    -------
    type[ShieldingParams]

    Raises
    ------
    ParserError
        If no plugin has that name.
    """
    _discover()
    for (plugin_name, _), plugin in _REGISTER.items():
        if plugin_name == name.lower():
            return plugin.config_model()

    raise ParserError(f"Unknown shielding type {name!r}")


def shielding_function(name: str, reaction: str | None) -> ShieldingFunction:
    """The shielding plugin of type *name* for the serialized *reaction*.

    A plugin bound to *reaction* wins over a global one of the same name.

    Parameters
    ----------
    name : str
        The ``shielding.type`` value (case-insensitive).
    reaction : str or None
        The serialized reaction, or ``None`` for a global plugin only.

    Returns
    -------
    ShieldingFunction

    Raises
    ------
    ParserError
        If no plugin of that type applies to *reaction*.
    """
    _discover()
    key = (name.lower(), reaction)
    if key not in _REGISTER:
        # Fall back to a global function registered with reaction=None.
        key = (name.lower(), None)

    if key not in _REGISTER:
        raise ParserError(
            f"Invalid shielding function type {name} or reaction {reaction} is not supported"
        )

    return _REGISTER[key]


__all__ = [
    "register",
    "settings_model",
    "shielding_function",
    "shielding_functions",
    "ShieldingFunction",
    "ShieldingResult",
]
