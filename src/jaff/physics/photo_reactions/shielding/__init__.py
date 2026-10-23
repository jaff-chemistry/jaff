# ABOUTME: Registry of the line-shielding plugins, keyed by (name, reaction);
# ABOUTME: discovers the plugin subpackages on first lookup
from jaff.common._helper import import_subpackages
from jaff.errors import ParserError, RegistrationError

from ._base import ShieldingFunction

_REGISTER: dict[tuple[str, str | None], ShieldingFunction] = {}
_DISCOVERED: bool = False


def _register(cls: type[ShieldingFunction]) -> type[ShieldingFunction]:
    key = (cls.name, cls.reaction)
    if key in _REGISTER:
        raise RegistrationError(
            f"Key already exists for shielding function: {cls.name}, {cls.reaction}"
        )

    _REGISTER[key] = cls()
    return cls


def _discover_shielding() -> None:
    global _DISCOVERED
    if not _DISCOVERED:
        import_subpackages(__name__)
        _DISCOVERED = True


def _get_shielding_function(name: str, reaction: str | None):
    _discover_shielding()
    name = name.lower()
    key = (name, reaction)

    if key not in _REGISTER:
        # Fall back to a global function registered with reaction=None.
        key = (name, None)

    if key not in _REGISTER:
        raise ParserError(
            f"Invalid shielding function type {name} or reaction {reaction} is not supported"
        )

    return _REGISTER[key]


__all__ = [
    "_register",
    "ShieldingFunction",
]
