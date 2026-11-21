# ABOUTME: Registry: named, opt-in plugins of one kind keyed by (name, reaction), a
# ABOUTME: reaction's own plugin beating a global one; plugins imported on first lookup
"""A registry of named plugins of one kind (e.g. shielding functions)."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Generic, TypeVar

from ..errors import RegistrationError
from . import plugins
from ._helper import import_subpackages

T = TypeVar("T")

Key = tuple[str, str | None]
"""``(name, reaction)``: a plugin's lower-cased name and the serialized reaction it is
bound to, ``None`` for a global plugin (any reaction)."""


class Registry(Generic[T]):
    """Plugins of one kind, registered by a class decorator and found by name.

    A plugin class declares ``name`` (the value a config file selects it by,
    case-insensitive) and ``reaction``: the serialized reaction it is bound to, or
    ``None`` for a global plugin serving any reaction.  One instance is registered per
    class.  :meth:`find` prefers a reaction's own plugin and falls back to the global
    one of the same name.  The plugin packages, then the repo-root ``plugins`` folder
    (:func:`~jaff.common.plugins.load_plugins`), are imported once, on the first
    lookup, so their ``register`` decorators run.

    Parameters
    ----------
    kind : str
        What the plugins are, for messages (e.g. ``"shielding function"``).
    packages : Sequence[str]
        Packages whose subpackages hold the plugins (imported on first use).

    Examples
    --------
    >>> from jaff.common.registry import Registry
    >>> models = Registry("model", packages=[])
    >>> @models.register
    ... class Simple:
    ...     name, reaction = "simple", None
    >>> models.find("SIMPLE", "H2._PHOTON__H.H").name
    'simple'
    """

    def __init__(self, kind: str, packages: Sequence[str]) -> None:
        self.kind = kind
        self._packages = list(packages)
        self._plugins: dict[Key, T] = {}
        self._discovered = False

    def register(self, cls: type[Any]) -> type[Any]:
        """Class decorator registering one instance of *cls* under ``(name, reaction)``.

        Parameters
        ----------
        cls : type
            The plugin class, with ``name`` and ``reaction`` attributes.

        Returns
        -------
        type
            *cls*, unchanged.

        Raises
        ------
        RegistrationError
            If a plugin with the same ``(name, reaction)`` is already registered.
        """
        key = (cls.name.lower(), cls.reaction)
        if key in self._plugins:
            raise RegistrationError(
                f"Key already exists for {self.kind}: {cls.name}, {cls.reaction}"
            )

        self._plugins[key] = cls()

        return cls

    def find(self, name: str, reaction: str | None) -> T | None:
        """The plugin *name* for the serialized *reaction*, or ``None``.

        A plugin bound to *reaction* wins over a global one of the same name; a plugin
        bound to another reaction never applies.

        Parameters
        ----------
        name : str
            The plugin's name (case-insensitive).
        reaction : str or None
            The serialized reaction, or ``None`` for a global plugin only.

        Returns
        -------
        object or None
        """
        self._discover()
        own = self._plugins.get((name.lower(), reaction))

        return own if own is not None else self._plugins.get((name.lower(), None))

    def items(self) -> dict[Key, T]:
        """Every registered plugin, keyed by ``(name, reaction)`` (a copy)."""
        self._discover()

        return dict(self._plugins)

    def registered(self) -> dict[Key, T]:
        """The plugins registered so far, without importing the plugin packages (a copy).

        For checks made while plugins register themselves, when importing the
        packages would re-enter them.
        """
        return dict(self._plugins)

    def _discover(self) -> None:
        """Import the plugin packages once, so their ``register`` decorators run."""
        if self._discovered:
            return

        self._discovered = True
        for package in self._packages:
            import_subpackages(package)
        plugins.load_plugins()
