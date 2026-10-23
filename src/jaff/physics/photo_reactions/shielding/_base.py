"""Base interface for a line-shielding function plugin."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

from pydantic import ValidationError
from sympy import Expr

from ....config.schema import ShieldingParams
from ....errors import ParserError

if TYPE_CHECKING:
    from ....core.network import Network
    from ....core.reaction import Reaction


class ShieldingFunction(ABC):
    """A single photo-reaction line-shielding model.

    Subclasses live in their own module under ``shielding`` and register
    themselves via the :func:`~.register` decorator, keyed by their ``name``
    (the ``shielding.type`` value used in the network config).  The lookup is a
    dict, so resolving a type to its function is O(1) and independent of import
    order or any directory layout.

    Class attributes
    ----------------
    name : str
        Shielding-type identifier, matched case-insensitively against the
        ``shielding.type`` config value.
    reaction : str | None
        Serialized reaction this function is bound to (a *local* model), or
        ``None`` for a *global* model that applies to any reaction.  A function
        is resolved by the pair ``(reaction, name)`` with a fallback to
        ``(None, name)``, so a local model only matches its own reaction.

    A model with settings overrides :meth:`config_model` to declare them.
    """

    name: str
    reaction: str | None = None

    @staticmethod
    def config_model() -> type[ShieldingParams]:
        """Return the model of this shielding type's settings (keys but ``type``).

        The default has no keys, so a shielding type without settings rejects
        any.
        """
        return ShieldingParams

    def params(self, reaction: "Reaction") -> ShieldingParams:
        """Validate *reaction*'s shielding settings against :meth:`config_model`.

        Raises
        ------
        ParserError
            If a setting is unknown or has an invalid value.
        """
        assert reaction.photo is not None and reaction.photo.shielding is not None

        try:
            return self.config_model().model_validate(reaction.photo.shielding.props)
        except ValidationError as e:
            raise ParserError(
                f"Invalid {self.name} shielding settings for {reaction}:\n{e}"
            ) from e

    @abstractmethod
    def get_shielding(self, reaction: "Reaction", network: "Network") -> Expr:
        """Return the dimensionless shielding factor for *reaction*."""
        ...
