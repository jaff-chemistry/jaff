# ABOUTME: ShieldingFunction: abstract base of the line-shielding plugins,
# ABOUTME: registered by (name, reaction) and selected by shielding.type
"""Base interface for a line-shielding function plugin."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

from sympy import Expr

from .config import ShieldingParams

if TYPE_CHECKING:
    from ....core.reaction import Reaction


@dataclass(frozen=True)
class ShieldingResult:
    """What a shielding plugin computes for one reaction.

    Attributes
    ----------
    factor : sympy.Expr
        The dimensionless shielding factor multiplying the reaction's band
        rates.
    tables : dict[str, Any]
        Data tables the generated code interpolates over, keyed by file name
        (e.g. ``shielding_<reaction>.hdf5``); empty when the model needs none.
        Kept as ``reaction.photo.shielding_tables`` and written out by
        :meth:`~jaff.core.network.Network.write_shielding_tables`.
    """

    factor: Expr
    tables: dict[str, Any] = field(default_factory=dict)


class ShieldingFunction(ABC):
    """A single photo-reaction line-shielding model.

    Subclasses live in their own module under ``shielding`` and register
    themselves via the :func:`~jaff.physics.photo_reactions.shielding.register`
    decorator, keyed by their ``name``
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
        is resolved by the pair ``(name, reaction)`` with a fallback to
        ``(name, None)``, so a local model only matches its own reaction.

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
        """Return *reaction*'s shielding settings, typed by :meth:`config_model`.

        They were validated when the config was read
        (``reaction.config.photo.shielding.params``).
        """
        photo = reaction.photo
        assert photo is not None and photo.config.shielding is not None
        params = photo.config.shielding.params
        assert isinstance(params, self.config_model())

        return params

    @abstractmethod
    def get_shielding(self, reaction: "Reaction") -> ShieldingResult:
        """Return the shielding factor for *reaction* and any data tables it needs."""
        ...
