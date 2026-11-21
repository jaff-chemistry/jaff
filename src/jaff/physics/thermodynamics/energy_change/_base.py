# ABOUTME: EnergyChangeFunction: base of the per-reaction energy-change plugins (energy one
# ABOUTME: reaction event adds to the gas, erg; + heats, - cools), selected by name
"""Base interface of an energy-change plugin."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

from sympy import Expr

from .config import EnergyChangeParams

if TYPE_CHECKING:
    from ....core.reaction import Reaction


class EnergyChangeFunction(ABC):
    """A model of the energy one reaction event adds to the gas.

    A plugin subclasses it, sets ``name`` (the ``energy_change`` table's ``type``) and
    ``reaction`` (the serialized reaction it is bound to, or ``None`` for any reaction),
    and registers itself with
    :func:`~jaff.physics.thermodynamics.energy_change.register`.  It is used only when
    a reaction's ``thermodynamics.energy_change`` table selects it, and a ``.jfunc``
    ``deltaE<N>`` takes preference over it.

    Class attributes
    ----------------
    name : str
        The ``type`` value selecting it (case-insensitive).
    reaction : str | None
        The serialized reaction it is bound to; ``None`` for any reaction.

    A model with settings overrides :meth:`config_model` to declare them.
    """

    name: str
    reaction: str | None = None

    @staticmethod
    def config_model() -> type[EnergyChangeParams]:
        """Return the model of this type's settings (its table's keys but ``type``).

        The default has no keys, so a type without settings rejects any.
        """
        return EnergyChangeParams

    def params(self, reaction: Reaction) -> EnergyChangeParams:
        """Return *reaction*'s settings for this model, typed by :meth:`config_model`.

        They were validated when the config was read.
        """
        thermodynamics = reaction.config.thermodynamics
        assert thermodynamics is not None and thermodynamics.energy_change is not None

        return thermodynamics.energy_change.params

    @abstractmethod
    def energy_change(self, reaction: Reaction) -> Expr:
        """The energy one event of *reaction* adds to the gas (erg; negative cools).

        An expression in the network's symbols (e.g. ``tgas``, densities).
        """
