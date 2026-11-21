# ABOUTME: AbsorptionCoefficientFunction: base of the absorption coefficient plugins (a
# ABOUTME: spectrum in ph_E, cm^-1), selected by name, summed with the .jfunc one
"""Base interface of an absorption coefficient plugin."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

from sympy import Expr

from .config import AbsorptionCoefficientParams

if TYPE_CHECKING:
    from ....core.network import Network


class AbsorptionCoefficientFunction(ABC):
    """An absorber: part of the network's absorption coefficient.

    The fraction of the radiation absorbed per unit path length.  A plugin subclasses
    it, sets ``name`` (the ``absorption_coefficient`` entry's ``type``) and registers
    itself with ``register`` (:mod:`jaff.physics.photo_reactions.absorption_coefficient`).
    It is used only when an entry of ``[network.radiation] absorption_coefficient``
    selects it; every entry's spectrum is added to the ``.jfunc``
    ``globalAbsorptionCoefficient`` before the total is averaged over each radiation
    group.

    Class attributes
    ----------------
    name : str
        The ``type`` value selecting it (case-insensitive).
    reaction : None
        Always ``None``: the source belongs to no reaction.

    A model with settings overrides :meth:`config_model` to declare them.
    """

    name: str
    reaction: None = None

    @staticmethod
    def config_model() -> type[AbsorptionCoefficientParams]:
        """Return the model of this type's settings (its table's keys but ``type``).

        The default has no keys, so a type without settings rejects any.
        """
        return AbsorptionCoefficientParams

    @abstractmethod
    def absorption_coefficient(
        self, net: Network, params: AbsorptionCoefficientParams
    ) -> Expr:
        """The absorption coefficient spectrum (cm⁻¹) as an expression in ``ph_E`` (eV).

        Parameters
        ----------
        net : Network
            The network; build the expression from its symbols (``net.symbols``,
            ``net.symbols.ph_E`` for the photon energy).
        params : AbsorptionCoefficientParams
            This entry's settings, typed by :meth:`config_model` (one type may be
            listed several times with other settings).

        Returns
        -------
        Expr
        """
