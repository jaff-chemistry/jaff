# ABOUTME: HeatingCoolingRateFunction: base of the global heating/cooling-rate plugins (a rate
# ABOUTME: in erg cm^-3 s^-1; + heats, - cools), selected by name, summed with the .jfunc one
"""Base interface of a global heating/cooling-rate plugin."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

from sympy import Expr

from .config import HeatingCoolingRateParams

if TYPE_CHECKING:
    from ....core.network import Network


class HeatingCoolingRateFunction(ABC):
    """A global (non-reactive) heating/cooling mechanism of the gas.

    A plugin subclasses it, sets ``name`` (the ``heating_cooling_rate`` entry's
    ``type``) and registers itself with
    :func:`~jaff.physics.thermodynamics.heating_cooling_rate.register`.  It is used only
    when an entry of ``[network.thermodynamics] heating_cooling_rate`` selects it; every
    entry's rate is added to the ``.jfunc`` ``heatingCoolingRate``.

    Class attributes
    ----------------
    name : str
        The ``type`` value selecting it (case-insensitive).
    reaction : None
        Always ``None``: the mechanism belongs to no reaction.

    A model with settings overrides :meth:`config_model` to declare them.
    """

    name: str
    reaction: None = None

    @staticmethod
    def config_model() -> type[HeatingCoolingRateParams]:
        """Return the model of this type's settings (its table's keys but ``type``).

        The default has no keys, so a type without settings rejects any.
        """
        return HeatingCoolingRateParams

    @abstractmethod
    def heating_cooling_rate(
        self, net: Network, params: HeatingCoolingRateParams
    ) -> Expr:
        """The rate the mechanism adds to the gas energy (erg cm⁻³ s⁻¹; negative cools).

        Parameters
        ----------
        net : Network
            The network; build the expression from its symbols (``net.symbols``).
        params : HeatingCoolingRateParams
            This entry's settings, typed by :meth:`config_model` (one type may be
            listed several times with other settings).

        Returns
        -------
        Expr
        """
