# ABOUTME: EmissivityFunction: base of the emissivity plugins (a spectrum in ph_E,
# ABOUTME: erg cm^-3 s^-1 eV^-1), selected by name, summed with the .jfunc globalEmissivity
"""Base interface of an emissivity plugin."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

from sympy import Expr

from .config import EmissivityParams

if TYPE_CHECKING:
    from ....core.network import Network


class EmissivityFunction(ABC):
    """A source of the network's volume emissivity.

    Radiation energy emitted per unit volume, time and photon energy, integrated over
    solid angle.

    A plugin subclasses it, sets ``name`` (the ``emissivity`` entry's ``type``) and
    registers itself with :func:`~jaff.physics.photo_reactions.emissivity.register`.  It
    is used only when an entry of ``[network.radiation] emissivity`` selects it; every
    entry's spectrum is added to the ``.jfunc`` ``globalEmissivity`` before the total is
    integrated over each radiation group.

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
    def config_model() -> type[EmissivityParams]:
        """Return the model of this type's settings (its table's keys but ``type``).

        The default has no keys, so a type without settings rejects any.
        """
        return EmissivityParams

    @abstractmethod
    def emissivity(self, net: Network, params: EmissivityParams) -> Expr:
        """The emissivity spectrum (erg cm⁻³ s⁻¹ eV⁻¹) as an expression in ``ph_E`` (eV).

        Parameters
        ----------
        net : Network
            The network; build the expression from its symbols (``net.symbols``,
            ``net.symbols.ph_E`` for the photon energy).
        params : EmissivityParams
            This entry's settings, typed by :meth:`config_model` (one type may be
            listed several times with other settings).

        Returns
        -------
        Expr
        """
