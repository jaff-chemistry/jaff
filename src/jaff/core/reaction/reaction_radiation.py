# ABOUTME: ReactionRadiation: every reaction's radiation component, held as Reaction.rad:
# ABOUTME: the radiation it emits (spectral yield) and that yield per radiation group

from __future__ import annotations

from typing import TYPE_CHECKING, List, cast

from sympy import Expr

if TYPE_CHECKING:
    from .reaction import Reaction


class ReactionRadiation:
    """A reaction's coupling to the radiation field, held as ``Reaction.rad``.

    Built by :class:`~jaff.core.reaction.Reaction` for every reaction (any
    reaction type can emit).  The emission is the reaction's ``.jfunc``
    ``spectralYield<N>`` function; with a radiation field, :meth:`compute`
    integrates it over each radiation group.  It has no settings of its own.

    Attributes
    ----------
    reaction : Reaction
        The reaction it belongs to.
    group_yields : list[Expr] or None
        :attr:`spectral_yield` integrated over each radiation group, in
        ``net.radiation.groups`` order, in the units of the group density per
        event (erg in ``"u"`` mode, photons in ``"nph"`` mode); see
        :meth:`~jaff.physics.Radiation.integrate_over_groups`.  ``None``
        without radiation or without a yield.

    Examples
    --------
    >>> from jaff import Network
    >>> net = Network("h_photoionization")
    >>> rad = net.reactions[0].rad
    >>> rad.spectral_yield is None, rad.group_yields is None
    (True, True)
    """

    def __init__(self, reaction: Reaction):
        """Bind the component to *reaction*; :meth:`compute` fills the rest.

        Parameters
        ----------
        reaction : Reaction
            The reaction.
        """
        self.reaction: Reaction = reaction
        self.group_yields: List[Expr] | None = None

    @property
    def spectral_yield(self) -> Expr | None:
        """Radiation energy the reaction emits per event per unit photon energy.

        In erg eV⁻¹, integrated over solid angle, as an expression in ``ph_E``
        (eV): the reaction's ``.jfunc`` ``spectralYield<N>`` function.

        Returns
        -------
        sympy.Expr or None
            The yield, or ``None`` when the reaction emits nothing.
        """
        function = self.reaction.aux_funcs.spectralyield

        return cast(Expr, function.expr) if function is not None else None

    def compute(self) -> None:
        """Integrate :attr:`spectral_yield` over the radiation groups.

        Sets :attr:`group_yields`; does nothing without a radiation field or
        without a yield.
        """
        radiation = self.reaction.net.radiation
        spectral_yield = self.spectral_yield
        if radiation is None or spectral_yield is None:
            return

        self.group_yields = radiation.integrate_over_groups(spectral_yield)
