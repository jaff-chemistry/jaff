# ABOUTME: OdeSystem: the state vector a generated solver integrates and its symbolic
# ABOUTME: right-hand side, assembled from the network's components

"""The ODE system of a :class:`~jaff.core.Network`: its state and right-hand side."""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING, Literal

from sympy import Basic, ImmutableMatrix

from ...physics import Radiation

if TYPE_CHECKING:
    from .network import Network

StateKind = Literal[
    "species", "temperature", "energy", "radiation_density", "radiation_flux"
]

THERMAL_MODES: tuple[str, ...] = ("none", "dedt", "dtdt")
"""The thermal equations an :class:`OdeSystem` can carry: none, ``dE/dt`` or ``dT/dt``."""


@dataclass(frozen=True)
class StateVariable:
    """One entry of an :class:`OdeSystem` state vector.

    Attributes
    ----------
    kind : str
        ``"species"``, ``"temperature"``, ``"energy"``, ``"radiation_density"``
        or ``"radiation_flux"``.
    symbol : sympy.Basic
        The symbol the equations are written in: ``nden[i]``, ``tgas`` (for both
        thermal kinds: the energy equation is written in temperature),
        ``photden[i]`` / ``radeden[i]`` or ``rflux[i]``.
    index : int
        The species index or band index (``0`` for the thermal entry).
    """

    kind: StateKind
    symbol: Basic
    index: int


class OdeSystem:
    """The ODE system a generated solver integrates: a state vector and its rhs.

    Built by :meth:`~jaff.core.Network.ode_system`.  The state holds every
    species density, then the thermal variable (none, temperature or internal
    energy), then two moments per radiation band laid out by *rad_order* (see
    :meth:`~jaff.physics.Radiation.ordered_index`).  :attr:`rhs` has one row per
    state entry, in the same order, taken from the components that own the
    physics.

    Examples
    --------
    >>> from jaff import Network
    >>> net = Network("h_photoionization", radiation={"bands": [13.6, 20.0]})
    >>> system = net.ode_system(thermal="dtdt", radiation=True)
    >>> [var.kind for var in system.state]
    ['species', 'species', 'species', 'temperature', 'radiation_density', 'radiation_flux']
    >>> system.rhs.shape
    (6, 1)
    """

    def __init__(
        self,
        network: Network,
        thermal: str = "none",
        energy: str = "volumetric",
        radiation: bool = False,
        rad_order: int = 0,
    ) -> None:
        """Validate the choices of equations and lay out the state.

        Parameters
        ----------
        network : Network
            The network whose equations the system holds.
        thermal : str, optional
            Thermal equation: ``"none"`` (default), ``"dedt"`` (internal energy,
            in the *energy* form) or ``"dtdt"`` (temperature).
        energy : str, optional
            Form of the evolved internal energy with ``"dedt"``:
            ``"volumetric"`` (default), ``"specific"``, ``"per_particle"`` or
            ``"molar"`` (see
            :meth:`~jaff.physics.thermodynamics.InternalEnergy.normaliser`).
        radiation : bool, optional
            Include the radiation-moment equations.  Default ``False``.
        rad_order : int, optional
            Layout of the radiation equations, ``0`` (default) to ``3``; used
            only with *radiation*.

        Raises
        ------
        ValueError
            If *thermal* or (with ``"dedt"``) *energy* is not a valid choice, if
            *radiation* is requested on a network without a radiation field, or
            if *rad_order* is not a valid layout.
        """
        if thermal not in THERMAL_MODES:
            raise ValueError(
                f"Invalid thermal mode {thermal!r}; "
                f"valid modes are: {', '.join(THERMAL_MODES)}"
            )
        if thermal == "dedt":
            network.thermodynamics.eos.normaliser(energy)  # validates the form
        if radiation and network.radiation is None:
            raise ValueError("radiation=True, but the network has no radiation field")
        if radiation and rad_order not in Radiation.ORDERS:
            raise ValueError(
                f"Invalid radiation order {rad_order!r}; "
                f"supported orders are: {', '.join(map(str, Radiation.ORDERS))}"
            )

        self.network: Network = network
        self.thermal: str = thermal
        self.energy: str = energy
        self.radiation: bool = radiation
        self.rad_order: int = rad_order
        self.state: tuple[StateVariable, ...] = self._layout()

    def _layout(self) -> tuple[StateVariable, ...]:
        """The state vector: species, thermal slot, radiation by *rad_order*."""
        net = self.network
        ndens = net.symbols.ndens
        state = [StateVariable("species", ndens[i], i) for i in range(net.species.count)]
        if self.thermal != "none":
            kind: StateKind = "energy" if self.thermal == "dedt" else "temperature"
            state.append(StateVariable(kind, net.symbols.tgas, 0))
        if not self.radiation:
            return tuple(state)

        rad = net.radiation
        assert rad is not None  # checked in __init__
        moments: list[StateVariable | None] = [None] * (2 * rad.nbands)
        for i in range(rad.nbands):
            ei, fi = rad.ordered_index(i, self.rad_order)
            moments[ei] = StateVariable("radiation_density", rad.den[i], i)
            moments[fi] = StateVariable("radiation_flux", rad.flux[i], i)

        return tuple(state) + tuple(m for m in moments if m is not None)

    @cached_property
    def rhs(self) -> ImmutableMatrix:
        """Time derivative of the state, one row per :attr:`state` entry (n × 1).

        Species rows are :meth:`~jaff.core.Reactions.odes`; the thermal row is
        :attr:`~jaff.physics.Thermodynamics.dTdt_tot` or ``dEdt_tot`` in the
        *energy* form; the radiation rows are :meth:`~jaff.physics.Radiation.odes`
        in *rad_order*.
        """
        net = self.network
        rows: list[Basic] = list(net.reactions.odes(net.species))
        thermo = net.thermodynamics
        if self.thermal == "dedt":
            rows.append(thermo.dEdt_tot.normaliser(self.energy))
        elif self.thermal == "dtdt":
            rows.append(thermo.dTdt_tot)
        if self.radiation:
            assert net.radiation is not None  # checked in __init__
            rows.extend(net.radiation.odes(self.rad_order))

        return ImmutableMatrix(rows)
