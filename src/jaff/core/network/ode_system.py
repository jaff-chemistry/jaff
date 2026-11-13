# ABOUTME: OdeSystem: the state vector a generated solver integrates and its symbolic
# ABOUTME: right-hand side and Jacobian, assembled from the network's components

"""The ODE system of a :class:`~jaff.core.Network`: its state, rhs and Jacobian."""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING, Literal, cast

from sympy import Basic, Expr, ImmutableMatrix, ImmutableSparseMatrix, Indexed

from ...common.derivative import derivative
from ...io import jaff_progress
from ...physics import Radiation
from ..reaction.rate import check_form

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
    """The ODE system a generated solver integrates: its state, rhs and Jacobian.

    Built by :meth:`~jaff.core.Network.ode_system`.  The state holds every
    species density, then the thermal variable (none, temperature or internal
    energy), then two moments per radiation band laid out by *rad_order* (see
    :meth:`~jaff.physics.Radiation.ordered_index`).  :attr:`rhs` has one row per
    state entry, in the same order, taken from the components that own the
    physics; :attr:`jacobian` is its derivative with respect to the state.

    *form* chooses how the reaction rate coefficients appear (see
    :class:`~jaff.core.reaction.Rate`): ``"expanded"`` writes each one out,
    ``"reduced"`` writes reaction *r*'s as the symbol ``krate[r]``, defined in
    :attr:`rate_coefficients`, so that generated code evaluates each coefficient
    once and shares it between the rhs and the Jacobian.

    Examples
    --------
    >>> from jaff import Network
    >>> net = Network("h_photoionization", radiation={"bands": [13.6, 20.0]})
    >>> system = net.ode_system(thermal="dtdt", radiation=True)
    >>> [var.kind for var in system.state]
    ['species', 'species', 'species', 'temperature', 'radiation_density', 'radiation_flux']
    >>> system.rhs.shape
    (6, 1)

    In the reduced form the rows hold the ``krate`` symbols, and the Jacobian's
    temperature column takes ``∂k/∂T`` by the chain rule:

    >>> system = net.ode_system(thermal="dtdt", form="reduced")
    >>> system.rhs[0]
    -krate[0]*nden[0] + krate[1]*nden[1]*nden[2]
    >>> system.rate_coefficients[net.symbols.krates[1]]
    2.63e-13/(tgas/10000.0)**0.7
    >>> system.jacobian[0, 3]
    -1.16159247118804e-10*nden[1]*nden[2]/tgas**1.7
    """

    def __init__(
        self,
        network: Network,
        thermal: str = "none",
        energy: str = "volumetric",
        radiation: bool = False,
        rad_order: int = 0,
        form: str = "expanded",
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
        form : {"expanded", "reduced"}, optional
            How the rate coefficients appear in the equations: written out
            (``"expanded"``, default) or as their ``krate[r]`` symbols
            (``"reduced"``, defined in :attr:`rate_coefficients`).

        Raises
        ------
        ValueError
            If *thermal* or (with ``"dedt"``) *energy* is not a valid choice, if
            *radiation* is requested on a network without a radiation field, or
            if *rad_order* is not a valid layout, or if *form* is unknown.
        """
        check_form(form)
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
        self.form: str = form
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

        Species rows are :meth:`~jaff.core.Kinetics.odes`; the thermal row is
        :meth:`~jaff.physics.Thermodynamics.dTdt_tot` or ``dEdt_tot`` in the
        *energy* form; the radiation rows are :meth:`~jaff.physics.Radiation.odes`
        in *rad_order*; all in :attr:`form`.
        """
        net = self.network
        form = self.form
        rows: list[Basic] = net.kinetics.odes(form)
        thermo = net.thermodynamics
        if self.thermal == "dedt":
            rows.append(thermo.dEdt_tot(form).normaliser(self.energy))
        elif self.thermal == "dtdt":
            rows.append(thermo.dTdt_tot(form))
        if self.radiation:
            assert net.radiation is not None  # checked in __init__
            rows.extend(net.radiation.odes(self.rad_order, form))

        return ImmutableMatrix(rows)

    @cached_property
    def rate_coefficients(self) -> dict[Indexed, Expr]:
        """Definitions of the ``krate[r]`` symbols of the reduced form.

        ``{krate[r]: rate coefficient of reaction r}`` for every reaction, in
        catalogue order; empty in the expanded form, whose equations hold no
        ``krate`` symbol.
        """
        if self.form == "expanded":
            return {}

        return {r.rate.symbol: r.rate.coefficient for r in self.network.reactions}

    @cached_property
    def parameters(self) -> tuple[Basic, ...]:
        """The system's free parameters: the network's variables that are not state.

        The network's free symbols
        (:attr:`~jaff.core.network.NetworkSymbols.variables`, across rates, energy
        and radiation sources) minus this system's state variables, sorted by name; a
        generated solver needs a value for each.  ``tgas`` is a parameter only when
        the thermal mode is ``"none"``; otherwise it is state.  The list is
        network-wide: a symbol used only by the thermal or radiation equations is
        listed even when the system leaves those out.

        Examples
        --------
        >>> from jaff import Network
        >>> net = Network("COthin")
        >>> net.ode_system().parameters
        (av, crate, tdust, tgas)
        >>> net.ode_system(thermal="dtdt").parameters
        (av, crate, tdust)
        """
        state = {var.symbol for var in self.state}
        names = (s for s in self.network.symbols.variables if s not in state)

        return tuple(sorted(names, key=str))

    @cached_property
    def _rate_derivatives(self) -> dict[Indexed, dict[Basic, Expr]]:
        """``∂k_r/∂y`` per ``krate[r]``, for the state symbols *y* that k_r contains.

        Found from each coefficient's free symbols (unknown-function arguments
        included); a coefficient free of the state gets an empty entry.  Each
        derivative is taken once and shared by every row holding ``krate[r]``,
        written with ``krate[r]`` itself where cheaper
        (:func:`~jaff.common.derivative.derivative`).
        """
        state = {var.symbol for var in self.state}

        return {
            k: {
                y: derivative(coefficient, y, value=k)
                for y in coefficient.free_symbols & state
            }
            for k, coefficient in self.rate_coefficients.items()
        }

    def _rate_terms(self, row: Expr) -> list[tuple[dict[Basic, Expr], Expr]]:
        """``(∂k_r/∂y by y, ∂row/∂k_r)`` for each state-dependent ``krate[r]`` in *row*."""
        derivatives = self._rate_derivatives
        if not derivatives:
            return []

        return [
            (derivatives[k], row.diff(k))
            for k in row.atoms(Indexed)
            if derivatives.get(k)
        ]

    @staticmethod
    def _total_diff(
        row: Expr, y: Basic, terms: list[tuple[dict[Basic, Expr], Expr]]
    ) -> Expr:
        """``∂row/∂y |_k + Σ_r (∂row/∂k_r)(∂k_r/∂y)``: the chain rule over the rates."""
        derivative = row.diff(y)
        for dk_dy, drow_dk in terms:
            if y in dk_dy:
                derivative += drow_dk * dk_dy[y]

        return derivative

    @cached_property
    def jacobian(self) -> ImmutableSparseMatrix:
        """Derivative of :attr:`rhs` with respect to the state (n × n, sparse).

        Entry ``(i, j)`` is ``∂ rhs[i] / ∂ state[j]``; only non-zero entries are
        stored.  Species and radiation columns are plain derivatives.  The
        thermal column is ``∂rhs/∂T`` with ``"dtdt"``; with ``"dedt"`` the
        equations are written in ``T`` but the state holds the energy ``e``,
        so the chain rule through the EOS
        (:meth:`~jaff.physics.Thermodynamics.temperature_derivatives`) gives the
        energy column ``(∂rhs/∂T)(∂T/∂e)`` and adds ``(∂rhs/∂T)(∂T/∂n_j)`` to
        each species column.

        In the reduced form a row depends on the state also through the
        ``krate[r]`` it holds, so every derivative takes the chain rule over them:

        ```
        ∂rhs_i/∂y_j = ∂rhs_i/∂y_j |_k + Σ_r (∂rhs_i/∂k_r) · (∂k_r/∂y_j)
        ```

        summed over the reactions whose coefficient contains ``y_j`` (most often
        ``tgas``; a photo rate also holds the radiation densities).  Entries keep
        the ``krate`` symbols; ``∂k_r/∂y_j`` is ``krate[r] · ∂ln k_r/∂y_j`` where
        that is cheaper than writing it out
        (:func:`~jaff.common.derivative.derivative`).
        """
        rows = cast(list[Expr], self.rhs.flat())
        columns = {
            var.symbol: j
            for j, var in enumerate(self.state)
            if var.kind not in ("energy", "temperature")
        }
        entries: dict[tuple[int, int], Basic] = {}
        row_terms = []
        for i, row in enumerate(
            jaff_progress.track(rows, description="Differentiating the rhs")
        ):
            terms = self._rate_terms(row)
            row_terms.append(terms)
            # Only the state symbols a row contains, directly or through a rate
            # coefficient, give non-zero entries
            for sym in row.free_symbols.union(*(dk_dy for dk_dy, _ in terms)):
                j = columns.get(sym)
                if j is not None:
                    entries[i, j] = self._total_diff(row, sym, terms)
        if self.thermal != "none":
            self._add_thermal_column(rows, row_terms, entries)

        n = len(self.state)
        return ImmutableSparseMatrix(
            n, n, {key: expr for key, expr in entries.items() if expr != 0}
        )

    def _add_thermal_column(
        self,
        rows: list[Expr],
        row_terms: list[list[tuple[dict[Basic, Expr], Expr]]],
        entries: dict[tuple[int, int], Basic],
    ) -> None:
        """Add the thermal column (and the ``dedt`` species corrections) to *entries*."""
        net = self.network
        t = net.species.count  # the thermal slot follows the species
        tgas = net.symbols.tgas
        df_dT = [
            self._total_diff(row, tgas, terms)
            for row, terms in zip(
                jaff_progress.track(rows, description="Differentiating by T"),
                row_terms,
                strict=True,
            )
        ]
        if self.thermal == "dtdt":
            for i, d in enumerate(df_dT):
                entries[i, t] = d
            return

        dT_de, dT_dn = net.thermodynamics.temperature_derivatives(self.energy)
        for i, d in enumerate(df_dT):
            if d == 0:
                continue
            entries[i, t] = d * dT_de
            for j, dT_dn_j in enumerate(dT_dn):
                entries[i, j] = entries.get((i, j), 0) + d * dT_dn_j
