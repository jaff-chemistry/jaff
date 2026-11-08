# ABOUTME: net.analysis: numeric values for a network's unknowns (parameters,
# ABOUTME: species densities, unknown functions); enable/disable, pop, get_spec
"""Analysis mode: one place for the numeric values of a network's unknowns."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from types import MappingProxyType
from typing import TYPE_CHECKING, Any

import numpy as np
from sympy import Indexed, Symbol
from sympy.core.function import AppliedUndef

from .tables import Densities, Entry, Functions, Parameters

if TYPE_CHECKING:
    from ..network import Network

__all__ = ["Analysis", "AnalysisSpec"]


@dataclass(frozen=True)
class AnalysisSpec:
    """Read-only snapshot of the analysis tables.

    Attributes
    ----------
    parameters : mapping of str to Entry
        Entry of each parameter.
    densities : mapping of str to Entry
        Entry of each species.
    functions : mapping of str to Entry
        Entry of each unknown function.
    """

    parameters: Mapping[str, Entry]
    densities: Mapping[str, Entry]
    functions: Mapping[str, Entry]


class Analysis:
    """Numeric values for a network's unknowns, as ``net.analysis``.

    Disabled until :meth:`enable`, which builds three tables from the network:
    :attr:`parameters` (the symbols of the rate coefficients), :attr:`densities`
    (by species name) and :attr:`functions` (the unknown functions). Each entry
    holds a value and an optional legend label.

    Parameters
    ----------
    network : Network
        The network whose unknowns the tables hold.

    Examples
    --------
    >>> from jaff import Network
    >>> net = Network("networks/demos/water_co.jet")
    >>> net.analysis.enable(parameters={"tgas": 100.0, "av": 1.0})
    >>> net.analysis.parameters.values()
    {'av': 1.0, 'tgas': 100.0}
    """

    def __init__(self, network: Network):
        self.network = network
        self.parameters: Parameters | None = None
        self.densities: Densities | None = None
        self.functions: Functions | None = None

    @property
    def enabled(self) -> bool:
        """Whether analysis mode is on (the tables exist)."""
        return self.parameters is not None

    def enable(
        self,
        parameters: Mapping[str, Any] | None = None,
        densities: Mapping[str, Any] | Sequence[float] | np.ndarray | None = None,
        functions: Mapping[str, Any] | None = None,
    ) -> None:
        """Turn analysis mode on and set the given entries.

        Builds the tables from the network when disabled, then updates them,
        keeping entries set before.

        Parameters
        ----------
        parameters : mapping, optional
            Per parameter, a value or ``{"value": ..., "label": ...}``.
        densities : mapping or array, optional
            Per species name, a value or ``{"value": ..., "label": ...}``; or an
            array in species order.
        functions : mapping, optional
            Per unknown function, a value (number, expression of parameters, or
            callable without arguments) or ``{"value": ..., "label": ...}``.

        Raises
        ------
        KeyError
            If a name is not one of the network's.
        ValueError
            If a value is not allowed.
        """
        if not self.enabled:
            self._build()
        assert (
            self.parameters and self.densities is not None and self.functions is not None
        )
        self.parameters.update(parameters or {})
        if densities is not None:
            self.densities.update(densities)
        self.functions.update(functions or {})

    def disable(self) -> None:
        """Turn analysis mode off, dropping the tables and their values."""
        self.parameters = self.densities = self.functions = None

    def pop(self) -> None:
        """Clear every entry of every table; analysis mode stays on."""
        for table in self.tables():
            table.pop()

    def get_spec(self) -> AnalysisSpec | None:
        """Read-only snapshot of the tables, or ``None`` when disabled."""
        if not self.enabled:
            return None
        parameters, densities, functions = (
            MappingProxyType(t.get()) for t in self.tables()
        )

        return AnalysisSpec(parameters, densities, functions)

    def tables(self) -> tuple[Parameters, Densities, Functions]:
        """The parameter, density and function tables.

        Raises
        ------
        ValueError
            If analysis mode is off.
        """
        if self.parameters is None or self.densities is None or self.functions is None:
            raise ValueError("Analysis mode is off: call net.analysis.enable() first")

        return self.parameters, self.densities, self.functions

    def _build(self) -> None:
        symbols: set[str] = set()
        functions: set[str] = set()
        for reaction in self.network.reactions:
            expr = reaction.rate_coefficient
            indexed = {a.base.label for a in expr.atoms(Indexed)}
            symbols |= {
                s.name for s in expr.free_symbols - indexed if isinstance(s, Symbol)
            }
            functions |= {f.func.__name__ for f in expr.atoms(AppliedUndef)}
        species = self.network.species
        self.parameters = Parameters(sorted(symbols))
        self.densities = Densities(
            [s.name for s in species], {s.name: s.latex() for s in species}
        )
        self.functions = Functions(sorted(functions), self.parameters)
