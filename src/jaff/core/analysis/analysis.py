# ABOUTME: net.analysis: numeric values for a network's unknowns (parameters,
# ABOUTME: species densities, unknown functions); enable/disable, evaluate, get_spec
"""Analysis mode: one place for the numeric values of a network's unknowns."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, overload

import numpy as np
from sympy import Basic, Indexed, MatrixBase, Symbol, sympify
from sympy.core.function import AppliedUndef

from ._evaluate import evaluate
from .tables import Densities, Entry, Functions, FunctionValue, Parameters

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

    :meth:`evaluate` gives the value of any of the network's expressions (a rate
    coefficient, the rates, an ODE system's rhs or Jacobian, ...), with values given
    to it merged over the tables'.

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
    >>> net.analysis.evaluate(net.reactions[2].rate.coefficient)
    4.382088543149259e-07
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

    @overload
    def evaluate(
        self,
        exprs: Basic | float,
        parameters: Mapping[str, Any] | None = None,
        densities: Mapping[str, Any] | Sequence[float] | np.ndarray | None = None,
        functions: Mapping[str, FunctionValue] | None = None,
    ) -> float: ...

    @overload
    def evaluate(
        self,
        exprs: Sequence[Basic | float] | MatrixBase,
        parameters: Mapping[str, Any] | None = None,
        densities: Mapping[str, Any] | Sequence[float] | np.ndarray | None = None,
        functions: Mapping[str, FunctionValue] | None = None,
    ) -> np.ndarray: ...

    def evaluate(
        self,
        exprs: Any,
        parameters: Mapping[str, Any] | None = None,
        densities: Mapping[str, Any] | Sequence[float] | np.ndarray | None = None,
        functions: Mapping[str, FunctionValue] | None = None,
    ) -> float | np.ndarray:
        """Evaluate *exprs* numerically, with the given values over the tables'.

        Every symbol needs a value, a density ``nden[i]`` the density of species *i*,
        and every unknown function (one the ``.jfunc`` file uses but does not define,
        such as an interpolated shielding table) a value.  Each argument is merged
        over the matching table when analysis mode is on; without it only the given
        values are used.  Values the expressions do not use are ignored.

        Parameters
        ----------
        exprs : Expr, sequence of Expr, or Matrix
            What to evaluate: one expression (e.g. ``reaction.rate.coefficient``),
            a sequence of them (e.g. ``net.reactions.rates()``) or a matrix (e.g.
            an :class:`~jaff.core.OdeSystem`'s ``rhs`` or ``jacobian``).  Rates
            and equations must be in the expanded form: a ``krate[r]`` symbol has
            no value.
        parameters : mapping, optional
            A number per parameter symbol, e.g. ``{"tgas": 100.0, "av": 1.0}``.
        densities : mapping or array, optional
            Number density (cm⁻³) per species name, or an array in species order.
        functions : mapping, optional
            Per unknown function name, a number, an expression of the parameters
            or a callable without arguments returning either; the function's
            arguments in the network are not used.

        Returns
        -------
        float or numpy.ndarray
            A float for one expression, a 1-D array for a sequence, an array of
            the matrix's shape for a matrix (zeros where a sparse matrix stores
            nothing).

        Raises
        ------
        ValueError
            Naming every parameter, density and function without a value, all
            expressions together; or if an expression holds ``krate`` symbols.

        Examples
        --------
        >>> from jaff import Network
        >>> net = Network("networks/demos/water_co.jet")
        >>> net.analysis.evaluate(net.reactions[2].rate.coefficient, {"tgas": 100.0})
        4.382088543149259e-07
        >>> parameters = {"tgas": 100.0, "crate": 1.3e-17, "av": 1.0, "chi": 1.0}
        >>> net.analysis.evaluate(net.reactions.rate_coefficients(), parameters)[
        ...     [0, 2, 16]
        ... ]
        array([1.26100000e-17, 4.38208854e-07, 8.14930913e-12])
        >>> net.analysis.evaluate(net.reactions.rate_coefficients(), {"tgas": 100.0})
        Traceback (most recent call last):
        ...
        ValueError: Cannot evaluate: no value for parameter av, chi, crate
        """
        values = (parameters, densities, functions)
        if isinstance(exprs, MatrixBase):
            entries = exprs.todok()
            result = np.zeros(exprs.shape)
            numbers = evaluate(self.network, list(entries.values()), *values)
            for (i, j), number in zip(entries, numbers):
                result[i, j] = number

            return result

        if isinstance(exprs, (Basic, int, float)):
            return float(evaluate(self.network, [sympify(exprs)], *values)[0])

        return evaluate(self.network, [sympify(e) for e in exprs], *values)

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
            expr = reaction.rate.coefficient
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
