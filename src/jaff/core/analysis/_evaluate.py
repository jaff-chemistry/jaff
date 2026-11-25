# ABOUTME: Numeric evaluation of rate expressions: parameters, species densities and
# ABOUTME: unknown-function values merged over net.analysis; missing names reported
"""Evaluate a network's rate expressions with numeric values."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, cast

import numpy as np
from sympy import Basic, Expr, Float, Indexed, Symbol, sympify
from sympy.core.function import AppliedUndef

from .tables import FunctionValue

if TYPE_CHECKING:
    from ..network import Network

#: The indexed symbol of the species number densities.
NDEN = "nden"

#: The indexed symbol of the rate coefficients of the reduced form.
KRATE = "krate"


@dataclass(frozen=True)
class Inputs:
    """The values one evaluation uses.

    Attributes
    ----------
    parameters : dict of str to float
        Value of each parameter symbol.
    densities : dict of int to float
        Number density of each species, by species index.
    functions : dict of str to Expr
        Value of each unknown function, as an expression of the parameters.
    names : dict of int to str
        Species name of each index, for messages.
    """

    parameters: dict[str, float]
    densities: dict[int, float]
    functions: dict[str, Expr]
    names: dict[int, str]


def _densities(net: Network, given: Any) -> dict[str, float]:
    """Densities by species name, from a mapping or an array in species order."""
    if given is None:
        return {}
    if isinstance(given, Mapping):
        return {str(k): float(v) for k, v in given.items()}
    values = np.asarray(given, dtype=float)
    species = [s.name for s in net.species]
    if values.shape != (len(species),):
        raise ValueError(
            f"A density array needs one value per species ({len(species)}), "
            f"got shape {values.shape}"
        )

    return dict(zip(species, values.tolist()))


def _function(name: str, value: FunctionValue) -> Expr:
    """The value of one unknown function as an expression without unknown functions."""
    given: Any = value() if callable(value) and not isinstance(value, Basic) else value
    expr = cast(Expr, sympify(given))
    unknown = sorted(f.func.__name__ for f in expr.atoms(AppliedUndef))
    if unknown:
        raise ValueError(
            f"function {name!r}: its value may not use unknown functions "
            f"({', '.join(unknown)})"
        )

    return expr


def resolve(
    net: Network,
    parameters: Mapping[str, Any] | None = None,
    densities: Mapping[str, Any] | Sequence[float] | np.ndarray | None = None,
    functions: Mapping[str, FunctionValue] | None = None,
) -> Inputs:
    """The values to evaluate with: the given ones over ``net.analysis``'s.

    Each argument is merged over the matching analysis table when analysis
    mode is on; without it only the given values are used.

    Raises
    ------
    ValueError
        If *parameters* holds ``nden`` (densities go in *densities*) or a
        density array has the wrong length.
    """
    if parameters is not None and NDEN in parameters:
        raise ValueError(
            "Densities are not parameters: pass them as densities= (by species "
            "name, or an array in species order), not as parameters['nden']"
        )
    analysis = net.analysis
    stored_p: dict[str, Any] = {}
    stored_d: dict[str, Any] = {}
    stored_f: dict[str, Any] = {}
    if analysis.enabled:
        tables = analysis.tables()
        stored_p, stored_d, stored_f = (t.values() for t in tables)
    merged_p = stored_p | dict(parameters or {})
    merged_d = stored_d | _densities(net, densities)
    merged_f = stored_f | dict(functions or {})
    index = {s.name: s.index for s in net.species}

    return Inputs(
        parameters={k: float(v) for k, v in merged_p.items()},
        densities={index[k]: float(v) for k, v in merged_d.items() if k in index},
        functions={k: _function(k, v) for k, v in merged_f.items()},
        names={i: n for n, i in index.items()},
    )


def _with_functions(expr: Expr, inputs: Inputs) -> Expr:
    """*expr* with every unknown function that has a value replaced by it."""
    return cast(
        Expr,
        expr.replace(
            lambda e: isinstance(e, AppliedUndef) and e.func.__name__ in inputs.functions,
            lambda e: inputs.functions[e.func.__name__],
        ),
    )


def _needs(expr: Expr) -> tuple[set[str], set[int], set[str]]:
    """Parameters, density indices and unknown functions *expr* needs."""
    functions = {f.func.__name__ for f in expr.atoms(AppliedUndef)}
    # A function without a value is reported; its arguments are never used.
    bare = expr.replace(lambda e: isinstance(e, AppliedUndef), lambda e: Float(1))
    indexed = bare.atoms(Indexed)
    labels = {a.base.label for a in indexed}
    densities = {int(a.indices[0]) for a in indexed if str(a.base.label) == NDEN}
    parameters = {s.name for s in bare.free_symbols - labels if isinstance(s, Symbol)} | {
        str(a.base.label) for a in indexed if str(a.base.label) != NDEN
    }

    return parameters, densities, functions


def check(exprs: Iterable[Expr], inputs: Inputs) -> list[Expr]:
    """*exprs* with the functions replaced; one error names every missing value.

    Raises
    ------
    ValueError
        Naming every parameter, species density and function without a value.
    """
    ready = [_with_functions(e, inputs) for e in exprs]
    parameters: set[str] = set()
    densities: set[int] = set()
    functions: set[str] = set()
    for expr in ready:
        p, d, f = _needs(expr)
        parameters |= p - set(inputs.parameters)
        densities |= d - set(inputs.densities)
        functions |= f

    parts = []
    if parameters:
        parts.append(f"no value for parameter {', '.join(sorted(parameters))}")
    if densities:
        names = sorted(
            (inputs.names.get(i, f"nden[{i}]") for i in densities), key=str.lower
        )
        parts.append(f"no density for {', '.join(names)}")
    if functions:
        parts.append(f"no value for function {', '.join(sorted(functions))}")
    if parts:
        raise ValueError(f"Cannot evaluate: {'; '.join(parts)}")

    return ready


def number(expr: Expr, inputs: Inputs) -> float:
    """The value of *expr*, already through :func:`check`."""
    mapping: dict[Basic, Float] = {
        a: Float(inputs.densities[int(a.indices[0])])
        for a in expr.atoms(Indexed)
        if str(a.base.label) == NDEN
    }
    mapping |= {
        s: Float(inputs.parameters[s.name])
        for s in expr.free_symbols
        if isinstance(s, Symbol) and s.name in inputs.parameters
    }

    return float(expr.xreplace(mapping).evalf())


def evaluate(
    net: Network,
    exprs: Sequence[Expr],
    parameters: Mapping[str, Any] | None = None,
    densities: Mapping[str, Any] | Sequence[float] | np.ndarray | None = None,
    functions: Mapping[str, FunctionValue] | None = None,
) -> np.ndarray:
    """Values of *exprs*, with the given values merged over ``net.analysis``'s.

    Raises
    ------
    ValueError
        Naming every parameter, species density and function without a value, or
        if an expression holds ``krate`` symbols.
    """
    if not exprs:
        return np.zeros(0)
    if any(str(a.base.label) == KRATE for e in exprs for a in sympify(e).atoms(Indexed)):
        raise ValueError(
            "Cannot evaluate krate symbols: evaluate the expanded form (rate.expanded, "
            "or substitute an OdeSystem's rate_coefficients)"
        )
    inputs = resolve(net, parameters, densities, functions)
    ready = check(exprs, inputs)

    return np.array([number(e, inputs) for e in ready], dtype=float)
