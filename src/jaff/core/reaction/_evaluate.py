# ABOUTME: Numeric evaluation of rate expressions: symbol values, nden arrays and
# ABOUTME: unknown (.jfunc) functions supplied by the caller; missing names reported

from __future__ import annotations

from typing import Any, Callable, Iterable, Mapping

from sympy import Expr, Float, Indexed, Symbol
from sympy.core.function import AppliedUndef

#: A value for an unknown function: called with its numeric arguments, or a constant.
FunctionValue = Callable[..., float] | float


def required_names(expr: Expr) -> tuple[set[str], set[str]]:
    """The symbol names and unknown-function names *expr* needs values for.

    An indexed density ``nden[3]`` needs ``nden`` (an array).
    """
    symbols = {s.name for s in expr.free_symbols if isinstance(s, Symbol)}
    functions = {f.func.__name__ for f in expr.atoms(AppliedUndef)}

    return symbols, functions


def check_missing(
    exprs: Iterable[Expr],
    values: Mapping[str, Any],
    functions: Mapping[str, FunctionValue],
) -> None:
    """Raise one error naming every symbol and function without a value.

    Raises
    ------
    ValueError
        If any expression needs a name missing from *values* or *functions*.
    """
    missing_symbols: set[str] = set()
    missing_functions: set[str] = set()
    for expr in exprs:
        symbols, funcs = required_names(expr)
        missing_symbols |= symbols - set(values)
        missing_functions |= funcs - set(functions)

    parts = []
    if missing_symbols:
        parts.append(f"no value for {', '.join(sorted(missing_symbols))}")
    if missing_functions:
        parts.append(f"no function value for {', '.join(sorted(missing_functions))}")
    if parts:
        raise ValueError(f"Cannot evaluate the rate coefficients: {'; '.join(parts)}")


def _call(function: AppliedUndef, functions: Mapping[str, FunctionValue]) -> Float:
    """The value of one unknown function whose arguments are already numbers."""
    given = functions[function.func.__name__]
    if not callable(given):
        return Float(given)

    return Float(given(*(float(Float(arg)) for arg in function.args)))


def evaluate_expression(
    expr: Expr,
    values: Mapping[str, Any],
    functions: Mapping[str, FunctionValue] | None = None,
) -> float:
    """Evaluate *expr* numerically.

    Parameters
    ----------
    expr : Expr
        The expression, e.g. a rate coefficient.
    values : mapping
        A number per symbol name (``"tgas"``, ``"av"`` ...); an indexed symbol
        such as ``nden`` takes an array, indexed as in the expression.
    functions : mapping, optional
        Per unknown function name, a callable taking the function's numeric
        arguments, or a constant. Nested calls are evaluated inside out.

    Returns
    -------
    float
        The value.

    Raises
    ------
    ValueError
        If a symbol or function has no value.
    """
    functions = functions or {}
    check_missing([expr], values, functions)

    mapping: dict[Expr, Float] = {}
    labels = set()
    for atom in expr.atoms(Indexed):
        index = tuple(int(i) for i in atom.indices)
        array = values[atom.base.label.name]
        mapping[atom] = Float(array[index[0]] if len(index) == 1 else array[index])
        labels.add(atom.base.label)
    for symbol in expr.free_symbols:
        if isinstance(symbol, Symbol) and symbol not in labels:
            mapping[symbol] = Float(values[symbol.name])

    numeric = expr.xreplace(mapping)
    numeric = numeric.replace(
        lambda e: isinstance(e, AppliedUndef), lambda e: _call(e, functions)
    )

    return float(numeric.evalf())
