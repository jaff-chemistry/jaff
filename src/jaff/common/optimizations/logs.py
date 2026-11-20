# ABOUTME: shared_logs(): non-integer powers of a base used with several exponents become
# ABOUTME: exp(f*L) * b**n with one shared L = log(b) and one exp per fractional part f
"""Shared logarithms replacing ``pow()`` calls of a common base."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
from typing import Any, cast

from sympy import (
    Basic,
    Dummy,
    Expr,
    Indexed,
    Pow,
    Symbol,
    exp,
    floor,
    log,
    sqrt,
)

from .optimize import bottom_up

DIGITS = 12
"""Fractional parts equal to this many decimals share one exponential (``1.3 - 1`` and
``2.3 - 2`` differ in the last bit)."""


def shared_logs(
    exprs: Sequence[Basic], definitions: Mapping[Any, Basic]
) -> tuple[list[Basic], dict[Any, Basic]]:
    """*exprs* and *definitions* with the powers of each base used with several
    non-integer exponents rewritten over one shared ``L = log(base)``.

    ``pow(T, a)`` costs about 2.4 times ``exp`` (glibc), and a network raises
    ``tgas`` to dozens of exponents (GOW: 60 in its rhs and Jacobian).  Each such
    base gets a new symbol ``L`` defined as ``log(base)``, and:

    - for a symbol base, ``b**a`` (numeric ``a``) is ``P * b**n`` with ``n = floor(a)``
      and ``P`` a new symbol defined as ``exp(f*L)`` for ``f = a - n``, one per
      distinct ``f`` (``sqrt(b)`` for ``f = 1/2``); ``b**n`` is then a product
      (:func:`~jaff.common.optimizations.integer_power`), so ``T**1.3``, ``T**2.3`` and
      ``T**-0.7`` cost one ``exp`` and a few multiplications;
    - otherwise ``b**a`` is ``exp(a*L)`` (a compound base's integer power would stay a
      ``pow()`` call).

    The new symbols are returned as definitions, so
    :func:`~jaff.common.optimizations.shared_cse` computes each once (and guards
    them when only branches read them); as symbols, SymPy cannot fold
    ``exp(a*log(b))`` back into a power.  The relative error is about ``|a·L|`` units
    in the last place, within the error of the expressions themselves.  Numeric
    bases are left to :func:`~jaff.common.optimizations.exp_power`.

    Parameters
    ----------
    exprs : Sequence[Basic]
        The expressions.
    definitions : Mapping
        Symbol to the expression giving its value (e.g. ``krate[r]``); may be empty.

    Returns
    -------
    exprs : list of Basic
        The rewritten expressions.
    definitions : dict
        *definitions* rewritten, plus the new logarithms and exponentials.

    Examples
    --------
    >>> import sympy as sp
    >>> from jaff.common.optimizations import shared_logs
    >>> t = sp.Symbol("tgas")
    >>> exprs, defs = shared_logs([t**1.3, t**2.3, t**-0.5], {})
    >>> exprs
    [_P0*tgas, _P0*tgas**2, _P1/tgas]
    >>> list(defs.values())
    [log(tgas), exp(0.3*_L0), sqrt(tgas)]
    """
    keys = list(definitions)
    values = [definitions[k] for k in keys]
    exponents = _exponents([*exprs, *values])
    shared = [b for b, es in exponents.items() if len(es) >= 2]
    if not shared:
        return list(exprs), dict(definitions)

    rewrite = _Rewrite(shared)
    new_exprs = [rewrite(e) for e in exprs]
    new_definitions = {k: rewrite(v) for k, v in zip(keys, values)}

    return new_exprs, {**new_definitions, **rewrite.definitions}


def _exponents(exprs: Sequence[Basic]) -> dict[Basic, set[Basic]]:
    """Each non-numeric base's distinct non-integer exponents."""
    exponents: dict[Basic, set[Basic]] = defaultdict(set)
    seen: set[Basic] = set()
    pending = list(exprs)
    while pending:  # each distinct subtree once
        node = pending.pop()
        if node in seen:
            continue
        seen.add(node)
        if isinstance(node, Pow) and not node.base.is_Number and not node.exp.is_Integer:
            exponents[node.base].add(node.exp)
        pending.extend(node.args)

    return exponents


class _Rewrite:
    """Rewrites powers of the shared bases, creating their definitions on first use."""

    def __init__(self, bases: Sequence[Basic]) -> None:
        self.bases = set(bases)
        self.definitions: dict[Basic, Basic] = {}
        self._logs: dict[Basic, Expr] = {}
        self._parts: dict[tuple[Basic, Any], Expr] = {}
        # One cache for every expression: shared subtrees are rewritten once.
        self._rewrite = bottom_up([self._rule])

    def __call__(self, expr: Basic) -> Basic:
        return self._rewrite(expr)

    def _rule(self, node: Basic) -> Basic:
        if isinstance(node, Pow) and node.base in self.bases and not node.exp.is_Integer:
            return self._power(node)

        return node

    def _power(self, power: Basic) -> Basic:
        base, a = (cast(Expr, arg) for arg in power.args)
        if not a.is_Number or not isinstance(base, (Symbol, Indexed)):
            return exp(a * self._log(base))

        n = int(floor(a))
        return self._part(base, a - n) * base**n

    def _log(self, base: Expr) -> Expr:
        """``L``, defined as ``log(base)``."""
        if base not in self._logs:
            symbol = Dummy(f"L{len(self._logs)}")
            self._logs[base] = symbol
            self.definitions[symbol] = log(base)

        return self._logs[base]

    def _part(self, base: Expr, f: Expr) -> Expr:
        """``P``, defined as ``base**f`` (``0 < f < 1``): ``sqrt`` or ``exp(f*L)``."""
        half = float(f) == 0.5
        key = (base, "half" if half else round(float(f), DIGITS))
        if key not in self._parts:
            symbol = Dummy(f"P{len(self._parts)}")
            self._parts[key] = symbol
            self.definitions[symbol] = sqrt(base) if half else exp(f * self._log(base))

        return self._parts[key]
