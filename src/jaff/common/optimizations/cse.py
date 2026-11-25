# ABOUTME: shared_cse(): one CSE over expressions and symbol definitions (e.g. krate[r]),
# ABOUTME: each definition its own temporary, ordered, pruned, no aliases or conditions
"""Common-subexpression elimination shared between expressions and the symbols they use."""

from __future__ import annotations

from collections.abc import Iterator, Mapping, Sequence
from typing import Any

from sympy import Basic, Expr, Indexed, Symbol, cse

from .guard import guard
from sympy.logic.boolalg import Boolean

Temporaries = list[tuple[Basic, Expr]]


def shared_cse(
    exprs: Sequence[Basic],
    definitions: Mapping[Any, Basic],
    symbols: Iterator[Symbol],
    optimizations: str | None = None,
) -> tuple[Temporaries, list[Expr]]:
    """One CSE over *exprs* and the *definitions* of symbols they use.

    *definitions* maps symbols appearing in *exprs* (e.g. the ``krate[r]`` of a
    reduced ODE system) to the expressions giving their values.  Their values and
    *exprs* go through one :func:`sympy.cse`, so sub-expressions are shared between
    them, and each definition becomes a temporary of its own.  The temporaries are
    then:

    - ordered so each follows the ones it reads (stable otherwise); a definition
      may use another definition's symbol;
    - stripped of aliases: a temporary whose value is a symbol, an indexed symbol
      or a number is substituted where it is used;
    - without conditions: a temporary whose value is boolean (``tgas > 700``,
      ``a & b``) is substituted where it is used too, so every temporary is a
      number; typed languages (Fortran) declare all temporaries alike, and a
      comparison is cheap to repeat;
    - pruned of temporaries nothing uses.

    With no definitions and no ``Piecewise`` this is :func:`sympy.cse` followed by the
    pruning.  The result is exact: substituting the temporaries back, last to first,
    into the reduced expressions gives the values of *exprs* with *definitions*
    substituted.

    Parameters
    ----------
    exprs : Sequence[Basic]
        The expressions to compute.
    definitions : Mapping
        Symbol to the expression giving its value; may be empty.
    symbols : Iterator[Symbol]
        Names for the temporaries, consumed in order (e.g.
        ``sympy.numbered_symbols("cse")``); after pruning, some names may be
        skipped.
    optimizations : str, optional
        Passed to :func:`sympy.cse` (e.g. ``"basic"``).

    Returns
    -------
    temporaries : list of (Symbol, Expr)
        The temporaries in evaluation order.
    reduced : list of Expr
        One expression per entry of *exprs*, reading the temporaries; no
        definition symbol is left in either.

    Examples
    --------
    >>> import sympy as sp
    >>> from jaff.common.optimizations import shared_cse
    >>> t, x = sp.symbols("tgas x")
    >>> k = sp.IndexedBase("krate")
    >>> temporaries, reduced = shared_cse(
    ...     [k[0] * x + k[1], k[0] * x**2],
    ...     {k[0]: sp.exp(-50 / t), k[1]: 2 * sp.exp(-50 / t)},
    ...     sp.numbered_symbols("cse"),
    ... )
    >>> temporaries
    [(cse0, exp(-50/tgas)), (cse2, 2*cse0)]
    >>> reduced
    [cse0*x + cse2, cse0*x**2]

    ``krate[0]``'s own temporary, ``cse1 = cse0``, was an alias and is gone.
    """
    keys = list(definitions)
    if not exprs and not keys:
        return [], []

    temporaries, reduced = cse(
        [*(definitions[k] for k in keys), *exprs],
        symbols=symbols,
        optimizations=optimizations,
    )
    own = {k: next(symbols) for k in keys}
    defined = [*temporaries, *zip((own[k] for k in keys), reduced[: len(keys)])]
    defined = [(s, e.xreplace(own)) for s, e in defined]
    outputs = [e.xreplace(own) for e in reduced[len(keys) :]]

    defined, outputs = _without_aliases(_ordered(defined), outputs)

    return prune(guard(defined, outputs), outputs), outputs


def prune(temporaries: Temporaries, exprs: Sequence[Basic]) -> Temporaries:
    """*temporaries* without those that *exprs* do not use, directly or through others.

    :func:`sympy.cse` with ``optimizations="basic"`` can leave temporaries that are
    only read by temporaries no longer read themselves.  The order is kept.

    Parameters
    ----------
    temporaries : list of (Symbol, Expr)
        Temporaries in evaluation order.
    exprs : Sequence[Basic]
        The expressions reading them.

    Returns
    -------
    list of (Symbol, Expr)
    """
    value = dict(temporaries)
    names = set(value)
    used: set[Basic] = set()
    pending = [s for e in exprs for s in e.free_symbols & names]
    while pending:
        symbol = pending.pop()
        if symbol in used:
            continue
        used.add(symbol)
        pending.extend(value[symbol].free_symbols & names)

    return [(s, e) for s, e in temporaries if s in used]


def _ordered(temporaries: Temporaries) -> Temporaries:
    """*temporaries* ordered so each follows the ones it reads, otherwise stably."""
    value = dict(temporaries)
    position = {s: n for n, (s, _) in enumerate(temporaries)}
    done: set[Basic] = set()
    ordered: Temporaries = []
    for root, _ in temporaries:
        stack = [(root, False)]
        while stack:
            symbol, ready = stack.pop()
            if symbol in done:
                continue
            if ready:
                done.add(symbol)
                ordered.append((symbol, value[symbol]))
                continue
            stack.append((symbol, True))
            reads = (value[symbol].free_symbols & value.keys()) - done
            # Pushed last-first, so the earliest dependency is placed first.
            stack.extend((d, False) for d in sorted(reads, key=lambda s: -position[s]))

    return ordered


def _without_aliases(
    temporaries: Temporaries, exprs: list[Expr]
) -> tuple[Temporaries, list[Expr]]:
    """Substitute every temporary whose value is a symbol, a number or a boolean."""
    alias: dict[Basic, Basic] = {}
    kept: Temporaries = []
    for symbol, value in temporaries:
        value = value.xreplace(alias)
        if value.is_Atom or isinstance(value, (Indexed, Boolean)):
            alias[symbol] = value
        else:
            kept.append((symbol, value))

    return kept, [e.xreplace(alias) for e in exprs]
