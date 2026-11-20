# ABOUTME: optimize(): the per-node rules (RULES) in one bottom-up pass right before
# ABOUTME: printing; PRE_CSE_RULES for before CSE; bottom_up(): the memoized pass itself
"""The per-node rewrites applied to an expression right before it is printed."""

from __future__ import annotations

from collections.abc import Callable, Sequence

from sympy import Basic

from .powers import exp_power, integer_power
from .split_log import split_log

Rule = Callable[[Basic], Basic]
"""A rewrite of one node: the node itself when it does not apply."""

RULES: tuple[Rule, ...] = (exp_power, integer_power)
"""The rules :func:`optimize` applies at each node, in order, right before printing.

``exp_power`` comes first: ``10**x`` becomes ``exp(...)``, no longer a power for
``integer_power``."""

PRE_CSE_RULES: tuple[Rule, ...] = (split_log,)
"""The rules ``Codegen`` applies at each node before CSE (with :func:`bottom_up`), for
rewrites that must see inside sub-expressions CSE would turn into temporaries."""


def optimize(expr: Basic, rules: Sequence[Rule] = RULES) -> Basic:
    """*expr* with each of *rules* applied at every node, in one bottom-up pass.

    A node's arguments are rewritten first, then the node is rebuilt from them and
    each rule applied to it in turn.  Each distinct subtree is rewritten once, however
    often it appears.  The rules (:data:`RULES`) only change how the expression is
    written, never its value; some return unevaluated forms that a later SymPy
    rewrite would undo, so optimize last, right before printing (``Codegen`` does).

    Parameters
    ----------
    expr : Basic
        The expression.
    rules : Sequence of callables, optional
        Node rewrites, each returning its node unchanged when it does not apply.

    Returns
    -------
    Basic

    Examples
    --------
    >>> import sympy as sp
    >>> from jaff.common.optimizations import optimize
    >>> x = sp.Symbol("x")
    >>> sp.cxxcode(optimize(3 * x**3 + 10 ** (x / 2)))
    'std::exp(1.151292546497023*x) + 3*((x*x)*x)'
    """
    return bottom_up(rules)(expr)


def bottom_up(rules: Sequence[Rule]) -> Callable[[Basic], Basic]:
    """A rewriter applying *rules* at every node, bottom-up, as :func:`optimize`.

    Its cache of rewritten subtrees is kept between calls, so rewriting several
    expressions that share subtrees visits each distinct subtree once.
    """
    done: dict[Basic, Basic] = {}

    def visit(node: Basic) -> Basic:
        if node in done:
            return done[node]

        rewritten = node
        if node.args:
            args = tuple(visit(a) for a in node.args)
            if args != node.args:
                rewritten = node.func(*args)
        for rule in rules:
            rewritten = rule(rewritten)

        done[node] = rewritten
        return rewritten

    return visit
