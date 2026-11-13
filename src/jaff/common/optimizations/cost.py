# ABOUTME: cost(): weighted operation count of an expression for generated code, so
# ABOUTME: pow/exp/log and function calls weigh more than arithmetic when rewriting
"""How expensive an expression is to evaluate in generated code."""

from __future__ import annotations

from sympy import Add, Basic, Pow, count_ops

ARITHMETIC: dict[str, int] = {"ADD": 1, "SUB": 1, "MUL": 1, "NEG": 1, "DIV": 4, "POW": 2}
"""Weight of each cheap operation, by :func:`sympy.count_ops` name.

Addition, subtraction, multiplication and negation cost 1; a division 4; a power
with an integer exponent 2 (compilers turn it into multiplications)."""

BRANCH: dict[str, int] = {
    name: 1
    for name in ("PIECEWISE", "ITE", "LT", "LE", "GT", "GE", "EQ", "NE", "AND", "OR")
}
"""Weight of a branch or a comparison."""

CALL = 20
"""Weight of every other operation: ``exp``, ``log``, a power with a non-integer
exponent (``std::pow``, ``sqrt``) and any function call, e.g. an interpolation table."""


def cost(expr: Basic) -> int:
    """Weighted operation count of *expr*, a rough cost of evaluating it.

    :func:`sympy.count_ops` with each operation weighted: arithmetic is cheap
    (:data:`ARITHMETIC`), branches cost 1 (:data:`BRANCH`), and ``exp``, ``log``,
    non-integer powers and function calls cost :data:`CALL`.  Plain
    ``count_ops`` counts a ``std::pow`` like an addition, which makes a rewrite
    that removes one look no better.

    Parameters
    ----------
    expr : Basic
        The expression.

    Returns
    -------
    int

    Examples
    --------
    >>> import sympy as sp
    >>> from jaff.common.optimizations import cost
    >>> tgas, k = sp.symbols("tgas k")
    >>> cost(-3.3e-10 * tgas**-1.8), cost(-0.8 * k / tgas)
    (25, 6)
    """
    total = 0
    for term in Add.make_args(count_ops(expr, visual=True)):
        count, op = term.as_coeff_Mul()
        if op == 1:  # count_ops gives 0 for an atom
            continue

        name = str(op)
        total += int(count) * ARITHMETIC.get(name, BRANCH.get(name, CALL))

    # count_ops names every power POW; one with a non-integer exponent is a call.
    real_powers = sum(1 for p in expr.atoms(Pow) if not p.exp.is_Integer)

    return total + real_powers * (CALL - ARITHMETIC["POW"])
