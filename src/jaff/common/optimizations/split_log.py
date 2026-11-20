# ABOUTME: split_log(): log(c*exp(b)*R) as log(c) + b + log(R), so an exp that underflows
# ABOUTME: at low temperature no longer makes a log(0); a pre-CSE rule (PRE_CSE_RULES)
"""Logarithms of products with exponential factors, split."""

from __future__ import annotations

import math
from typing import cast

from sympy import Add, Basic, Expr, Float, Mul, exp, log


def split_log(expr: Basic) -> Basic:
    """``log(c*exp(b1)*...*R)`` as ``log(c) + b1 + ... + log(R)``; else *expr*.

    Rates are often ``c*exp(-E/T)`` and appear inside logarithms: in a network's
    own formulas (GOW interpolates ``log`` of rates) or in derivatives (a power
    with a rate as its exponent gives ``log(k)``).  Below a few tens of kelvin
    ``exp(-E/T)`` underflows to 0 and the logarithm of the product is ``log(0)``:
    Python raises, C gives ``-inf`` and then ``NaN``.  Split, each exponential
    contributes its argument and nothing underflows.  The value is unchanged
    wherever the original was finite.

    Only exponential factors and a positive numeric coefficient (folded into a
    float) come out; the rest ``R`` stays in one logarithm, as SymPy does not
    know that, e.g., ``tgas`` is positive.  ``log(exp(b))`` is ``b``.  A product
    without an exponential, or with a non-positive coefficient (whose log would
    be complex), is left alone.  Run before CSE (:data:`PRE_CSE_RULES`): after it,
    the argument may be a temporary.

    Parameters
    ----------
    expr : Basic
        One node of an expression.

    Returns
    -------
    Basic

    Examples
    --------
    >>> import sympy as sp
    >>> from jaff.common.optimizations import split_log
    >>> t = sp.Symbol("tgas")
    >>> split_log(sp.log(3.52e-9 * sp.exp(-43900 / t)))
    -19.4648048473364 - 43900/tgas
    """
    if not isinstance(expr, log):
        return expr

    arg = expr.args[0]
    if isinstance(arg, exp):
        return arg.args[0]
    if not isinstance(arg, Mul):
        return expr

    exps = [a for a in arg.args if isinstance(a, exp)]
    coefficient = Mul(*(a for a in arg.args if a.is_Number))
    if not exps or not coefficient.is_positive:
        return expr

    rest = Mul(*(a for a in arg.args if not a.is_Number and not isinstance(a, exp)))
    terms = [cast(Expr, e.args[0]) for e in exps]
    if coefficient != 1:
        terms.insert(0, Float(math.log(float(coefficient))))
    if rest != 1:
        terms.append(log(rest))

    return Add(*terms)
