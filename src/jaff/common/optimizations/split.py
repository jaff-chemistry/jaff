# ABOUTME: split_log (log(c*exp(b)*R) as log(c) + b + log(R)) and split_power ((c*exp(b)*R)**e
# ABOUTME: as c**e*exp(e*b)*R**e): exponential factors (exp, 10**x) out of logs and powers
"""Exponential factors split out of logarithms and powers."""

from __future__ import annotations

import math
from typing import cast

from sympy import Add, Basic, Expr, Float, Mul, Pow, exp, log


def split_log(expr: Basic) -> Basic:
    """``log(c*E1*...*R)`` as ``log(c) + ln E1 + ... + log(R)``; else *expr*.

    The ``E`` are exponential factors: ``exp(b)`` (``ln E = b``) or ``a**x`` with
    ``a`` a positive number and ``x`` not a number (``ln E = x·ln a``).  Rates are
    often ``c*exp(-E/T)`` or ``10**(...)`` and appear inside logarithms: in a
    network's own formulas (GOW interpolates ``log`` of rates) or in derivatives (a
    power with a rate as its exponent gives ``log(k)``).  Below a few tens of
    kelvin the exponential underflows to 0 and the logarithm of the product is
    ``log(0)``: Python raises, C gives ``-inf`` and then ``NaN``.  Split, each
    exponential contributes its exponent and nothing underflows.  The value is
    unchanged wherever the original was finite.

    Only exponential factors and a positive numeric coefficient (folded into a
    float) come out; the rest ``R`` stays in one logarithm, as SymPy does not know
    that, e.g., ``tgas`` is positive.  A product without an exponential factor, or
    with a non-positive coefficient (whose log would be complex), is left alone.
    Run before CSE (:data:`~jaff.common.optimizations.PRE_CSE_RULES`): after it,
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

    parts = _parts(expr.args[0])
    if parts is None:
        return expr

    coefficient, exponent, rest = parts
    terms = [exponent]
    if coefficient != 1:
        terms.insert(0, Float(math.log(float(coefficient))))
    if rest != 1:
        terms.append(log(rest))

    return Add(*terms)


def split_power(expr: Basic) -> Basic:
    """``(c*E1*...*R)**e`` (``e`` not an integer) as ``c**e * exp(e*(ln E1 + ...)) *
    R**e``; else *expr*.

    The ``E`` are exponential factors as in :func:`split_log`.  A rate raised to a
    power that depends on the state, such as COthin's H2 dissociation
    ``(1.18e-10*exp(-69500/T))**(1/(1 + n/n_cr))``, differentiates as ``B**e·(e'·log B
    + e·B'/B)``: in ``B'/B`` the ``exp(-E/T)`` of ``B'`` and the ``1/exp(-E/T)`` of
    ``1/B`` are not cancelled, and the Jacobian holds an ``exp(+E/T)`` that overflows
    at low temperature.  Split before differentiating, the derivative holds
    ``e·(E/T**2)`` instead.  Valid because ``c > 0`` and every ``E > 0``; the value
    is unchanged.

    Only exponential factors and a positive numeric coefficient come out; the rest
    stays one power.  A base without an exponential factor, with a non-positive
    coefficient, or an integer exponent is left alone.

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
    >>> from jaff.common.optimizations import split_power
    >>> t, e = sp.symbols("tgas e")
    >>> split_power((1.18e-10 * sp.exp(-69500 / t)) ** e)
    1.18e-10**e*exp(-69500*e/tgas)
    """
    if not isinstance(expr, Pow) or expr.exp.is_Integer:
        return expr

    parts = _parts(expr.base)
    if parts is None:
        return expr

    coefficient, exponent, rest = parts
    e = expr.exp

    return coefficient**e * exp(e * exponent) * rest**e


def _parts(expr: Basic) -> tuple[Expr, Expr, Expr] | None:
    """``(c, Σ ln E, R)`` of the product *expr*, or ``None`` without an exponential
    factor or with a non-positive numeric coefficient ``c``."""
    factors = Mul.make_args(cast(Expr, expr))
    exponents = [_exponent(f) for f in factors]
    if all(x is None for x in exponents):
        return None

    coefficient = Mul(*(f for f in factors if f.is_Number))
    if not coefficient.is_positive:
        return None

    rest = Mul(*(f for f, x in zip(factors, exponents) if x is None and not f.is_Number))

    return coefficient, Add(*(x for x in exponents if x is not None)), rest


def _exponent(factor: Basic) -> Expr | None:
    """``ln factor`` if *factor* is ``exp(b)`` or ``a**x`` (``a`` a positive number,
    ``x`` not a number), else ``None``."""
    if isinstance(factor, exp):
        return cast(Expr, factor.args[0])

    if isinstance(factor, Pow):
        base, x = factor.args
        if base.is_Number and base.is_positive and base != 1 and not x.is_Number:
            return x * Float(math.log(float(base)))

    return None
