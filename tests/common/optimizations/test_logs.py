# ABOUTME: Tests for optimizations.shared_logs: non-integer powers of a base used with
# ABOUTME: several exponents become exp(f*L) * b**n with one shared L = log(b)

import math

import pytest
import sympy as sp

from jaff.common.optimizations import shared_cse, shared_logs

T, X, Y = sp.symbols("tgas x y")
KRATE = sp.IndexedBase("krate")


def _logs(definitions):
    """The logarithms shared_logs added: symbol -> log(base)."""
    return {k: v for k, v in definitions.items() if isinstance(v, sp.log)}


def _value(expr, definitions, t=350.0, x=1.7, y=0.4):
    for _ in range(3):
        expr = expr.xreplace(definitions)
    return float(expr.xreplace({T: t, X: x, Y: y}))


def _real_powers(expr, base):
    return [
        p
        for p in sp.preorder_traversal(expr)
        if isinstance(p, sp.Pow) and p.base == base and not p.exp.is_Integer
    ]


def test_base_with_several_exponents_shares_one_log():
    exprs = [T**1.3 * X, T**-0.7 + Y, T**0.42]
    reduced, definitions = shared_logs(exprs, {})

    assert list(_logs(definitions).values()) == [sp.log(T)]
    assert not any(_real_powers(e, T) for e in reduced)
    for got, want in zip(reduced, exprs, strict=True):
        assert _value(got, definitions) == pytest.approx(_value(want, {}), rel=1e-14)


def test_base_with_one_exponent_is_unchanged():
    exprs = [T**1.3 * X, T**1.3 + T**2]

    assert shared_logs(exprs, {}) == (exprs, {})


def test_integer_powers_do_not_count():
    exprs = [T**1.3 * X, T**2, T**-3]

    assert shared_logs(exprs, {}) == (exprs, {})


def test_exponents_with_the_same_fractional_part_share_one_exponential():
    exprs = [T**1.3, T**2.3, T**-0.7]
    reduced, definitions = shared_logs(exprs, {})
    exps = [v for v in definitions.values() if isinstance(v, sp.exp)]

    assert len(exps) == 1
    assert not any(
        isinstance(e, sp.exp) for r in reduced for e in sp.preorder_traversal(r)
    )
    for got, want in zip(reduced, exprs, strict=True):
        assert _value(got, definitions) == pytest.approx(_value(want, {}), rel=1e-14)


def test_half_fractional_part_uses_sqrt():
    exprs = [T**1.5, T**-0.5, T**0.3]
    reduced, definitions = shared_logs(exprs, {})

    assert sp.sqrt(T) in definitions.values()
    for got, want in zip(reduced, exprs, strict=True):
        assert _value(got, definitions) == pytest.approx(_value(want, {}), rel=1e-14)


def test_compound_base_is_a_single_exponential():
    base = T + 1
    exprs = [base**1.3, base**2.7]
    reduced, definitions = shared_logs(exprs, {})

    assert list(_logs(definitions).values()) == [sp.log(base)]
    assert not any(
        isinstance(p, sp.Pow) and p.base == base
        for r in reduced
        for p in sp.preorder_traversal(r)
    )
    for got, want in zip(reduced, exprs, strict=True):
        assert _value(got, definitions) == pytest.approx(_value(want, {}), rel=1e-14)


def test_symbolic_exponent_is_an_exponential():
    exprs = [T**X, T**1.3]
    reduced, definitions = shared_logs(exprs, {})
    (log,) = _logs(definitions)

    assert reduced[0] == sp.exp(X * log)


def test_definition_values_are_rewritten_too():
    definitions = {KRATE[0]: 2e-10 * T**0.7, KRATE[1]: 3e-11 * T**-1.2}
    reduced, new = shared_logs([KRATE[0] * X + KRATE[1]], definitions)

    assert not _real_powers(new[KRATE[0]], T) and not _real_powers(new[KRATE[1]], T)
    assert _value(reduced[0], new) == pytest.approx(
        _value(reduced[0], definitions), rel=1e-14
    )


def test_shared_cse_computes_the_log_once():
    exprs = [T**1.3 * X, T**-0.7 * Y, T**0.42 + X]
    reduced, definitions = shared_logs(exprs, {})
    temporaries, out = shared_cse(reduced, definitions, sp.numbered_symbols("cse"))
    logs = [v for _, v in temporaries if v == sp.log(T)]

    assert len(logs) == 1
    for got, want in zip(out, exprs, strict=True):
        for symbol, value in reversed(temporaries):
            got = got.xreplace({symbol: value})
        assert float(got.xreplace({T: 350.0, X: 1.7, Y: 0.4})) == pytest.approx(
            float(want.xreplace({T: 350.0, X: 1.7, Y: 0.4})), rel=1e-14
        )


def test_value_at_extreme_temperatures():
    exprs = [T**13.536556, T**-2.5, T**0.042 * sp.exp(-2e-5 * T)]
    reduced, definitions = shared_logs(exprs, {})

    for t in (10.0, 1e5):
        for got, want in zip(reduced, exprs, strict=True):
            assert _value(got, definitions, t=t) == pytest.approx(
                _value(want, {}, t=t), rel=1e-13
            )
    assert math.isfinite(_value(reduced[0], definitions, t=1e5))
