# ABOUTME: Tests for optimizations.shared_cse: one CSE over expressions and symbol
# ABOUTME: definitions, topologically ordered, pruned, without aliases, always exact

import pytest
import sympy as sp

from jaff.common.optimizations import shared_cse

T, X, Y = sp.symbols("tgas x y")
KRATE = sp.IndexedBase("krate")
K0, K1, K2 = KRATE[0], KRATE[1], KRATE[2]

DEFINITIONS = {
    K0: 2e-10 * T**0.5 * sp.exp(-50 / T),
    K1: 3e-11 * sp.exp(-50 / T),
    K2: sp.Float(1e-9),
}
EXPRS = [
    K0 * X * Y + K1 * X - K2 * Y,
    -K0 * X * Y + K2 * Y,
    K0 * Y + K1,
    K0 * X,
    K0 * Y * (-0.5 / T),
]


def _cse(exprs=EXPRS, definitions=DEFINITIONS, **kwargs):
    return shared_cse(exprs, definitions, sp.numbered_symbols("cse"), **kwargs)


def _expand(temporaries, expr):
    """*expr* with every temporary substituted back, last to first."""
    for symbol, value in reversed(temporaries):
        expr = expr.xreplace({symbol: value})

    return expr


def _value(expr):
    return float(expr.xreplace({T: 300.0, X: 1.7, Y: 0.4}))


def test_substituting_back_gives_the_inputs():
    temporaries, reduced = _cse()

    assert len(reduced) == len(EXPRS)
    for got, want in zip(reduced, EXPRS, strict=True):
        assert _value(_expand(temporaries, got)) == pytest.approx(
            _value(want.xreplace(DEFINITIONS)), rel=1e-12
        )


def test_every_temporary_is_defined_before_it_is_used():
    temporaries, reduced = _cse()
    names = {symbol for symbol, _ in temporaries}
    defined: set = set()

    for symbol, value in temporaries:
        assert value.free_symbols & names <= defined, symbol
        defined.add(symbol)


def test_no_definition_symbol_is_left():
    temporaries, reduced = _cse()

    assert not any(e.has(KRATE) for _, e in temporaries)
    assert not any(e.has(KRATE) for e in reduced)


def test_each_varying_definition_is_one_temporary():
    temporaries, reduced = _cse()
    values = {_value(_expand(temporaries, symbol)) for symbol, _ in temporaries}

    for k in (K0, K1):
        want = _value(DEFINITIONS[k])
        assert any(v == pytest.approx(want, rel=1e-12) for v in values), k


def test_every_temporary_is_used():
    temporaries, reduced = _cse()
    used = set().union(*(e.free_symbols for e in reduced))
    used |= set().union(*(e.free_symbols for _, e in temporaries))

    assert all(symbol in used for symbol, _ in temporaries)


def test_no_temporary_is_an_alias():
    # krate[2] is a number: it is substituted, not kept as `cseN = 1e-9`.
    temporaries, _ = _cse()

    assert not any(
        value.is_Atom or isinstance(value, sp.Indexed) for _, value in temporaries
    )


def test_identical_definitions_share_one_temporary():
    rate = 2e-10 * T**0.5 * sp.exp(-50 / T)
    temporaries, reduced = _cse([K0 * X + K1 * Y, K0 * Y], {K0: rate, K1: rate})

    assert reduced[0].free_symbols - {X, Y} == reduced[1].free_symbols - {Y}
    assert _value(_expand(temporaries, reduced[0])) == pytest.approx(
        _value(rate * (X + Y)), rel=1e-12
    )


def test_nested_definitions_are_ordered():
    # A definition reading another definition's symbol.
    definitions = {K1: 2 * K0 * sp.exp(-1 / T), K0: T**0.5 * sp.exp(-50 / T)}
    temporaries, reduced = _cse([K1 * X + K0 * Y], definitions)
    names = {symbol for symbol, _ in temporaries}
    defined: set = set()

    for symbol, value in temporaries:
        assert value.free_symbols & names <= defined
        defined.add(symbol)
    assert _value(_expand(temporaries, reduced[0])) == pytest.approx(
        _value((K1 * X + K0 * Y).xreplace(definitions).xreplace(definitions)),
        rel=1e-12,
    )


def test_temporaries_come_from_the_given_names():
    temporaries, _ = shared_cse(EXPRS, DEFINITIONS, sp.numbered_symbols("tmp", start=5))

    assert all(str(symbol).startswith("tmp") for symbol, _ in temporaries)
    assert len({symbol for symbol, _ in temporaries}) == len(temporaries)


def test_without_definitions_is_plain_cse():
    exprs = [sp.exp(-50 / T) * X + Y, sp.exp(-50 / T) * Y, (X + Y) ** 2 * T]
    temporaries, reduced = shared_cse(exprs, {}, sp.numbered_symbols("cse"))
    expected = sp.cse(exprs, symbols=sp.numbered_symbols("cse"))

    assert (temporaries, reduced) == (expected[0], expected[1])


def test_empty():
    assert shared_cse([], {}, sp.numbered_symbols("cse")) == ([], [])
