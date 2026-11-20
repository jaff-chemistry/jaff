# ABOUTME: Tests for optimizations.guard: a CSE temporary read only inside Piecewise
# ABOUTME: branches is computed only under their (input-only, simplified) conditions

import pytest
import sympy as sp

from jaff.common.optimizations import guard

T, X = sp.symbols("tgas x")
C0, C1, C2 = sp.symbols("cse0 cse1 cse2")
ZERO = sp.Float(0.0)


def _condition(temporaries, symbol):
    """The guard of *symbol*'s value, ``True`` if it has none."""
    value = dict(temporaries)[symbol]
    if isinstance(value, sp.Piecewise) and value.args[-1] == (ZERO, True):
        return value.args[0].cond

    return sp.true


def _expand(temporaries, expr):
    """*expr* with every temporary substituted back, last to first."""
    for symbol, value in reversed(temporaries):
        expr = expr.xreplace({symbol: value})

    return expr


def test_unconditional_temporary_is_unchanged():
    temporaries = [(C0, sp.exp(X))]
    exprs = [C0 + sp.Piecewise((C0, T < 500), (1, True))]

    assert guard(temporaries, exprs) == temporaries


def test_branch_only_temporary_gets_its_branch_condition():
    temporaries = [(C0, sp.exp(X))]
    guarded = guard(temporaries, [sp.Piecewise((C0, T < 500), (1, True))])

    assert guarded == [(C0, sp.Piecewise((sp.exp(X), T < 500), (ZERO, True)))]


def test_later_branch_needs_the_earlier_conditions_false():
    exprs = [sp.Piecewise((1, T < 10), (C0, T < 500), (2, True))]

    assert _condition(guard([(C0, sp.exp(X))], exprs), C0) == sp.And(T >= 10, T < 500)


def test_default_branch_needs_every_condition_false():
    exprs = [sp.Piecewise((1, T < 10), (C0, True))]

    assert _condition(guard([(C0, sp.exp(X))], exprs), C0) == (T >= 10)


def test_several_branches_give_the_union_of_their_conditions():
    exprs = [
        sp.Piecewise((C0, T < 100), (1, True)),
        sp.Piecewise((C0, T > 1000), (2, True)),
    ]

    assert _condition(guard([(C0, sp.exp(X))], exprs), C0) == sp.Or(T < 100, T > 1000)


def test_nested_piecewise_conditions_combine():
    inner = sp.Piecewise((C0, T > 50), (1, True))
    exprs = [sp.Piecewise((inner, T < 500), (2, True))]

    assert _condition(guard([(C0, sp.exp(X))], exprs), C0) == sp.And(T > 50, T < 500)


def test_redundant_bounds_are_simplified():
    exprs = [sp.Piecewise((C0, (T < 100) & (T < 6000)), (1, True))]

    assert _condition(guard([(C0, sp.exp(X))], exprs), C0) == (T < 100)


def test_ite_condition_is_simplified():
    exprs = [sp.Piecewise((C0, sp.ITE(T < 100, T < 10, T > 500)), (1, True))]

    assert _condition(guard([(C0, sp.exp(X))], exprs), C0) == sp.Or(T < 10, T > 500)


def test_condition_on_several_variables_is_kept_as_is():
    exprs = [sp.Piecewise((C0, (T < 500) & (X > 1)), (1, True))]

    assert _condition(guard([(C0, sp.exp(X))], exprs), C0) == sp.And(T < 500, X > 1)


def test_temporaries_read_by_a_guarded_one_inherit_its_condition():
    temporaries = [(C0, sp.exp(X)), (C1, 2 * C0)]
    guarded = guard(temporaries, [sp.Piecewise((C1, T < 500), (1, True))])

    assert _condition(guarded, C0) == (T < 500)
    assert _condition(guarded, C1) == (T < 500)


def test_ite_branches_are_conditional():
    exprs = [sp.Piecewise((1, sp.ITE(T < 100, C0 < 1, C1 < 1)), (2, True))]
    guarded = guard([(C0, sp.exp(X)), (C1, sp.exp(2 * X))], exprs)

    assert _condition(guarded, C0) == (T < 100)
    assert _condition(guarded, C1) == (T >= 100)


def test_and_or_operands_after_the_first_printed_are_conditional():
    # Printed code short-circuits: `(tgas < 10) or (cse0*x < 1)` reads cse0 only
    # when tgas >= 10.
    condition = sp.Or(T < 10, C0 * X < 1)
    exprs = [sp.Piecewise((1, condition), (2, True))]

    assert sp.pycode(condition) == "(tgas < 10) or (cse0*x < 1)"
    assert _condition(guard([(C0, sp.exp(X))], exprs), C0) == (T >= 10)


def test_first_printed_operand_is_unconditional():
    # Printers sort the operands: here cse0 comes first, so it is always read.
    condition = sp.Or(T < 10, C0 < 1)
    temporaries = [(C0, sp.exp(X))]

    assert sp.pycode(condition) == "(cse0 < 1) or (tgas < 10)"
    assert guard(temporaries, [sp.Piecewise((1, condition), (2, True))]) == temporaries


def test_conditions_reading_temporaries_are_weakened():
    # The fit is tested in the condition selecting its own branch: guarding it by
    # that condition would read it before it is computed.  Only `tgas >= 10` stays.
    fit = C0
    exprs = [sp.Piecewise((1, sp.Or(T < 10, fit + X < 1e-100)), (fit * X, True))]
    guarded = guard([(fit, 10 ** (70 * sp.log(T) ** 8))], exprs)

    assert _condition(guarded, fit) == (T >= 10)


def test_guards_read_no_temporary():
    temporaries = [(C0, sp.exp(X)), (C1, C0 + 1), (C2, sp.exp(-X))]
    exprs = [
        sp.Piecewise((C2, C1 > 2), (C1, T < 300), (0, True)),
        sp.Piecewise((C2 * X, sp.And(T > 10, C0 < 5)), (1, True)),
    ]
    names = {symbol for symbol, _ in temporaries}

    for symbol, _ in temporaries:
        assert not _condition(guard(temporaries, exprs), symbol).free_symbols & names


@pytest.mark.parametrize("t", [5.0, 50.0, 250.0, 400.0, 2000.0])
@pytest.mark.parametrize("x", [-1.0, 0.5, 3.0])
def test_outputs_are_unchanged(t, x):
    temporaries = [(C0, sp.exp(X)), (C1, C0 + 1), (C2, sp.exp(-X) * C0)]
    exprs = [
        C2 + sp.Piecewise((C1, T < 300), (C0 * 2, T < 1000), (0, True)),
        sp.Piecewise((C2 * X, sp.And(T > 10, C1 < 5)), (1, True)),
        sp.Piecewise((1, sp.Or(T < 100, C1 > 2)), (C1 * X, True)),
    ]
    guarded = guard(temporaries, exprs)
    point = {T: t, X: x}

    for expr in exprs:
        want = float(_expand(temporaries, expr).xreplace(point))
        assert float(_expand(guarded, expr).xreplace(point)) == pytest.approx(want)


def test_no_piecewise_changes_nothing():
    temporaries = [(C0, sp.exp(X)), (C1, C0 + T)]

    assert guard(temporaries, [C1 * 2, C0]) == temporaries
