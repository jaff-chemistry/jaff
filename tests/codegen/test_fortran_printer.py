# ABOUTME: Fortran code printer: user-defined function calls keep integer arguments,
# ABOUTME: every other constant (Piecewise branch values too) prints as real

import sympy as sp
from sympy.core.function import UndefinedFunction

from jaff.codegen._languages import Language

x = sp.Symbol("x")
foo: UndefinedFunction = sp.Function("foo")  # type: ignore


def _fcode(expr) -> str:
    return Language("fortran").code_gen(expr, strict=False, allow_unknown_functions=True)


def test_user_function_keeps_integer_argument():
    assert _fcode(foo(1, sp.Float(2.5))) == "foo(1, 2.5d0)"


def test_user_function_evaluates_other_arguments_as_reals():
    assert _fcode(foo(sp.Rational(1, 2), x + 1)) == "foo(0.5d0, x + 1.0d0)"


def test_intrinsic_argument_still_printed_as_real():
    assert _fcode(sp.sin(x + 1)) == "sin(x + 1.0d0)"


def test_piecewise_integer_branch_prints_as_real():
    # merge() needs both values of the same type: 0 against a real is an error
    expr = sp.Piecewise((x**2, x > 700), (0, True))
    assert _fcode(expr) == "merge(x**2, 0.0d0, x > 700)"


def test_nested_piecewise_integer_branches_print_as_reals():
    expr = sp.Piecewise((1, x < 1), (x, x < 10), (2, True))
    assert _fcode(expr) == "merge(1.0d0, merge(x, 2.0d0, x < 10), x < 1)"
