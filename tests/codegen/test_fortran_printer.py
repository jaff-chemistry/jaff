# ABOUTME: Fortran code printer: user-defined function calls keep integer arguments
# ABOUTME: (e.g. the photorates index) while every other constant still prints as real

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
