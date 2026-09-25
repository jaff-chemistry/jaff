# ABOUTME: Regression tests for formal-argument binding in AuxiliaryFunctionParser.
# ABOUTME: Guards against sequential-substitution corruption of nested func calls.

from pathlib import Path

import pytest
import sympy as sp

from jaff.core.parsers.auxiliary_func import AuxiliaryFunctionParser
from jaff.errors import ParserError


def _write(tmp_path: Path, text: str) -> Path:
    path = tmp_path / "network.jfunc"
    path.write_text(text)
    return path


def test_swapped_arguments(tmp_path: Path) -> None:
    """A callee invoked with its own formals swapped must not collapse.

    ``swapped(x, y)`` calls ``subtract(y, x)``; binding must be simultaneous
    so the result is ``y - x`` rather than ``0``.
    """
    path = _write(
        tmp_path,
        "@function subtract(x,y)\n"
        "return x-y\n"
        "@function swapped(x,y)\n"
        "return subtract(y,x)\n",
    )
    funcs = AuxiliaryFunctionParser(path).get_dict()
    x, y = sp.symbols("x y")
    assert sp.simplify(funcs["swapped"]["def"] - (y - x)) == 0


def test_actual_arg_contains_other_formal(tmp_path: Path) -> None:
    """An actual-arg expression using another formal must bind exactly once."""
    path = _write(
        tmp_path,
        "@function subtract(x,y)\n"
        "return x-y\n"
        "@function combine(a,b)\n"
        "return subtract(a+b, b)\n",
    )
    funcs = AuxiliaryFunctionParser(path).get_dict()
    a, b = sp.symbols("a b")
    # subtract(a+b, b) = (a+b) - b = a
    assert sp.simplify(funcs["combine"]["def"] - a) == 0


def test_too_few_arguments(tmp_path: Path) -> None:
    """A call with fewer actuals than formals is rejected, not silently bound."""
    path = _write(
        tmp_path,
        "@function subtract(x,y)\n"
        "return x-y\n"
        "@function bad(t)\n"
        "return subtract(t)\n",
    )
    with pytest.raises(ParserError, match="subtract"):
        AuxiliaryFunctionParser(path)


def test_too_many_arguments(tmp_path: Path) -> None:
    """A call with more actuals than formals is rejected, not truncated."""
    path = _write(
        tmp_path,
        "@function subtract(x,y)\n"
        "return x-y\n"
        "@function bad(t)\n"
        "return subtract(t,1,99)\n",
    )
    with pytest.raises(ParserError, match="subtract"):
        AuxiliaryFunctionParser(path)


def test_zero_arg_func_called_with_arg(tmp_path: Path) -> None:
    """A zero-arg function called with an argument is rejected."""
    path = _write(
        tmp_path,
        "@function nullary()\n"
        "return 42\n"
        "@function bad(t)\n"
        "return nullary(t)\n",
    )
    with pytest.raises(ParserError, match="nullary"):
        AuxiliaryFunctionParser(path)
