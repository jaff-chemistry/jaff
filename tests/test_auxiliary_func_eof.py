# ABOUTME: Regression tests for end-of-file state validation in AuxiliaryFunctionParser.
# ABOUTME: Guards against open function blocks and pending continuations at EOF.

from pathlib import Path

import pytest
import sympy as sp

from jaff.core.parsers.auxiliary_func import AuxiliaryFunctionParser
from jaff.errors import ParserError


def _write(tmp_path: Path, text: str) -> Path:
    path = tmp_path / "network.jfunc"
    path.write_text(text)
    return path


def test_function_block_missing_return(tmp_path: Path) -> None:
    """A function block left open at EOF (no return) is rejected."""
    path = _write(
        tmp_path,
        "@function foo(x)\n"
        "y = 2*x\n",
    )
    with pytest.raises(ParserError, match="foo"):
        AuxiliaryFunctionParser(path)


def test_pending_continuation_at_eof(tmp_path: Path) -> None:
    """A directive left mid-continuation at EOF is rejected, not dropped."""
    path = _write(
        tmp_path,
        "@var a = 1 + \\\n",
    )
    with pytest.raises(ParserError):
        AuxiliaryFunctionParser(path)


def test_valid_return_no_trailing_newline(tmp_path: Path) -> None:
    """A complete function whose final return lacks a trailing newline parses."""
    path = _write(tmp_path, "@function foo(x)\nreturn 2*x")  # no newline at EOF
    funcs = AuxiliaryFunctionParser(path).get_dict()
    x = sp.Symbol("x")
    assert sp.simplify(funcs["foo"]["def"] - 2 * x) == 0
