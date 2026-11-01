# ABOUTME: Tests for unit converters in .jfunc files: calls to jaff.physics.units
# ABOUTME: functions are replaced after parsing; names clashing with them raise

from pathlib import Path

import pytest
import sympy as sp

from jaff.core.parsers.auxiliary_func import AuxiliaryFunctionParser
from jaff.errors import ParserError
from jaff.physics import units


def _parse(tmp_path: Path, text: str) -> dict:
    path = tmp_path / "network.jfunc"
    path.write_text(text)

    return AuxiliaryFunctionParser(path).get_dict()


def test_converter_lower_case_names_are_unique():
    assert len({name.lower() for name in units.__all__}) == len(units.__all__)


def test_converter_call_is_replaced(tmp_path):
    funcs = _parse(tmp_path, "@function e(lam)\n@return nm_to_eV(lam)\n")
    lam = sp.Symbol("lam")

    assert funcs["e"].expr == units.nm_to_eV(lam)


def test_converter_on_a_number(tmp_path):
    funcs = _parse(tmp_path, "@function e()\n@return 2 * eV_to_erg(1.5)\n")

    assert float(funcs["e"].expr) == pytest.approx(2 * units.eV_to_erg(1.5))


def test_name_case_does_not_matter(tmp_path):
    funcs = _parse(tmp_path, "@function e(x)\n@return EV_TO_ERG(x)\n")

    assert funcs["e"].expr == units.eV_to_erg(sp.Symbol("x"))


def test_nested_converters(tmp_path):
    funcs = _parse(tmp_path, "@function e(x)\n@return nm_to_eV(um_to_nm(x))\n")
    x = sp.Symbol("x")

    assert funcs["e"].expr == units.nm_to_eV(units.um_to_nm(x))


def test_converter_in_a_var_and_a_local(tmp_path):
    funcs = _parse(
        tmp_path,
        "@var E0 = nm_to_eV(500)\n@function e(x)\n    t = eV_to_K(x)\n@return E0 + t\n",
    )
    x = sp.Symbol("x")

    assert sp.simplify(funcs["e"].expr - (units.nm_to_eV(500) + units.eV_to_K(x))) == 0


def test_converter_in_a_called_function(tmp_path):
    funcs = _parse(
        tmp_path,
        "@function lam_e(x)\n@return nm_to_eV(x)\n"
        "@function twice(x)\n@return 2 * lam_e(x)\n",
    )

    assert funcs["twice"].expr == 2 * units.nm_to_eV(sp.Symbol("x"))


def test_other_undefined_functions_stay(tmp_path):
    funcs = _parse(tmp_path, "@function k(t)\n@return rate_interp1d(t)\n")

    assert funcs["k"].expr == sp.Function("rate_interp1d")(sp.Symbol("t"))


def test_wrong_argument_count_raises(tmp_path):
    with pytest.raises(ParserError, match="nm_to_ev"):
        _parse(tmp_path, "@function e(x, y)\n@return nm_to_eV(x, y)\n")


@pytest.mark.parametrize(
    "text",
    [
        "@function nm_to_eV(x)\n@return x\n",
        "@function NM_TO_EV(x)\n@return x\n",
        "@var eV_to_erg = 1.6e-12\n@function e(x)\n@return x\n",
        "@function e(x)\n    km_to_cm = 1e5\n@return x * km_to_cm\n",
    ],
    ids=["function", "function-case", "var", "local"],
)
def test_name_clashing_with_a_converter_raises(tmp_path, text):
    with pytest.raises(ParserError, match="unit conversion"):
        _parse(tmp_path, text)
