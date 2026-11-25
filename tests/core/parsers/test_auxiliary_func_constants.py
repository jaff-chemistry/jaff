# ABOUTME: Tests for built-in CGS constants in .jfunc files: names like kB and
# ABOUTME: c_light are replaced by their value; names clashing with them raise

from pathlib import Path

import pytest
import sympy as sp

from jaff.core.parsers.auxiliary_func import AuxiliaryFunctionParser
from jaff.core.parsers.auxiliary_func._engine import CONSTANTS
from jaff.errors import ParserError
from jaff.physics import constants as c

#: Each built-in constant and its CGS value from jaff.physics.constants.
EXPECTED = {
    "c_light": c.c.cgs.value,
    "h_planck": c.h.cgs.value,
    "hbar": c.hbar.cgs.value,
    "kB": c.k_B.cgs.value,
    "G_grav": c.G.cgs.value,
    "e_charge": c.e.esu.value,
    "m_e": c.m_e.cgs.value,
    "m_p": c.m_p.cgs.value,
    "m_n": c.m_n.cgs.value,
    "amu": c.u.cgs.value,
    "N_A": c.N_A.value,
    "sigma_SB": c.sigma_sb.cgs.value,
    "a_rad": c.a_rad.cgs.value,
    "sigma_T": c.sigma_T.cgs.value,
    "a_bohr": c.a0.cgs.value,
    "alpha_fs": c.alpha.value,
    "Ry": c.Ry.to("erg").value,
    "pc": c.pc.cgs.value,
    "au": c.au.cgs.value,
    "M_sun": c.M_sun.cgs.value,
    "L_sun": c.L_sun.cgs.value,
    "R_sun": c.R_sun.cgs.value,
}


def _parse(tmp_path: Path, text: str) -> dict:
    path = tmp_path / "network.jfunc"
    path.write_text(text)

    return AuxiliaryFunctionParser(path).get_dict()


def test_table_is_the_cgs_values():
    assert set(CONSTANTS) == set(EXPECTED)
    for name, value in EXPECTED.items():
        assert CONSTANTS[name] == pytest.approx(value, rel=1e-12), name


def test_constant_in_a_return(tmp_path):
    funcs = _parse(tmp_path, "@function e(t)\n@return kB * t\n")

    t = sp.Symbol("t")

    assert funcs["e"].expr.free_symbols == {t}
    assert float(funcs["e"].expr.coeff(t)) == pytest.approx(EXPECTED["kB"])


def test_constant_in_a_var_and_a_local(tmp_path):
    funcs = _parse(
        tmp_path,
        "@var a_c = c_light * a_rad\n@function e(t)\n    m = 28 * amu\n@return a_c * m\n",
    )

    expected = EXPECTED["c_light"] * EXPECTED["a_rad"] * 28 * EXPECTED["amu"]
    assert float(funcs["e"].expr) == pytest.approx(expected)


def test_names_are_case_sensitive(tmp_path):
    funcs = _parse(tmp_path, "@function e(t)\n@return KB * t\n")

    assert sp.Symbol("KB") in funcs["e"].expr.free_symbols


def test_constants_and_converters_together(tmp_path):
    funcs = _parse(tmp_path, "@function e()\n@return eV_to_erg(1) / kB\n")

    assert float(funcs["e"].expr) == pytest.approx(c.e.si.value * 1e7 / EXPECTED["kB"])


@pytest.mark.parametrize(
    "text",
    [
        "@var kB = 1.38e-16\n@function e(t)\n@return t\n",
        "@function e(t)\n    pc = 3.086e18\n@return t / pc\n",
        "@function e(c_light)\n@return c_light\n",
        "@function m_p(t)\n@return t\n",
    ],
    ids=["var", "local", "argument", "function"],
)
def test_name_clashing_with_a_constant_raises(tmp_path, text):
    with pytest.raises(ParserError, match="built-in constant"):
        _parse(tmp_path, text)
