# ABOUTME: KROME BASE_GLOBALS: the n_X grammar, and KROME's built-in temperature shortcuts
# ABOUTME: (get_Tshortcut) available in rates without @var, overridable by @var

import math

import pytest
import sympy as sp

from jaff import Network
from jaff.core.parsers.network._formats.krome.parser import KromeParser
from tests._paths import REPO_ROOT

TGAS = 350.0
T0 = 288.0
TE = TGAS * 8.617343e-5
T32 = TGAS / 300.0
# KROME's get_Tshortcut (kromelib.py), evaluated at TGAS.
SHORTCUTS = {
    "logT": math.log10(TGAS),
    "lnT": math.log(TGAS),
    "Te": TE,
    "lnTe": math.log(TE),
    "T32": T32,
    "t3": T32,
    "t4": TGAS * 1e-4,
    "invT": 1 / TGAS,
    "invT32": 1 / T32,
    "invTgas": 1 / TGAS,
    "invTe": 1 / TE,
    "sqrTgas": math.sqrt(TGAS),
    "invsqrT32": 1 / math.sqrt(T32),
    "invsqrT": 1 / math.sqrt(TGAS),
    "sqrT32": math.sqrt(T32),
    "Tgas2": TGAS**2,
    "Tgas3": TGAS**3,
    "Tgas4": TGAS**4,
    "T0": T0,
    "T02": T0**2,
    "T03": T0**3,
    "T04": T0**4,
    "T0inv": 1 / T0,
}


def test_base_globals_use_n_x_grammar():
    bg = KromeParser.BASE_GLOBALS
    assert bg["get_hnuclei(n)"] == "n_H_nuc"
    assert bg["n(idx_h)"] == "n_H"
    assert bg["n(idx_h2)"] == "n_H2"
    assert bg["n_global(idx_h2)"] == "n_H2"
    # temperature/user shorthands unchanged
    assert bg["t32"] == "tgas/3e2"
    assert bg["user_av"] == "av"


def _rate(tmp_path, body):
    path = tmp_path / "net.dat"
    path.write_text("@format:idx,R,R,P,rate\n" + body + "\n", encoding="utf-8")
    (reaction,) = list(Network(str(path)).reactions)

    return reaction.rate.coefficient


@pytest.mark.parametrize(("name", "value"), SHORTCUTS.items())
def test_temperature_shortcut_is_defined(tmp_path, name, value):
    rate = _rate(tmp_path, f"1,H,H,H2,1d-10*{name}")
    tgas = sp.Symbol("tgas")

    assert rate.free_symbols <= {tgas}
    assert float(rate.subs(tgas, TGAS)) == pytest.approx(1e-10 * value, rel=1e-12)


def test_var_overrides_a_shortcut(tmp_path):
    rate = _rate(tmp_path, "@var:invT=2d0\n1,H,H,H2,1d-10*invT")

    assert float(rate) == pytest.approx(2e-10)


def test_bundled_popsicle_rates_use_no_undefined_shortcut():
    net = Network(
        str(REPO_ROOT / "networks" / "popsicle_semenov" / "react_popsicle_semenov.jet")
    )
    names = {
        str(s)
        for r in net.reactions
        if not isinstance(r.rate.coefficient, str)
        for s in r.rate.coefficient.free_symbols
    }

    assert not names & {"invt", "lnte", "invsqrt"}
