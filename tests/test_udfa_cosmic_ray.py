# ABOUTME: Regression tests for RATE22 cosmic-ray rate formulas (CP direct, CR induced).
# ABOUTME: Verifies reference values, beta preservation, and vanishing at zero crate.

from pathlib import Path

import pytest
import sympy as sp

from jaff.core.parsers.network._formats.udfa.reaction import UdfaReaction

REPO = Path(__file__).resolve().parent.parent
BUNDLED = REPO / "networks" / "rate22_final" / "rate22_final.rates.jet"

ZETA0 = 1.36e-17  # reference cosmic-ray ionisation rate (s^-1)
OMEGA = 0.5  # grain albedo


def _bundled_line(index: str) -> str:
    with open(BUNDLED, encoding="utf-8") as f:
        for line in f:
            if line.strip().split(":", 1)[0] == index:
                return line
    raise AssertionError(f"reaction {index} not found")


def _segment(index: str) -> dict:
    (seg,) = UdfaReaction().parse(_bundled_line(index), 1, {}, BUNDLED)
    return seg


def _evaluate(rate: str, *, crate: float, tgas: float = 300.0, av: float = 0.0) -> float:
    expr = sp.sympify(rate)
    subs = {sp.Symbol("crate"): crate, sp.Symbol("tgas"): tgas, sp.Symbol("av"): av}
    return float(expr.subs(subs))


def test_cp_direct_ionization_reference() -> None:
    """CP (row 821): k = alpha * zeta/zeta0, i.e. equals alpha at zeta = zeta0."""
    seg = _segment("821")
    assert seg["type"] == "cosmic_ray"
    # alpha = 2.30e-17
    assert _evaluate(seg["rate"], crate=ZETA0) == pytest.approx(2.30e-17, rel=1e-3, abs=1e-25)


def test_cp_vanishes_at_zero_crate() -> None:
    """Direct CR ionization must be disableable via crate=0."""
    seg = _segment("821")
    assert _evaluate(seg["rate"], crate=0.0) == pytest.approx(0.0, abs=1e-30)


def test_cr_induced_photo_reference() -> None:
    """CR (row 833): k = alpha*gamma/(1-omega) * zeta/zeta0 = 6.5e-15 at zeta0."""
    seg = _segment("833")
    assert seg["type"] == "cosmic_ray"
    expected = 1.30e-17 * 250.0 / (1.0 - OMEGA)  # zeta/zeta0 = 1
    assert expected == pytest.approx(6.5e-15, rel=1e-6)
    assert _evaluate(seg["rate"], crate=ZETA0) == pytest.approx(6.5e-15, rel=1e-3, abs=1e-25)


def test_cr_vanishes_at_zero_crate() -> None:
    """CR-induced photoreaction must vanish at crate=0."""
    seg = _segment("833")
    assert _evaluate(seg["rate"], crate=0.0) == pytest.approx(0.0, abs=1e-30)


def test_cr_beta_preserved() -> None:
    """CR (row 990, beta=1.17): temperature dependence (T/300)**beta is applied."""
    seg = _segment("990")
    at_300 = _evaluate(seg["rate"], crate=ZETA0, tgas=300.0)
    at_600 = _evaluate(seg["rate"], crate=ZETA0, tgas=600.0)
    # Ratio must be 2**beta, not 1 (which is what a dropped beta would give).
    assert at_600 / at_300 == pytest.approx(2.0**1.17, rel=1e-3)
    # Reference at T=300: alpha*gamma/(1-omega).  Tolerance is 1% because the
    # emitted coefficient is formatted to two significant figures.
    assert at_300 == pytest.approx(1.30e-17 * 105.0 / (1.0 - OMEGA), rel=1e-2, abs=1e-25)
