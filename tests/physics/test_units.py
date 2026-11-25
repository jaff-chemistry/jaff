# ABOUTME: Tests for jaff.physics.units: every conversion function against astropy,
# ABOUTME: round trips, and numbers / NumPy arrays / SymPy expressions as inputs

import astropy.units as u
import numpy as np
import pytest
import sympy as sp

from jaff.physics import units

EQUIV = u.spectral() + u.temperature_energy()

ENERGY = {
    "eV": u.eV,
    "keV": u.keV,
    "erg": u.erg,
    "J": u.J,
    "K": u.K,
    "Hz": u.Hz,
    "cm1": 1 / u.cm,
}
LENGTH = {
    "m": u.m,
    "km": u.km,
    "nm": u.nm,
    "um": u.um,
    "AA": u.AA,
    "au": u.au,
    "pc": u.pc,
    "kpc": u.kpc,
    "Mpc": u.Mpc,
    "ly": u.lyr,
}
TIME = {
    "minute": u.min,
    "hour": u.h,
    "day": u.d,
    "year": u.yr,
    "kyr": u.kyr,
    "Myr": u.Myr,
    "Gyr": u.Gyr,
}
MASS = {"kg": u.kg, "amu": u.u, "Msun": u.M_sun}
VELOCITY = {"km_s": u.km / u.s, "pc_Myr": u.pc / u.Myr}

#: Multiplicative conversions, as (from, to) unit-name pairs; both directions exist.
PAIRS = [
    *[("eV", "erg"), ("eV", "K"), ("K", "erg"), ("Hz", "eV"), ("Hz", "erg")],
    *[("cm1", "eV"), ("Hz", "cm1"), ("keV", "eV"), ("J", "erg")],
    *[(name, "cm") for name in LENGTH],
    *[("m", "km"), ("nm", "AA"), ("nm", "um"), ("um", "AA")],
    *[("au", "pc"), ("pc", "kpc"), ("kpc", "Mpc"), ("pc", "ly")],
    *[(name, "s") for name in TIME],
    *[("day", "year"), ("year", "Myr"), ("kyr", "Myr"), ("Myr", "Gyr")],
    *[(name, "g") for name in MASS],
    ("km_s", "pc_Myr"),
]
#: Photon wavelength <-> energy / frequency (E = hc / lambda, nu = c / lambda).
PHOTON = [(w, q) for w in ("nm", "AA", "um") for q in ("eV", "Hz")]

UNIT = {**ENERGY, **LENGTH, **TIME, **MASS, **VELOCITY}
UNIT |= {"cm": u.cm, "s": u.s, "g": u.g}


def _both(pairs):
    return [p for a, b in pairs for p in ((a, b), (b, a))]


def _expected(value, a, b):
    return (value * UNIT[a]).to(UNIT[b], equivalencies=EQUIV).value


def test_every_conversion_is_public():
    names = {f"{a}_to_{b}" for a, b in _both(PAIRS + PHOTON)}

    assert set(units.__all__) == names
    assert len(units.__all__) == len(names)


@pytest.mark.parametrize(("a", "b"), _both(PAIRS + PHOTON))
def test_matches_astropy(a, b):
    convert = getattr(units, f"{a}_to_{b}")

    assert convert(2.5) == pytest.approx(_expected(2.5, a, b), rel=1e-12)


@pytest.mark.parametrize(("a", "b"), PAIRS + PHOTON)
def test_round_trip(a, b):
    there = getattr(units, f"{a}_to_{b}")
    back = getattr(units, f"{b}_to_{a}")

    assert back(there(3.7)) == pytest.approx(3.7, rel=1e-12)


def test_photon_conversion_is_inverse():
    # E = hc / lambda: half the wavelength, twice the energy.
    assert units.nm_to_eV(250.0) == pytest.approx(2 * units.nm_to_eV(500.0))
    assert units.nm_to_eV(500.0) == pytest.approx(2.479683968, rel=1e-9)


def test_julian_year():
    assert units.year_to_s(1.0) == pytest.approx(365.25 * 86400)


def test_numpy_arrays():
    values = np.array([1.0, 2.0, 4.0])

    np.testing.assert_allclose(units.eV_to_erg(values), values * units.eV_to_erg(1.0))
    np.testing.assert_allclose(units.nm_to_eV(values), units.nm_to_eV(1.0) / values)


def test_sympy_expressions():
    x = sp.Symbol("x")

    assert isinstance(units.eV_to_erg(x), sp.Expr)
    assert float(units.eV_to_erg(x).subs(x, 2)) == pytest.approx(units.eV_to_erg(2.0))
    assert float(units.nm_to_eV(x).subs(x, 500)) == pytest.approx(units.nm_to_eV(500.0))


@pytest.mark.parametrize("name", units.__all__)
def test_documented(name):
    doc = getattr(units, name).__doc__

    assert doc and "Parameters" in doc and "Returns" in doc


def test_exported_from_physics():
    from jaff import physics

    assert "units" in physics.__all__
    assert physics.units is units
