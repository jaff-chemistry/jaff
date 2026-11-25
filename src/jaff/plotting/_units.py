# ABOUTME: Unit conversion (astropy) and mathtext axis labels for cross-section plots
# ABOUTME: Converts photon energy / wavelength and cross-section units from eV and cm²

"""Unit conversion and axis-label helpers for cross-section plots.

Photo cross sections are stored internally as photon energy in eV and cross
section in cm².  The plotting layer lets callers request other units; the pure
helpers here perform the conversion (via :mod:`astropy.units`) and render the
matching axis labels with mathtext.
"""

from __future__ import annotations

import astropy.units as u
import numpy as np

# Mapping from JAFF unit strings to astropy units.
_ENERGY_UNITS: dict[str, u.UnitBase] = {
    "eV": u.eV,
    "erg": u.erg,
    "nm": u.nm,
    "um": u.um,
}
_XSEC_UNITS: dict[str, u.UnitBase] = {
    "cm2": u.cm**2,
    "cm^2": u.cm**2,
    "Mb": u.Mbarn,
    "barn": u.barn,
}

# Mathtext rendering of unit strings for axis labels.
_UNIT_TEX: dict[str, str] = {
    "cm^2": r"cm$^2$",
    "cm2": r"cm$^2$",
    "Mb": "Mb",
    "barn": "barn",
    "eV": "eV",
    "erg": "erg",
    "nm": "nm",
    "um": r"$\mu$m",
}

_ENERGY_LABELS: dict[str, str] = {
    "eV": "Photon energy (eV)",
    "erg": "Photon energy (erg)",
    "nm": "Wavelength (nm)",
    "um": r"Wavelength ($\mu$m)",
}

#: Display labels for the photo cross-section processes.
PROCESS_LABELS: dict[str, str] = {
    "photo_absorption": "Photoabsorption",
    "photodecay": "Photodecay",
}


def convert_energy(
    value: float | np.ndarray, from_unit: str, to_unit: str
) -> float | np.ndarray:
    """Convert between photon energies and wavelengths via astropy.

    Energy <-> wavelength conversions use the :func:`astropy.units.spectral`
    equivalency.

    Parameters
    ----------
    value : float or numpy.ndarray
        Value(s) in *from_unit*.
    from_unit, to_unit : str
        Source and target units; one of ``"eV"``, ``"erg"``, ``"nm"``,
        ``"um"``.

    Returns
    -------
    float or numpy.ndarray
        Value(s) expressed in *to_unit*.

    Raises
    ------
    ValueError
        If either unit is not a recognised energy/wavelength unit.
    """
    if from_unit not in _ENERGY_UNITS:
        raise ValueError(f"Unknown energy unit: {from_unit}")
    if to_unit not in _ENERGY_UNITS:
        raise ValueError(f"Unknown energy unit: {to_unit}")

    q = np.asarray(value) * _ENERGY_UNITS[from_unit]
    return q.to(_ENERGY_UNITS[to_unit], equivalencies=u.spectral()).value


def convert_xsec(
    value: float | np.ndarray, from_unit: str, to_unit: str
) -> float | np.ndarray:
    """Convert a cross section between area units via astropy.

    Parameters
    ----------
    value : float or numpy.ndarray
        Value(s) in *from_unit*.
    from_unit, to_unit : str
        Source and target units; one of ``"cm^2"``/``"cm2"``, ``"Mb"``,
        ``"barn"``.

    Returns
    -------
    float or numpy.ndarray
        Value(s) expressed in *to_unit*.

    Raises
    ------
    ValueError
        If either unit is not a recognised cross-section unit.
    """
    if from_unit not in _XSEC_UNITS:
        raise ValueError(f"Unknown cross-section unit: {from_unit}")
    if to_unit not in _XSEC_UNITS:
        raise ValueError(f"Unknown cross-section unit: {to_unit}")

    q = np.asarray(value) * _XSEC_UNITS[from_unit]
    return q.to(_XSEC_UNITS[to_unit]).value


def fmt_unit(unit: str) -> str:
    """Render a unit string with mathtext superscripts where known."""
    return _UNIT_TEX.get(unit, unit)


def energy_label(unit: str) -> str:
    """Axis label for a photon energy/wavelength *unit*."""
    return _ENERGY_LABELS.get(unit, f"Photon energy ({fmt_unit(unit)})")


def xsec_label(unit: str) -> str:
    """Axis label for a cross-section *unit*."""
    return rf"Cross section $\sigma$ ({fmt_unit(unit)})"
