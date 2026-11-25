# ABOUTME: Unit conversion functions (energy, temperature, frequency, wavelength,
# ABOUTME: length, time, mass, velocity), factors from astropy; CGS-centred
"""Unit conversion functions.

Each function converts a value from one unit to another: ``eV_to_erg``,
``nm_to_eV``, ``pc_to_cm``, ``Myr_to_s`` ... The factors come from
:mod:`astropy.units` (CODATA constants, the same source as
:mod:`jaff.physics.constants`) and are computed once, on import.

A value may be a number, a NumPy array or a SymPy expression; the result is of
the same kind. Most conversions multiply by a factor. A photon's wavelength is
inversely related to its energy (E = h c / lambda) and frequency
(nu = c / lambda), so ``nm_to_eV``, ``eV_to_nm`` and the other photon
conversions divide a constant by the value instead.

Lengths, times and masses convert to and from CGS (cm, s, g), plus the pairs in
common use (``nm_to_AA``, ``pc_to_kpc``, ``year_to_Myr`` ...). ``year`` is the
Julian year (365.25 days), as in astropy.

Examples
--------
>>> from jaff.physics import units
>>> units.eV_to_K(1.0)
11604.518121550082
>>> units.nm_to_eV(500.0)
2.479683968664005
>>> units.pc_to_cm(1.0)
3.085677581491367e+18
"""

from __future__ import annotations

from typing import Any

import astropy.units as _u

#: A number, a NumPy array or a SymPy expression.
Value = Any

_EQUIV = _u.spectral() + _u.temperature_energy()


def _factor(src: Any, dst: Any) -> float:
    """The value of 1 *src* in *dst*."""
    return float((1 * src).to(dst, equivalencies=_EQUIV).value)


# Factors: the value of 1 unit of the first in the second.
_EV_TO_ERG = _factor(_u.eV, _u.erg)
_ERG_TO_EV = _factor(_u.erg, _u.eV)
_KEV_TO_EV = _factor(_u.keV, _u.eV)
_EV_TO_KEV = _factor(_u.eV, _u.keV)
_J_TO_ERG = _factor(_u.J, _u.erg)
_ERG_TO_J = _factor(_u.erg, _u.J)
_EV_TO_K = _factor(_u.eV, _u.K)
_K_TO_EV = _factor(_u.K, _u.eV)
_K_TO_ERG = _factor(_u.K, _u.erg)
_ERG_TO_K = _factor(_u.erg, _u.K)
_HZ_TO_EV = _factor(_u.Hz, _u.eV)
_EV_TO_HZ = _factor(_u.eV, _u.Hz)
_HZ_TO_ERG = _factor(_u.Hz, _u.erg)
_ERG_TO_HZ = _factor(_u.erg, _u.Hz)
_CM1_TO_EV = _factor(1 / _u.cm, _u.eV)
_EV_TO_CM1 = _factor(_u.eV, 1 / _u.cm)
_HZ_TO_CM1 = _factor(_u.Hz, 1 / _u.cm)
_CM1_TO_HZ = _factor(1 / _u.cm, _u.Hz)
_M_TO_CM = _factor(_u.m, _u.cm)
_CM_TO_M = _factor(_u.cm, _u.m)
_KM_TO_CM = _factor(_u.km, _u.cm)
_CM_TO_KM = _factor(_u.cm, _u.km)
_NM_TO_CM = _factor(_u.nm, _u.cm)
_CM_TO_NM = _factor(_u.cm, _u.nm)
_UM_TO_CM = _factor(_u.um, _u.cm)
_CM_TO_UM = _factor(_u.cm, _u.um)
_AA_TO_CM = _factor(_u.AA, _u.cm)
_CM_TO_AA = _factor(_u.cm, _u.AA)
_AU_TO_CM = _factor(_u.au, _u.cm)
_CM_TO_AU = _factor(_u.cm, _u.au)
_PC_TO_CM = _factor(_u.pc, _u.cm)
_CM_TO_PC = _factor(_u.cm, _u.pc)
_KPC_TO_CM = _factor(_u.kpc, _u.cm)
_CM_TO_KPC = _factor(_u.cm, _u.kpc)
_MPC_TO_CM = _factor(_u.Mpc, _u.cm)
_CM_TO_MPC = _factor(_u.cm, _u.Mpc)
_LY_TO_CM = _factor(_u.lyr, _u.cm)
_CM_TO_LY = _factor(_u.cm, _u.lyr)
_M_TO_KM = _factor(_u.m, _u.km)
_KM_TO_M = _factor(_u.km, _u.m)
_NM_TO_AA = _factor(_u.nm, _u.AA)
_AA_TO_NM = _factor(_u.AA, _u.nm)
_NM_TO_UM = _factor(_u.nm, _u.um)
_UM_TO_NM = _factor(_u.um, _u.nm)
_UM_TO_AA = _factor(_u.um, _u.AA)
_AA_TO_UM = _factor(_u.AA, _u.um)
_AU_TO_PC = _factor(_u.au, _u.pc)
_PC_TO_AU = _factor(_u.pc, _u.au)
_PC_TO_KPC = _factor(_u.pc, _u.kpc)
_KPC_TO_PC = _factor(_u.kpc, _u.pc)
_KPC_TO_MPC = _factor(_u.kpc, _u.Mpc)
_MPC_TO_KPC = _factor(_u.Mpc, _u.kpc)
_PC_TO_LY = _factor(_u.pc, _u.lyr)
_LY_TO_PC = _factor(_u.lyr, _u.pc)
_MINUTE_TO_S = _factor(_u.min, _u.s)
_S_TO_MINUTE = _factor(_u.s, _u.min)
_HOUR_TO_S = _factor(_u.h, _u.s)
_S_TO_HOUR = _factor(_u.s, _u.h)
_DAY_TO_S = _factor(_u.d, _u.s)
_S_TO_DAY = _factor(_u.s, _u.d)
_YEAR_TO_S = _factor(_u.yr, _u.s)
_S_TO_YEAR = _factor(_u.s, _u.yr)
_KYR_TO_S = _factor(_u.kyr, _u.s)
_S_TO_KYR = _factor(_u.s, _u.kyr)
_MYR_TO_S = _factor(_u.Myr, _u.s)
_S_TO_MYR = _factor(_u.s, _u.Myr)
_GYR_TO_S = _factor(_u.Gyr, _u.s)
_S_TO_GYR = _factor(_u.s, _u.Gyr)
_DAY_TO_YEAR = _factor(_u.d, _u.yr)
_YEAR_TO_DAY = _factor(_u.yr, _u.d)
_YEAR_TO_MYR = _factor(_u.yr, _u.Myr)
_MYR_TO_YEAR = _factor(_u.Myr, _u.yr)
_KYR_TO_MYR = _factor(_u.kyr, _u.Myr)
_MYR_TO_KYR = _factor(_u.Myr, _u.kyr)
_MYR_TO_GYR = _factor(_u.Myr, _u.Gyr)
_GYR_TO_MYR = _factor(_u.Gyr, _u.Myr)
_KG_TO_G = _factor(_u.kg, _u.g)
_G_TO_KG = _factor(_u.g, _u.kg)
_AMU_TO_G = _factor(_u.u, _u.g)
_G_TO_AMU = _factor(_u.g, _u.u)
_MSUN_TO_G = _factor(_u.M_sun, _u.g)
_G_TO_MSUN = _factor(_u.g, _u.M_sun)
_KM_S_TO_PC_MYR = _factor(_u.km / _u.s, _u.pc / _u.Myr)
_PC_MYR_TO_KM_S = _factor(_u.pc / _u.Myr, _u.km / _u.s)
# Photons: the energy / frequency at a wavelength of 1, divided by the value.
_HC_EV_NM = _factor(_u.nm, _u.eV)
_C_HZ_NM = _factor(_u.nm, _u.Hz)
_HC_EV_AA = _factor(_u.AA, _u.eV)
_C_HZ_AA = _factor(_u.AA, _u.Hz)
_HC_EV_UM = _factor(_u.um, _u.eV)
_C_HZ_UM = _factor(_u.um, _u.Hz)

__all__ = [
    "eV_to_erg",
    "erg_to_eV",
    "keV_to_eV",
    "eV_to_keV",
    "J_to_erg",
    "erg_to_J",
    "eV_to_K",
    "K_to_eV",
    "K_to_erg",
    "erg_to_K",
    "Hz_to_eV",
    "eV_to_Hz",
    "Hz_to_erg",
    "erg_to_Hz",
    "cm1_to_eV",
    "eV_to_cm1",
    "Hz_to_cm1",
    "cm1_to_Hz",
    "m_to_cm",
    "cm_to_m",
    "km_to_cm",
    "cm_to_km",
    "nm_to_cm",
    "cm_to_nm",
    "um_to_cm",
    "cm_to_um",
    "AA_to_cm",
    "cm_to_AA",
    "au_to_cm",
    "cm_to_au",
    "pc_to_cm",
    "cm_to_pc",
    "kpc_to_cm",
    "cm_to_kpc",
    "Mpc_to_cm",
    "cm_to_Mpc",
    "ly_to_cm",
    "cm_to_ly",
    "m_to_km",
    "km_to_m",
    "nm_to_AA",
    "AA_to_nm",
    "nm_to_um",
    "um_to_nm",
    "um_to_AA",
    "AA_to_um",
    "au_to_pc",
    "pc_to_au",
    "pc_to_kpc",
    "kpc_to_pc",
    "kpc_to_Mpc",
    "Mpc_to_kpc",
    "pc_to_ly",
    "ly_to_pc",
    "minute_to_s",
    "s_to_minute",
    "hour_to_s",
    "s_to_hour",
    "day_to_s",
    "s_to_day",
    "year_to_s",
    "s_to_year",
    "kyr_to_s",
    "s_to_kyr",
    "Myr_to_s",
    "s_to_Myr",
    "Gyr_to_s",
    "s_to_Gyr",
    "day_to_year",
    "year_to_day",
    "year_to_Myr",
    "Myr_to_year",
    "kyr_to_Myr",
    "Myr_to_kyr",
    "Myr_to_Gyr",
    "Gyr_to_Myr",
    "kg_to_g",
    "g_to_kg",
    "amu_to_g",
    "g_to_amu",
    "Msun_to_g",
    "g_to_Msun",
    "km_s_to_pc_Myr",
    "pc_Myr_to_km_s",
    "nm_to_eV",
    "eV_to_nm",
    "nm_to_Hz",
    "Hz_to_nm",
    "AA_to_eV",
    "eV_to_AA",
    "AA_to_Hz",
    "Hz_to_AA",
    "um_to_eV",
    "eV_to_um",
    "um_to_Hz",
    "Hz_to_um",
]

# =============================================================================
# Energy, temperature, frequency and wavenumber
# =============================================================================


def eV_to_erg(x: Value) -> Value:
    """Convert an energy from eV to erg.

    Multiplies by 1.602176634e-12.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The energy in eV.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The energy in erg.
    """
    return x * _EV_TO_ERG


def erg_to_eV(x: Value) -> Value:
    """Convert an energy from erg to eV.

    Multiplies by 6.241509074e+11.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The energy in erg.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The energy in eV.
    """
    return x * _ERG_TO_EV


def keV_to_eV(x: Value) -> Value:
    """Convert an energy from keV to eV.

    Multiplies by 1000.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The energy in keV.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The energy in eV.
    """
    return x * _KEV_TO_EV


def eV_to_keV(x: Value) -> Value:
    """Convert an energy from eV to keV.

    Multiplies by 0.001.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The energy in eV.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The energy in keV.
    """
    return x * _EV_TO_KEV


def J_to_erg(x: Value) -> Value:
    """Convert an energy from J to erg.

    Multiplies by 10000000.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The energy in J.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The energy in erg.
    """
    return x * _J_TO_ERG


def erg_to_J(x: Value) -> Value:
    """Convert an energy from erg to J.

    Multiplies by 1e-07.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The energy in erg.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The energy in J.
    """
    return x * _ERG_TO_J


def eV_to_K(x: Value) -> Value:
    """Convert an energy in eV to a temperature in K.

    Multiplies by 11604.51812 (E = k_B T).

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The energy in eV.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The temperature in K.
    """
    return x * _EV_TO_K


def K_to_eV(x: Value) -> Value:
    """Convert a temperature in K to an energy in eV.

    Multiplies by 8.617333262e-05 (E = k_B T).

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The temperature in K.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The energy in eV.
    """
    return x * _K_TO_EV


def K_to_erg(x: Value) -> Value:
    """Convert a temperature in K to an energy in erg.

    Multiplies by 1.380649e-16 (E = k_B T).

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The temperature in K.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The energy in erg.
    """
    return x * _K_TO_ERG


def erg_to_K(x: Value) -> Value:
    """Convert an energy in erg to a temperature in K.

    Multiplies by 7.242970516e+15 (E = k_B T).

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The energy in erg.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The temperature in K.
    """
    return x * _ERG_TO_K


def Hz_to_eV(x: Value) -> Value:
    """Convert a frequency in Hz to an energy in eV.

    Multiplies by 4.135667697e-15 (E = h nu).

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The frequency in Hz.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The energy in eV.
    """
    return x * _HZ_TO_EV


def eV_to_Hz(x: Value) -> Value:
    """Convert an energy in eV to a frequency in Hz.

    Multiplies by 2.417989242e+14 (E = h nu).

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The energy in eV.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The frequency in Hz.
    """
    return x * _EV_TO_HZ


def Hz_to_erg(x: Value) -> Value:
    """Convert a frequency in Hz to an energy in erg.

    Multiplies by 6.62607015e-27 (E = h nu).

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The frequency in Hz.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The energy in erg.
    """
    return x * _HZ_TO_ERG


def erg_to_Hz(x: Value) -> Value:
    """Convert an energy in erg to a frequency in Hz.

    Multiplies by 1.50919018e+26 (E = h nu).

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The energy in erg.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The frequency in Hz.
    """
    return x * _ERG_TO_HZ


def cm1_to_eV(x: Value) -> Value:
    """Convert a wavenumber in cm⁻¹ to an energy in eV.

    Multiplies by 0.0001239841984 (E = h c / lambda, with 1 / lambda the wavenumber).

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The wavenumber in cm⁻¹.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The energy in eV.
    """
    return x * _CM1_TO_EV


def eV_to_cm1(x: Value) -> Value:
    """Convert an energy in eV to a wavenumber in cm⁻¹.

    Multiplies by 8065.543937 (E = h c / lambda, with 1 / lambda the wavenumber).

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The energy in eV.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The wavenumber in cm⁻¹.
    """
    return x * _EV_TO_CM1


def Hz_to_cm1(x: Value) -> Value:
    """Convert a frequency in Hz to a wavenumber in cm⁻¹.

    Multiplies by 3.335640952e-11 (1 / lambda = nu / c).

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The frequency in Hz.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The wavenumber in cm⁻¹.
    """
    return x * _HZ_TO_CM1


def cm1_to_Hz(x: Value) -> Value:
    """Convert a wavenumber in cm⁻¹ to a frequency in Hz.

    Multiplies by 2.99792458e+10 (1 / lambda = nu / c).

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The wavenumber in cm⁻¹.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The frequency in Hz.
    """
    return x * _CM1_TO_HZ


# =============================================================================
# Length
# =============================================================================


def m_to_cm(x: Value) -> Value:
    """Convert a length from m to cm.

    Multiplies by 100.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in m.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in cm.
    """
    return x * _M_TO_CM


def cm_to_m(x: Value) -> Value:
    """Convert a length from cm to m.

    Multiplies by 0.01.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in cm.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in m.
    """
    return x * _CM_TO_M


def km_to_cm(x: Value) -> Value:
    """Convert a length from km to cm.

    Multiplies by 100000.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in km.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in cm.
    """
    return x * _KM_TO_CM


def cm_to_km(x: Value) -> Value:
    """Convert a length from cm to km.

    Multiplies by 1e-05.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in cm.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in km.
    """
    return x * _CM_TO_KM


def nm_to_cm(x: Value) -> Value:
    """Convert a length from nm to cm.

    Multiplies by 1e-07.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in nm.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in cm.
    """
    return x * _NM_TO_CM


def cm_to_nm(x: Value) -> Value:
    """Convert a length from cm to nm.

    Multiplies by 10000000.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in cm.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in nm.
    """
    return x * _CM_TO_NM


def um_to_cm(x: Value) -> Value:
    """Convert a length from µm to cm.

    Multiplies by 0.0001.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in µm.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in cm.
    """
    return x * _UM_TO_CM


def cm_to_um(x: Value) -> Value:
    """Convert a length from cm to µm.

    Multiplies by 10000.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in cm.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in µm.
    """
    return x * _CM_TO_UM


def AA_to_cm(x: Value) -> Value:
    """Convert a length from Å to cm.

    Multiplies by 1e-08.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in Å.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in cm.
    """
    return x * _AA_TO_CM


def cm_to_AA(x: Value) -> Value:
    """Convert a length from cm to Å.

    Multiplies by 100000000.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in cm.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in Å.
    """
    return x * _CM_TO_AA


def au_to_cm(x: Value) -> Value:
    """Convert a length from au to cm.

    Multiplies by 1.495978707e+13.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in au.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in cm.
    """
    return x * _AU_TO_CM


def cm_to_au(x: Value) -> Value:
    """Convert a length from cm to au.

    Multiplies by 6.684587122e-14.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in cm.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in au.
    """
    return x * _CM_TO_AU


def pc_to_cm(x: Value) -> Value:
    """Convert a length from pc to cm.

    Multiplies by 3.085677581e+18.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in pc.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in cm.
    """
    return x * _PC_TO_CM


def cm_to_pc(x: Value) -> Value:
    """Convert a length from cm to pc.

    Multiplies by 3.240779289e-19.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in cm.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in pc.
    """
    return x * _CM_TO_PC


def kpc_to_cm(x: Value) -> Value:
    """Convert a length from kpc to cm.

    Multiplies by 3.085677581e+21.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in kpc.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in cm.
    """
    return x * _KPC_TO_CM


def cm_to_kpc(x: Value) -> Value:
    """Convert a length from cm to kpc.

    Multiplies by 3.240779289e-22.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in cm.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in kpc.
    """
    return x * _CM_TO_KPC


def Mpc_to_cm(x: Value) -> Value:
    """Convert a length from Mpc to cm.

    Multiplies by 3.085677581e+24.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in Mpc.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in cm.
    """
    return x * _MPC_TO_CM


def cm_to_Mpc(x: Value) -> Value:
    """Convert a length from cm to Mpc.

    Multiplies by 3.240779289e-25.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in cm.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in Mpc.
    """
    return x * _CM_TO_MPC


def ly_to_cm(x: Value) -> Value:
    """Convert a length from light-years to cm.

    Multiplies by 9.460730473e+17.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in light-years.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in cm.
    """
    return x * _LY_TO_CM


def cm_to_ly(x: Value) -> Value:
    """Convert a length from cm to light-years.

    Multiplies by 1.057000834e-18.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in cm.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in light-years.
    """
    return x * _CM_TO_LY


def m_to_km(x: Value) -> Value:
    """Convert a length from m to km.

    Multiplies by 0.001.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in m.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in km.
    """
    return x * _M_TO_KM


def km_to_m(x: Value) -> Value:
    """Convert a length from km to m.

    Multiplies by 1000.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in km.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in m.
    """
    return x * _KM_TO_M


def nm_to_AA(x: Value) -> Value:
    """Convert a length from nm to Å.

    Multiplies by 10.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in nm.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in Å.
    """
    return x * _NM_TO_AA


def AA_to_nm(x: Value) -> Value:
    """Convert a length from Å to nm.

    Multiplies by 0.1.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in Å.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in nm.
    """
    return x * _AA_TO_NM


def nm_to_um(x: Value) -> Value:
    """Convert a length from nm to µm.

    Multiplies by 0.001.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in nm.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in µm.
    """
    return x * _NM_TO_UM


def um_to_nm(x: Value) -> Value:
    """Convert a length from µm to nm.

    Multiplies by 1000.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in µm.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in nm.
    """
    return x * _UM_TO_NM


def um_to_AA(x: Value) -> Value:
    """Convert a length from µm to Å.

    Multiplies by 10000.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in µm.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in Å.
    """
    return x * _UM_TO_AA


def AA_to_um(x: Value) -> Value:
    """Convert a length from Å to µm.

    Multiplies by 0.0001.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in Å.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in µm.
    """
    return x * _AA_TO_UM


def au_to_pc(x: Value) -> Value:
    """Convert a length from au to pc.

    Multiplies by 4.848136811e-06.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in au.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in pc.
    """
    return x * _AU_TO_PC


def pc_to_au(x: Value) -> Value:
    """Convert a length from pc to au.

    Multiplies by 206264.8062.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in pc.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in au.
    """
    return x * _PC_TO_AU


def pc_to_kpc(x: Value) -> Value:
    """Convert a length from pc to kpc.

    Multiplies by 0.001.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in pc.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in kpc.
    """
    return x * _PC_TO_KPC


def kpc_to_pc(x: Value) -> Value:
    """Convert a length from kpc to pc.

    Multiplies by 1000.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in kpc.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in pc.
    """
    return x * _KPC_TO_PC


def kpc_to_Mpc(x: Value) -> Value:
    """Convert a length from kpc to Mpc.

    Multiplies by 0.001.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in kpc.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in Mpc.
    """
    return x * _KPC_TO_MPC


def Mpc_to_kpc(x: Value) -> Value:
    """Convert a length from Mpc to kpc.

    Multiplies by 1000.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in Mpc.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in kpc.
    """
    return x * _MPC_TO_KPC


def pc_to_ly(x: Value) -> Value:
    """Convert a length from pc to light-years.

    Multiplies by 3.261563777.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in pc.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in light-years.
    """
    return x * _PC_TO_LY


def ly_to_pc(x: Value) -> Value:
    """Convert a length from light-years to pc.

    Multiplies by 0.3066013938.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The length in light-years.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The length in pc.
    """
    return x * _LY_TO_PC


# =============================================================================
# Time
# =============================================================================


def minute_to_s(x: Value) -> Value:
    """Convert a time from minutes to s.

    Multiplies by 60.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The time in minutes.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The time in s.
    """
    return x * _MINUTE_TO_S


def s_to_minute(x: Value) -> Value:
    """Convert a time from s to minutes.

    Multiplies by 0.01666666667.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The time in s.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The time in minutes.
    """
    return x * _S_TO_MINUTE


def hour_to_s(x: Value) -> Value:
    """Convert a time from hours to s.

    Multiplies by 3600.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The time in hours.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The time in s.
    """
    return x * _HOUR_TO_S


def s_to_hour(x: Value) -> Value:
    """Convert a time from s to hours.

    Multiplies by 0.0002777777778.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The time in s.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The time in hours.
    """
    return x * _S_TO_HOUR


def day_to_s(x: Value) -> Value:
    """Convert a time from days to s.

    Multiplies by 86400.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The time in days.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The time in s.
    """
    return x * _DAY_TO_S


def s_to_day(x: Value) -> Value:
    """Convert a time from s to days.

    Multiplies by 1.157407407e-05.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The time in s.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The time in days.
    """
    return x * _S_TO_DAY


def year_to_s(x: Value) -> Value:
    """Convert a time from Julian years to s.

    Multiplies by 31557600.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The time in Julian years.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The time in s.
    """
    return x * _YEAR_TO_S


def s_to_year(x: Value) -> Value:
    """Convert a time from s to Julian years.

    Multiplies by 3.168808781e-08.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The time in s.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The time in Julian years.
    """
    return x * _S_TO_YEAR


def kyr_to_s(x: Value) -> Value:
    """Convert a time from kyr to s.

    Multiplies by 3.15576e+10.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The time in kyr.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The time in s.
    """
    return x * _KYR_TO_S


def s_to_kyr(x: Value) -> Value:
    """Convert a time from s to kyr.

    Multiplies by 3.168808781e-11.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The time in s.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The time in kyr.
    """
    return x * _S_TO_KYR


def Myr_to_s(x: Value) -> Value:
    """Convert a time from Myr to s.

    Multiplies by 3.15576e+13.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The time in Myr.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The time in s.
    """
    return x * _MYR_TO_S


def s_to_Myr(x: Value) -> Value:
    """Convert a time from s to Myr.

    Multiplies by 3.168808781e-14.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The time in s.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The time in Myr.
    """
    return x * _S_TO_MYR


def Gyr_to_s(x: Value) -> Value:
    """Convert a time from Gyr to s.

    Multiplies by 3.15576e+16.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The time in Gyr.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The time in s.
    """
    return x * _GYR_TO_S


def s_to_Gyr(x: Value) -> Value:
    """Convert a time from s to Gyr.

    Multiplies by 3.168808781e-17.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The time in s.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The time in Gyr.
    """
    return x * _S_TO_GYR


def day_to_year(x: Value) -> Value:
    """Convert a time from days to Julian years.

    Multiplies by 0.002737850787.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The time in days.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The time in Julian years.
    """
    return x * _DAY_TO_YEAR


def year_to_day(x: Value) -> Value:
    """Convert a time from Julian years to days.

    Multiplies by 365.25.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The time in Julian years.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The time in days.
    """
    return x * _YEAR_TO_DAY


def year_to_Myr(x: Value) -> Value:
    """Convert a time from Julian years to Myr.

    Multiplies by 1e-06.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The time in Julian years.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The time in Myr.
    """
    return x * _YEAR_TO_MYR


def Myr_to_year(x: Value) -> Value:
    """Convert a time from Myr to Julian years.

    Multiplies by 1000000.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The time in Myr.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The time in Julian years.
    """
    return x * _MYR_TO_YEAR


def kyr_to_Myr(x: Value) -> Value:
    """Convert a time from kyr to Myr.

    Multiplies by 0.001.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The time in kyr.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The time in Myr.
    """
    return x * _KYR_TO_MYR


def Myr_to_kyr(x: Value) -> Value:
    """Convert a time from Myr to kyr.

    Multiplies by 1000.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The time in Myr.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The time in kyr.
    """
    return x * _MYR_TO_KYR


def Myr_to_Gyr(x: Value) -> Value:
    """Convert a time from Myr to Gyr.

    Multiplies by 0.001.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The time in Myr.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The time in Gyr.
    """
    return x * _MYR_TO_GYR


def Gyr_to_Myr(x: Value) -> Value:
    """Convert a time from Gyr to Myr.

    Multiplies by 1000.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The time in Gyr.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The time in Myr.
    """
    return x * _GYR_TO_MYR


# =============================================================================
# Mass
# =============================================================================


def kg_to_g(x: Value) -> Value:
    """Convert a mass from kg to g.

    Multiplies by 1000.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The mass in kg.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The mass in g.
    """
    return x * _KG_TO_G


def g_to_kg(x: Value) -> Value:
    """Convert a mass from g to kg.

    Multiplies by 0.001.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The mass in g.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The mass in kg.
    """
    return x * _G_TO_KG


def amu_to_g(x: Value) -> Value:
    """Convert a mass from atomic mass units to g.

    Multiplies by 1.660539069e-24.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The mass in atomic mass units.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The mass in g.
    """
    return x * _AMU_TO_G


def g_to_amu(x: Value) -> Value:
    """Convert a mass from g to atomic mass units.

    Multiplies by 6.022140754e+23.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The mass in g.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The mass in atomic mass units.
    """
    return x * _G_TO_AMU


def Msun_to_g(x: Value) -> Value:
    """Convert a mass from solar masses to g.

    Multiplies by 1.988409871e+33.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The mass in solar masses.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The mass in g.
    """
    return x * _MSUN_TO_G


def g_to_Msun(x: Value) -> Value:
    """Convert a mass from g to solar masses.

    Multiplies by 5.029144216e-34.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The mass in g.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The mass in solar masses.
    """
    return x * _G_TO_MSUN


# =============================================================================
# Velocity
# =============================================================================


def km_s_to_pc_Myr(x: Value) -> Value:
    """Convert a velocity from km s⁻¹ to pc Myr⁻¹.

    Multiplies by 1.022712165.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The velocity in km s⁻¹.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The velocity in pc Myr⁻¹.
    """
    return x * _KM_S_TO_PC_MYR


def pc_Myr_to_km_s(x: Value) -> Value:
    """Convert a velocity from pc Myr⁻¹ to km s⁻¹.

    Multiplies by 0.9777922217.

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The velocity in pc Myr⁻¹.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The velocity in km s⁻¹.
    """
    return x * _PC_MYR_TO_KM_S


# =============================================================================
# Photon wavelength, energy and frequency (inversely related)
# =============================================================================


def nm_to_eV(x: Value) -> Value:
    """Convert a photon wavelength in nm to its energy in eV.

    Divides 1239.841984 by *x* (E = h c / lambda).

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The photon wavelength in nm.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The photon energy in eV.
    """
    return _HC_EV_NM / x


def eV_to_nm(x: Value) -> Value:
    """Convert a photon energy in eV to its wavelength in nm.

    Divides 1239.841984 by *x* (E = h c / lambda).

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The photon energy in eV.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The photon wavelength in nm.
    """
    return _HC_EV_NM / x


def nm_to_Hz(x: Value) -> Value:
    """Convert a photon wavelength in nm to its frequency in Hz.

    Divides 2.99792458e+17 by *x* (nu = c / lambda).

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The photon wavelength in nm.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The photon frequency in Hz.
    """
    return _C_HZ_NM / x


def Hz_to_nm(x: Value) -> Value:
    """Convert a photon frequency in Hz to its wavelength in nm.

    Divides 2.99792458e+17 by *x* (nu = c / lambda).

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The photon frequency in Hz.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The photon wavelength in nm.
    """
    return _C_HZ_NM / x


def AA_to_eV(x: Value) -> Value:
    """Convert a photon wavelength in Å to its energy in eV.

    Divides 12398.41984 by *x* (E = h c / lambda).

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The photon wavelength in Å.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The photon energy in eV.
    """
    return _HC_EV_AA / x


def eV_to_AA(x: Value) -> Value:
    """Convert a photon energy in eV to its wavelength in Å.

    Divides 12398.41984 by *x* (E = h c / lambda).

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The photon energy in eV.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The photon wavelength in Å.
    """
    return _HC_EV_AA / x


def AA_to_Hz(x: Value) -> Value:
    """Convert a photon wavelength in Å to its frequency in Hz.

    Divides 2.99792458e+18 by *x* (nu = c / lambda).

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The photon wavelength in Å.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The photon frequency in Hz.
    """
    return _C_HZ_AA / x


def Hz_to_AA(x: Value) -> Value:
    """Convert a photon frequency in Hz to its wavelength in Å.

    Divides 2.99792458e+18 by *x* (nu = c / lambda).

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The photon frequency in Hz.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The photon wavelength in Å.
    """
    return _C_HZ_AA / x


def um_to_eV(x: Value) -> Value:
    """Convert a photon wavelength in µm to its energy in eV.

    Divides 1.239841984 by *x* (E = h c / lambda).

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The photon wavelength in µm.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The photon energy in eV.
    """
    return _HC_EV_UM / x


def eV_to_um(x: Value) -> Value:
    """Convert a photon energy in eV to its wavelength in µm.

    Divides 1.239841984 by *x* (E = h c / lambda).

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The photon energy in eV.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The photon wavelength in µm.
    """
    return _HC_EV_UM / x


def um_to_Hz(x: Value) -> Value:
    """Convert a photon wavelength in µm to its frequency in Hz.

    Divides 2.99792458e+14 by *x* (nu = c / lambda).

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The photon wavelength in µm.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The photon frequency in Hz.
    """
    return _C_HZ_UM / x


def Hz_to_um(x: Value) -> Value:
    """Convert a photon frequency in Hz to its wavelength in µm.

    Divides 2.99792458e+14 by *x* (nu = c / lambda).

    Parameters
    ----------
    x : float, numpy.ndarray or sympy.Expr
        The photon frequency in Hz.

    Returns
    -------
    float, numpy.ndarray or sympy.Expr
        The photon wavelength in µm.
    """
    return _C_HZ_UM / x
