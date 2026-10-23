# pyright: reportAttributeAccessIssue=false
"""Physical and astronomical constants as :class:`astropy.units.Quantity`.

This module is the single source of truth for physical constants in JAFF.
Every constant is an :class:`astropy.units.Quantity` carrying both a value and
a unit, sourced from :mod:`astropy.constants` wherever astropy ships it (the
CODATA / IAU releases of the installed astropy; CODATA 2022 and IAU 2015 with
astropy 8), and computed from those primitives otherwise.  Names match the
:mod:`astropy.constants` spelling (``k_B``, ``m_p``, ``M_sun``, ``N_A``, ``u``,
``R`` ...).  The values listed with each constant are its CGS view with
astropy 8.

Notes
-----
**Unit systems.** The legacy ``cgs`` / ``si`` / ``gaussian`` / ``natural``
constant tables are gone.  Because every constant is now a unit-aware
``Quantity``, the unit system is selected *per quantity* at the call site with
``.cgs``, ``.si`` or ``.to(...)``; take ``.value`` for a bare float.  The
elementary charge ``e`` keeps astropy's electromagnetic-unit views: use
``e.esu`` (Gaussian/CGS-ESU, Franklin) or ``e.si`` (Coulomb); it has no
``.cgs`` view.

**Migration** from the old ``Constants`` tables:

| Old | New |
| --- | --- |
| ``constants.cgs.c`` | ``constants.c.cgs`` |
| ``constants.cgs.k_b`` | ``constants.k_B.cgs`` (astropy spelling) |
| ``constants.si.m_e`` | ``constants.m_e.si`` |
| ``constants.cgs.ev_to_erg`` | ``(1 * u.eV).to("erg")`` (``u`` = :mod:`astropy.units`) |
| ``constants.cgs.kb_ev`` | ``constants.k_B.to("eV / K")`` |

The ``gaussian`` and ``natural`` systems and the trivial conversion factors
(``ev_to_erg``, ``kb_ev``, ``Ry_hc``) were removed; derive them on demand with
``.to(...)``.

Examples
--------
>>> from jaff.physics import constants as c
>>> c.c.cgs  # speed of light in CGS
<Quantity 2.99792458e+10 cm / s>
>>> c.c.to("m / s")  # ... in SI (c.c.si is the astropy Constant itself)
<Quantity 2.99792458e+08 m / s>
>>> c.k_B.to("eV / K")  # ... in any explicit unit
<Quantity 8.61733326e-05 eV / K>
>>> float(c.m_p.cgs.value)  # bare float when a number is required
1.67262192595e-24
>>> float(c.e.esu.value)  # Gaussian/CGS-ESU (Franklin)
4.803204712570263e-10
>>> float(c.e.si.value)  # SI (Coulomb)
1.602176634e-19
"""

from __future__ import annotations

import astropy.constants as _ac
import astropy.units as _u

# =============================================================================
# Fundamental physical constants (astropy / CODATA)
# =============================================================================
c = _ac.c
"""Speed of light in vacuum; 2.99792458e10 cm s⁻¹."""
h = _ac.h
"""Planck constant; 6.62607015e-27 erg s."""
hbar = _ac.hbar
"""Reduced Planck constant (h / 2pi); 1.05457182e-27 erg s."""
G = _ac.G
"""Newtonian gravitational constant; 6.6743e-8 cm³ g⁻¹ s⁻²."""
k_B = _ac.k_B
"""Boltzmann constant; 1.380649e-16 erg K⁻¹."""
e = _ac.e
"""Elementary charge; 4.80320471e-10 Fr (``.esu``) or 1.602176634e-19 C (``.si``)."""
m_e = _ac.m_e
"""Electron rest mass; 9.10938371e-28 g."""
m_p = _ac.m_p
"""Proton rest mass; 1.67262193e-24 g."""
m_n = _ac.m_n
"""Neutron rest mass; 1.6749275e-24 g."""
u = _ac.u
"""Atomic mass unit (1 u); 1.66053907e-24 g."""
me_mp = (m_e / m_p).decompose()
"""Electron-to-proton mass ratio (dimensionless); 5.4462e-4."""

# =============================================================================
# Astronomical constants
# =============================================================================
M_sun = _ac.M_sun
"""Solar mass; 1.98840987e33 g."""
R_sun = _ac.R_sun
"""Solar radius (nominal); 6.957e10 cm."""
L_sun = _ac.L_sun
"""Solar luminosity (nominal); 3.828e33 erg s⁻¹."""
pc = _ac.pc
"""Parsec; 3.08567758e18 cm."""
kpc = _ac.kpc
"""Kiloparsec; 3.08567758e21 cm."""
Mpc = 1.0 * _u.Mpc
"""Megaparsec (stored as ``1 Mpc``); 3.08567758e24 cm."""
au = _ac.au
"""Astronomical unit; 1.49597871e13 cm."""
ly = 1.0 * _u.lyr
"""Light-year (stored as ``1 lyr``); 9.46073047e17 cm."""
H0 = 67.4 * _u.km / _u.s / _u.Mpc
"""Hubble constant (Planck 2018); 67.4 km s⁻¹ Mpc⁻¹ = 2.18428524e-18 s⁻¹."""

# =============================================================================
# Astrochemistry / radiation (computed from the primitives above)
# =============================================================================
sigma_sb = _ac.sigma_sb
"""Stefan-Boltzmann constant; 5.67037442e-5 erg cm⁻² s⁻¹ K⁻⁴."""
a_rad = (4.0 * sigma_sb / c).to("erg / (cm3 K4)")
"""Radiation constant (a = 4 sigma / c); 7.56573325e-15 erg cm⁻³ K⁻⁴."""
alpha = _ac.alpha
"""Fine-structure constant (dimensionless, ~1/137); 7.29735256e-3."""
sigma_T = _ac.sigma_T
"""Thomson scattering cross section; 6.65245871e-25 cm²."""
lambda_C = (h / (m_e * c)).to("cm")
"""Compton wavelength of the electron (h / m_e c); 2.42631024e-10 cm."""
a0 = _ac.a0
"""Bohr radius; 5.29177211e-9 cm."""
Ry = (_ac.Ryd * h * c).to("eV")
"""Rydberg energy (Ryd h c); 13.6056931 eV (stored in eV)."""
r_e = (alpha**2 * a0).to("cm")
"""Classical electron radius (alpha^2 * a0); 2.81794032e-13 cm."""
gyro_coeff = (e.esu / m_e).to("Fr / g")
"""Electron gyromagnetic coefficient (e / m_e) in CGS-ESU; 5.27280973e17 Fr g⁻¹."""

# =============================================================================
# Astrochemistry-specific constants
# =============================================================================
mu_H = m_p
"""Mean mass of a hydrogen atom (equal to ``m_p``); 1.67262193e-24 g."""
mu_H2 = 2.0 * m_p
"""Mean mass of a hydrogen molecule (``2 m_p``); 3.34524385e-24 g."""
nH_ref = (1.0 / m_p).to("1 / g")
"""Reference hydrogen number per unit mass (1 / m_p); 5.9786374e23 g⁻¹."""

# =============================================================================
# Temperature
# =============================================================================
T_cmb = 2.725 * _u.K
"""CMB temperature; 2.725 K."""
T_1ev = (1.0 * _u.eV / k_B).to("K")
"""Temperature equivalent of 1 eV (k_B T = 1 eV); 11604.5181 K."""

# =============================================================================
# Cross-section reference values
# =============================================================================
barn = 1.0 * _u.barn
"""1 barn; 1e-24 cm² (stored in barn)."""
mbarn = 1.0e6 * _u.barn
"""1 megabarn (1e6 barn); 1e-18 cm² (stored in barn)."""

# =============================================================================
# Reference
# =============================================================================
N_A = _ac.N_A
"""Avogadro constant; 6.02214076e23 mol⁻¹."""
R = _ac.R
"""Molar gas constant (k_B * N_A); 8.31446262e7 erg K⁻¹ mol⁻¹."""

__all__ = [
    "c",
    "h",
    "hbar",
    "G",
    "k_B",
    "e",
    "m_e",
    "m_p",
    "m_n",
    "u",
    "me_mp",
    "M_sun",
    "R_sun",
    "L_sun",
    "pc",
    "kpc",
    "Mpc",
    "au",
    "ly",
    "H0",
    "sigma_sb",
    "a_rad",
    "alpha",
    "sigma_T",
    "lambda_C",
    "a0",
    "Ry",
    "r_e",
    "gyro_coeff",
    "mu_H",
    "mu_H2",
    "nH_ref",
    "T_cmb",
    "T_1ev",
    "barn",
    "mbarn",
    "N_A",
    "R",
]
