---
tags:
    - Api
icon: phosphor/atom
---

# jaff.physics

Physical constants and photochemical cross-section lookup for astrochemical calculations.

## Classes

| Class            | Description                                                       |
| ---------------- | ----------------------------------------------------------------- |
| `Photochemistry` | Photo cross-section lookup from the bundled databases (see below) |
| [`RadiationProps`](radiation/index.md#radiationprops) | Radiation-field configuration (bands, per-band spectral index, mode, speed of light) |
| [`Radiation`](radiation/index.md#radiation_1) | Band collection built by `Network` (`net.radiation`); piecewise spectral profiles |
| [`RadiationGroup`](radiation/index.md#radiationgroup) | One band: edges, `profile_idx`, `photden`, `eavg`, per-reaction rate data |

## Submodules

| Submodule                         | Description                                                           |
| --------------------------------- | --------------------------------------------------------------------- |
| [`constants`](constants/index.md) | Physical & astronomical constants as `astropy` Quantities (see below) |

## Photochemistry methods

`jaff.physics.Photochemistry` resolves a reaction's cross sections and
shielding factors from `jaff.db`. Constructing it downloads the cross-section
and line-shielding data files on first use (cached thereafter), so instantiate
once and reuse:

```python
from jaff.physics import Photochemistry

photo = Photochemistry()
photo.get_xsec(rxn)              # XsecsProps from the tabulated databases
photo.get_verner_xsec(rxn)       # analytic Verner σ(E) (sympy) or None
Photochemistry.shielding(rxn, net)  # symbolic shielding factor (sympy)
```

| Method                          | Returns               | Description                                                                                                              |
| ------------------------------- | --------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| `get_xsec(reaction, required=False)` | `XsecsProps or None` | Tabulated (Leiden / NORAD: `photon_energy` in eV, `photo_absorption` / `photodecay` in cm²) or Verner symbolic (`photodecay_expr`) cross sections, per `pi_database`; carries `database` |
| `get_verner_xsec(reaction)`     | `sympy.Basic or None` | Analytic Verner (1996) σ(E) expression (symbol `E` in eV, σ in cm²)                                                     |
| `shielding(reaction, network)`  | `sympy.Expr`          | Dimensionless line-shielding factor; dispatches to the global/local shielding function named by the reaction metadata   |

## Unit systems

`jaff.physics.constants` exposes physical constants as module-level
[`astropy.units.Quantity`](https://docs.astropy.org/en/stable/units/quantity.html)
objects. There are no separate `cgs`/`si`/`gaussian`/`natural` tables — select
the unit system **per quantity** at the call site with `.cgs`, `.si`, or
`.to(...)`. See [`constants`](constants/index.md) for the full list.

## Example

```python
from jaff.physics import constants as c

# Per-quantity unit selection
c.c.cgs              # speed of light [cm/s]
c.k_B.cgs            # Boltzmann constant [erg/K]
c.c.si               # speed of light [m/s]
c.k_B.to("eV / K")   # Boltzmann constant [eV/K]

# Elementary charge keeps astropy's EM-unit views
c.e.esu              # Gaussian / CGS-ESU [Fr]
c.e.si               # SI [C]

# Bare float (in the chosen unit) when a number is required
float(c.m_e.cgs.value)   # electron mass [g]
```
