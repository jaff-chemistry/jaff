---
tags:
    - Api
    - Network
---

# thermodynamics

`#!python net.thermodynamics`

The `Thermodynamics` object of a network: its equation of state and the heating
rates and temperature rates derived from it. All attributes are built lazily on
first access and cached.

```python
from jaff import Network
from jaff.physics import EosProps

net = Network("COthin", eos_props=EosProps("ideal", gamma=5.0 / 3.0))
th = net.thermodynamics
th.eos.specific      # erg/g
th.dEdt_tot          # chemical + extra heating/cooling
th.dTdt_tot          # dT/dt including the composition term
```

| Attribute       | Type            | Description                                                                              |
| --------------- | --------------- | ---------------------------------------------------------------------------------------- |
| `eos`           | `InternalEnergy`| Internal energy for the network's `eos_props` (ideal gas, $\gamma = 1.6666666666667$, by default) |
| `dEdt_chemical` | `DEDt`          | Chemical heating/cooling rate $\sum_r dE_r F_r$ (erg cm⁻³ s⁻¹)                           |
| `dEdt_extra`    | `DEDt`          | Non-reactive rate from the `heatingcoolingrate` auxiliary function (else `0`)            |
| `dEdt_tot`      | `DEDt`          | `dEdt_chemical + dEdt_extra`                                                             |
| `dTdt_chemical` | `sympy.Expr`    | $(\dot{E}_\mathrm{chem} - \sum_i \partial E/\partial n_i\, \dot{n}_i) / (\partial E/\partial T)$ (K s⁻¹) |
| `dTdt_extra`    | `sympy.Expr`    | $\dot{E}_\mathrm{extra} / (\partial E/\partial T)$ (K s⁻¹)                              |
| `dTdt_tot`      | `sympy.Expr`    | `dTdt_chemical + dTdt_extra`, built as one fraction                                      |

`DEDt` wraps a volumetric rate and offers the same normalised forms as
`InternalEnergy` through `normaliser(form)`; `DEDt` objects (and `InternalEnergy`
objects) of the same class support `+` and `-`.

The code generator uses `eos` to form the temperature column of the Jacobian via
the chain rule
$\partial \dot{x} / \partial e = (\partial \dot{x} / \partial T) / (\partial e / \partial T)$.

## eos

`#!python net.thermodynamics.eos`

**Returns**

_InternalEnergy_
: Symbolic internal energy exposing the forms below (CGS units).

| Property       | Expression                    | Units      |
| -------------- | ----------------------------- | ---------- |
| `volumetric`   | $E$                           | erg cm⁻³   |
| `specific`     | $E / \rho$                    | erg g⁻¹    |
| `per_particle` | $E / n_\mathrm{tot}$          | erg        |
| `molar`        | $N_A\, E / n_\mathrm{tot}$    | erg mol⁻¹  |

Each form is the volumetric energy divided by `InternalEnergy.normaliser(form)`
(`1`, $\rho$, $n_\mathrm{tot}$ or $n_\mathrm{tot}/N_A$), which the code
generator also uses to normalise `dE/dt`.

## EosProps

`#!python EosProps(type, **kwargs)`

Omitted parameters take the type's default (`ideal`: `gamma = 1.6666666666667`).
Validated on construction: an unknown `type`, a missing or unexpected
parameter, or an adiabatic index $\le 1$ raises `ValueError`; a non-numeric
index or a non-`dict` `gamma_map` raises `TypeError`.

| `type`                          | Parameters                                   | Volumetric energy                                                    |
| ------------------------------- | -------------------------------------------- | -------------------------------------------------------------------- |
| `ideal`                         | `gamma` (float > 1, default 1.6666666666667) | $E = \dfrac{n_\mathrm{tot}\, k_B\, T_\mathrm{gas}}{\gamma - 1}$      |
| `multi_gamma`                   | `default_gamma` (float > 1), `gamma_map` (dict[str, float > 1]) | per-species sum; `gamma_map` keyed by species name, falling back to `default_gamma` |
| `fermi_degenerate`              | —                                            | not implemented                                                      |
| `relativistic_fermi_degenerate` | —                                            | not implemented                                                      |
