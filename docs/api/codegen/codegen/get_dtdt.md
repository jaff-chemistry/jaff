---
tags:
    - Api
    - Code-generation
---

# get_dtdt

`#!python get_dtdt()`

Generates code for the gas-temperature time derivative (`dT/dt`). A thin printer over [`thermodynamics.dTdt_tot`](../../core/network/thermodynamics.md), rendered as a target-language expression:

$$
\dot{T} = \frac{\dot{E}_\mathrm{tot} - \sum_i \frac{\partial E}{\partial n_i}\,\dot{n}_i}{\partial E / \partial T}
$$

The composition term accounts for reactions changing the particle numbers that share the thermal energy. It is selected in the RHS and Jacobian generators with `thermal="dtdt"` (template modifier `THERMAL dtdt`).

**Returns**

_str_
: Temperature-rate code string (single target-language expression, no assignment or line terminator).
