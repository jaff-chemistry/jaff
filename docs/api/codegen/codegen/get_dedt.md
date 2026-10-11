---
tags:
    - Api
    - Code-generation
---

# get_dedt

`#!python get_dedt(energy="volumetric")`

Generates code for the internal energy time derivative (`dE/dt`). A thin printer over [`thermodynamics.dEdt_tot`](../../core/network/thermodynamics.md): it renders `net.thermodynamics.dEdt_tot.normaliser(energy)` as a target-language expression, where `normaliser(energy)` divides the volumetric rate by the normaliser of the chosen form (`den` below). For the temperature rate see [get_dtdt](get_dtdt.md).

**Parameters**

**energy** : _str, optional_
: Evolved internal-energy form.

    - `"volumetric"` → `den = 1` (erg/cm³/s).
    - `"specific"` → `den = ρ = Σ m_i · nden[i]` (erg/g/s).
    - `"per_particle"` → `den = n_tot = Σ nden[i]` (erg/s per particle).
    - `"molar"` → `den = n_tot / N_A` (erg/mol/s).

    Default `"volumetric"`. Raises `ValueError` for any other value.

**Returns**

_str_
: Energy-equation code string (single target-language expression, no assignment or line terminator).

The normalisation is applied through the quotient rule, so for the non-volumetric forms the result also contains the term from the changing normaliser (e.g. `d(E/ρ)/dt`), not just `dE/dt / den`.
