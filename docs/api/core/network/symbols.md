---
tags:
    - Api
    - Network
---

# symbols

`#!python net.symbols`

`NetworkSymbols` instance owning the network's canonical SymPy symbols, density
expressions, introspection sets and symbol standardization.

**Fixed symbols** (class attributes — also usable as `NetworkSymbols.<name>` without a network)

| Attribute    | Symbol       | Meaning                                |
| ------------ | ------------ | -------------------------------------- |
| `tgas`       | `tgas`       | Gas temperature [K]                    |
| `tdust`      | `tdust`      | Dust temperature [K]                   |
| `av`         | `av`         | Visual extinction [mag]                |
| `crate`      | `crate`      | Cosmic-ray ionisation rate [s⁻¹]       |
| `chi`        | `chi`        | UV field scaling                       |
| `chi_pe`     | `chi_pe`     | Photoelectric-band field (placeholder) |
| `zd`         | `Zd`         | Grain charge                           |
| `vdisp`      | `vdisp`      | Velocity dispersion [cm s⁻¹]           |
| `photorates` | `photorates` | Photo-rate placeholder function        |

`#!python ncol(name)` → `Symbol(f"ncol_{name}")` (column density).

**Densities** (cached)

| Attribute              | Type             | Description                                    |
| ---------------------- | ---------------- | ---------------------------------------------- |
| `ndens`                | `IndexedBase`    | `nden[i]` = number density of species *i*      |
| `ntot`                 | `sympy.Expr`     | `Σ_i nden[i]`                                  |
| `rho`                  | `sympy.Expr`     | `Σ_i m_i · nden[i]`                            |
| `n_hnuc`               | `sympy.Expr`     | Hydrogen-nuclei density (`n_H_nuc` token)      |
| `element_sum(element)` | `Expr \| None`   | Nucleus density of *element*; `None` if absent |

**Introspection** (cached `frozenset`, filled after loading)

| Attribute             | Description                                          |
| --------------------- | ---------------------------------------------------- |
| `variables`           | Free symbols across rates and energy/radiation terms |
| `interp_functions`    | Names of `*interp*` functions                        |
| `undefined_functions` | Names of other undefined functions                   |

**Methods**

`#!python free_symbols(expr)` — free symbols of `expr`, excluding `nden` entries.

`#!python standardize(expr)` — replace `ntot`, `n_X`, `n_e`, `n_X_nuc`, `rc_N` and
`chi_pe` with their network expressions (names are case-insensitive).

`#!python log_summary()` — log the network's variables, interpolation and undefined
functions (called once when the network is loaded).

```python
net = Network("networks/h_photoionization/h_photo.jet")

net.symbols.ntot            # nden[0] + nden[1] + nden[2]
net.symbols.variables       # frozenset({tgas})
net.symbols.tgas            # tgas
```
