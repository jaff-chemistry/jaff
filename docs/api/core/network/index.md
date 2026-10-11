---
tags:
    - Api
    - Network
---

# Network

`jaff.core.network.Network`

The `Network` class is the most important class in JAFF. It reads a reaction network file, auto-detects its format, validates mass and charge conservation, and assembles the full species and reaction catalogues along with stoichiometry matrices. It also handles optional radiation transport, photochemistry cross-sections, and auxiliary function files.

## Constructor

`#!python Network(fname, config=None, errors=False, label=None, funcfile=True, duplicate_policy=None, expand_nuclei=True, radiation_props=None, dust_props=None, use_proxy_photoreaction=False, eos_props=None)`

**Parameters**

**fname** : _str or Path_
: Path to a network file, or the name of a built-in network (a sub-directory of `networks/` containing a single `.jet` file). A built-in network name wins over a same-named path on disk. Supported formats: KIDA, UDFA, PRIZMO, KROME, UCLCHEM, a combination of the above and the `.jaff` file (Refer to [to_jaff](to_jaff.md) for more details).

**config** : _str, Path, or None, optional_
: Path to a TOML configuration file. When `None` (default), JAFF looks for `jaff.toml` in the network file's directory.

**errors** : _bool, optional_
: Exit on validation errors. Default `False`.

**label** : _str or None, optional_
: Network identifier. Defaults to the file stem.

**funcfile** : _bool, str, or Path, optional_
: Path to .jfunc auxiliary functions file. `True` (default) scans the network directory; `False` skips.

**duplicate_policy** : _str or None, optional_
: Resolve duplicate rate coefficients (same reaction, mechanism, and temperature range): `preserve-first`, `preserve-last`, or `error`. When `None` (default), the network's `jaff.toml` value is used, falling back to `preserve-first`.

**expand_nuclei** : _bool, optional_
: When `True` (default), an `n_<element>_nuc` symbol in rate expressions (e.g. `n_H_nuc`, `n_He_nuc`) is expanded to the element-nucleus density sum (weighted by atom count over all species bearing that element). When `False`, it is left as a free symbol `n<element>_nuc` (e.g. `nh_nuc`) instead. Only affects `n_<element>_nuc`; plain `n_X` always resolves to species `X`.

**radiation_props** : _RadiationProps or None, optional_
: Radiation-field configuration (bands, spectral index, mode, speed of light, background field). `None` (default) disables radiation transport.

**dust_props** : _DustProps or None, optional_
: Dust-module configuration (Rv, radiation reductions, photoelectric band edges). `None` (default) disables the dust module.

**use_proxy_photoreaction** : _bool, optional_
: Use proxy photo-reactions when computing cross-sections instead of bypassing them. Default `False`.

**eos_props** : _EosProps or None, optional_
: Equation-of-state configuration used by [`thermodynamics.eos`](thermodynamics.md). `None` (default) uses an ideal gas, `EosProps("ideal")` with `gamma = 1.6666666666667`.

**Raises**

_FileNotFoundError_
: If `fname` does not exist.

## Attributes

| Attribute         | Type                 | Description                                                                                                         |
| ----------------- | -------------------- | ------------------------------------------------------------------------------------------------------------------- |
| `label`           | `str`                | Human-readable network identifier; defaults to the source file stem                                                 |
| `filename`        | `Path`               | Resolved absolute path to the source network file                                                                   |
| `spec`            | `NetworkSpec`        | Normalized construction parameters (resolved `fname`, parsed `config` dict, `funcfile`, `aux_funcs`, ...)           |
| `species`         | `Species`            | Ordered catalogue of the network's core (real) species; special pseudo-species (`_PHOTON`, `_CR`, ...) are excluded |
| `reactions`       | `Reactions`          | Ordered catalogue of all reactions in the network                                                                   |
| `elements`        | `Elements`           | Element catalogue derived from all species; used for composition matrices                                           |
| `reactant_matrix` | `ndarray`            | Shape (n_reactions, n_species) stoichiometry matrix for reactants                                                   |
| `product_matrix`  | `ndarray`            | Shape (n_reactions, n_species) stoichiometry matrix for products                                                    |
| `mass_dict`       | `dict`               | Mapping from element symbol to mass properties, used for conservation checks                                        |
| `thermodynamics`  | `Thermodynamics`     | Thermal equations: `eos`, `dEdt_chemical`, `dEdt_extra`, `dEdt_tot`, `dTdt_chemical`, `dTdt_extra`, `dTdt_tot`; see [thermodynamics](thermodynamics.md) |
| `dRad_dt_extra`   | `sympy.Basic`        | Extra radiation moment source terms from `@function` definitions                                                    |
| `radiation`       | `Radiation or None`  | Radiation field object; `None` when no radiation bands are specified                                                |
| `symbols`         | `NetworkSymbols`     | Canonical symbols, densities (`ndens`, `ntot`, `rho`, `n_hnuc`), introspection sets and `standardize`; see [symbols](symbols.md) |
