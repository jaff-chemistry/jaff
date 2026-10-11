---
tags:
    - Api
    - Code-generation
---

# get_rhs_str

`#!python get_rhs_str(idx_offset=0, use_cse=True, cse_var="cse", ode_var="f", brac_format="", def_prefix="", assignment_op="", line_end="", thermal="dedt", energy="volumetric", radiation=False, rad_order=0)`

Generates the complete RHS code block (ODE system + thermal equation). The thermal equation (`thermal`) is assigned to `ode_var[n_species]` and the radiation equations are appended after that if radiation code generation in enabled.

**Parameters**

**idx_offset** : _int, optional_
: Base index for the output array subscripts. Default `0`. Negative values use the language default (`0` for C/C++/Python/Rust, `1` for Fortran/Julia/R).

**use_cse** : _bool, optional_
: Enable joint CSE across all RHS equations. Default `True`.

**cse_var** : _str, optional_
: CSE variable prefix. Default `"cse"`.

**ode_var** : _str, optional_
: Output array name. Default `"f"`.

**brac_format** : _str, optional_
: Override 1-D bracket pair. Empty string uses the language default (`"[]"` for most languages, `"()"` for Fortran). Valid values: `"()"`, `"{}"`, `"[]"`, `"<>"`. Default `""`.

**def_prefix** : _str, optional_
: Type declaration prefix for CSE temporaries. Empty string uses the language default (`"const double "` for C/C++, `"const f64 "` for Rust, `"const Float64 "` for Julia, empty for Python/Fortran/R). Default `""`.

**assignment_op** : _str, optional_
: Assignment operator override. Empty string uses the language default (`"="` for most, `"<-"` for R). Default `""`.

**line_end** : _str, optional_
: Line terminator override. Empty string uses the language default (`";"` for C/C++/Rust, empty for Python/Fortran/Julia/R). Default `""`.

**thermal** : _str, optional_
: Thermal equation appended after the species ODEs: `"dedt"` (default, total `dE/dt`, see [get_dedt](get_dedt.md)), `"dtdt"` (total `dT/dt`, see [get_dtdt](get_dtdt.md)) or `"none"` (no thermal row). Raises `ValueError` for any other value.

**energy** : _str, optional_
: Evolved internal-energy form when `thermal="dedt"`: `"volumetric"`, `"specific"`, `"per_particle"` or `"molar"` (see [get_dedt](get_dedt.md)). Default `"volumetric"`.

**radiation** : _bool, optional_
: Include radiation moment ODEs after the energy equation. Default `False`.

**rad_order** : _int, optional_
: Radiation moment closure order (passed to `get_indexed_radodes`). Valid values: `0`, `1`, `2`, `3`. Default `0`.

**Returns**

_str_
: Complete RHS code block including the energy equation.

### Example

For a network with two species and two reactions (binding energy `deltae0 = -7.17e-12` erg), default settings produce:

```python
f[0] = -1.5e-10*tgas**(-0.5)*y[H]*y[H2]
f[1] = 1.5e-10*tgas**(-0.5)*y[H]*y[H2]
f[2] = -1.076e-21*tgas**(-0.5)*y[H]*y[H2]
```

`f[2]` is the energy equation assigned to index `n_species`.

For C++ (`lang="cxx"`) with `def_prefix="const double "`:

```python
rhs = cg.get_rhs_str(def_prefix="const double ")
```

**Output**

```cpp
const double cse0 = pow(tgas, -0.5);
f[0] = -1.5e-10*cse0*y[H]*y[H2];
f[1] = 1.5e-10*cse0*y[H]*y[H2];
f[2] = -1.076e-21*cse0*y[H]*y[H2];
```
