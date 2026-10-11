---
tags:
    - Api
    - Code-generation
---

# get_jacobian_str

`#!python get_jacobian_str(thermal="none", idx_offset=0, use_cse=True, cse_var="cse", jac_var="J", matrix_format="", var_prefix="", assignment_op="", line_end="")`

Generates the complete Jacobian matrix code block. Result is cached after the first call.

**Parameters**

**thermal** : _str, optional_
: Thermal row/column included in the Jacobian: `"none"` (default), `"dedt"` (total `dE/dt`, see [get_dedt](get_dedt.md)) or `"dtdt"` (total `dT/dt`, see [get_dtdt](get_dtdt.md)). Raises `ValueError` for any other value.

**idx_offset** : _int, optional_
: Base index for row and column subscripts. Default `0`. Negative values use the language default (`0` for C/C++/Python/Rust, `1` for Fortran/Julia/R).

**use_cse** : _bool, optional_
: Apply common subexpression elimination. Default `True`.

**cse_var** : _str, optional_
: CSE variable prefix. Default `"cse"`.

**jac_var** : _str, optional_
: Jacobian matrix name. Default `"J"`.

**matrix_format** : _str, optional_
: Override 2-D bracket/separator format. Empty string uses the language default (`"]["` → `J[i][j]` for C/C++/Python/Rust, `"(,)"` → `J(i, j)` for Fortran, `"][" → J[i][j]` for Julia, `"[,]" → J[i, j]` for R). Valid values:

    | Format     | Output style |
    | ---------- | ------------ |
    | `"()"` / `"()()"`   | `J(i)(j)`    |
    | `"(,)"`    | `J(i, j)`    |
    | `"[]"` / `"[][]"`   | `J[i][j]`    |
    | `"[,]"`    | `J[i, j]`    |
    | `"{}"` / `"{}{}"`   | `J{i}{j}`    |
    | `"{,}"`    | `J{i, j}`    |
    | `"<>"` / `"<><>"`   | `J<i><j>`    |
    | `"<,>"`    | `J<i, j>`    |

    Default `""`.

**var_prefix** : _str, optional_
: Type declaration prefix for CSE temporaries. Empty string uses the language default (`"const double "` for C/C++, `"const f64 "` for Rust, `"const Float64 "` for Julia, empty for Python/Fortran/R). Default `""`.

**assignment_op** : _str, optional_
: Assignment operator override. Empty string uses the language default (`"="` for most, `"<-"` for R). Default `""`.

**line_end** : _str, optional_
: Line terminator override. Empty string uses the language default (`";"` for C/C++/Rust, empty for Python/Fortran/Julia/R). Default `""`.

**Returns**

_str_
: Jacobian code block.

### Example

For a network with two species and two reactions, default settings produce:

```python
cse0 =
J[0][0] = -1.5e-10*tgas**(-0.5)*y[H2]
J[0][1] = -1.5e-10*tgas**(-0.5)*y[H]
J[1][0] = 1.5e-10*tgas**(-0.5)*y[H2]
J[1][1] = 1.5e-10*tgas**(-0.5)*y[H]
```

Only non-zero Jacobian elements are emitted. `J[i][j]` is $\dfrac{\partial f_i}{\partial y_j}$.

For C++ (`lang="cxx"`) with `var_prefix="const double "`:

```python
jac = cg.get_jacobian_str(var_prefix="const double ")
```

**Output**

```cpp
const double cse0 = pow(tgas, -0.5);
J[0][0] = -1.5e-10*cse0*y[H2];
J[0][1] = -1.5e-10*cse0*y[H];
J[1][0] = 1.5e-10*cse0*y[H2];
J[1][1] = 1.5e-10*cse0*y[H];
```
