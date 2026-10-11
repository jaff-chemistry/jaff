---
tags:
    - Api
    - Code-generation
---

# get_indexed_jacobian

`#!python get_indexed_jacobian(thermal="none", use_cse=True, cse_var="cse")`

Computes the analytical Jacobian matrix for the chemical network $\left(\dfrac{\partial f_i}{\partial y_j}\right)$ using symbolic differentiation and optional CSE.

**Parameters**

**thermal** : _str, optional_
: Thermal row included: `"none"` (default), `"dedt"` or `"dtdt"`. Raises `ValueError` for any other value.

**use_cse** : _bool, optional_
: Apply common subexpression elimination. Default `True`.

**cse_var** : _str, optional_
: CSE variable prefix. Default `"cse"`.

**Returns**

_IndexedReturn_
: `expressions` contains `IndexedValue` objects with 2D indices `[i, j]`.
