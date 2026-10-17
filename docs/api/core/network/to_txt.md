---
tags:
    - Api
    - Network
---

# to_txt

`#!python to_txt(fname, label=None, T_min=None, T_max=None, nT=64, err_tol=0.01, rate_min=1e-30, rate_max=1e100, fast_log=False, include_all=False, verbose=False, nT_max=16385)`

Tabulates reaction rate coefficients on a temperature grid and writes to a quokka-compatible plain-text file. Parameters are identical to `to_hdf5`.

The text layout is:

```
1              # table dimensionality
<nReact>       # number of outputs per entry
2 or 3         # axis spacing (2 = log10, 3 = fast_log)
<nTemp>        # number of temperature points
<T_min> <T_max>
<coeff row 0>
...
```

**Parameters**

**fname** : _str or Path_
: Output file path. The `.txt` extension is recommended.

**label** : _str or None, optional_
: Dataset label written into the file header. Defaults to `self.label`.

**T_min** : _float or None, optional_
: Minimum temperature in Kelvin. Inferred from `reaction.tmin` values when not supplied.

**T_max** : _float or None, optional_
: Maximum temperature in Kelvin. Inferred from `reaction.tmax` values when not supplied.

**nT** : _int, optional_
: Initial number of temperature grid points. Adaptive refinement may increase this significantly. Default `64`.

**err_tol** : _float or None, optional_
: Maximum permitted relative interpolation error. `None` skips adaptive refinement and returns exactly `nT` points. Default `0.01`.

**rate_min** : _float, optional_
: Small positive floor used in the relative-error denominator to avoid division by zero for near-zero rates. Default `1e-30`.

**rate_max** : _float, optional_
: Upper clamp applied to each rate before storing; prevents overflow in log-rate arithmetic. Default `1e100`.

**fast_log** : _bool, optional_
: If `True`, temperature grid is uniform in `fast_log2` space (axis spacing code `3`) rather than `log10` space (code `2`). Default `False`.

**include_all** : _bool, optional_
: If `True`, include all reactions even if they contain `NaN` values or have a constant rate. Default `False`.

**verbose** : _bool, optional_
: Print tabulation progress. Default `False`.

**nT_max** : _int or None, optional_
: Largest grid the adaptive refinement may build. When the next refinement would exceed it, refinement stops with a warning even if `err_tol` has not been reached. `None` removes the cap. Default `16385`.
