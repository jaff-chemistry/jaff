---
tags:
    - Api
    - Network
---

# write_table

`#!python write_table(fname, label=None, T_min=None, T_max=None, nT=64, err_tol=0.01, rate_min=1e-30, rate_max=1e100, fast_log=False, format="auto", include_all=False, verbose=False, nT_max=16385)`

General-purpose version of `to_hdf5` / `to_txt`. Tabulates reaction rate coefficients on a temperature grid and writes to a text or HDF5 file. Format is inferred from the file extension when `format="auto"`.

**Parameters**

**fname** : _str or Path_
: Output file path. When `format="auto"`, `.txt` → text, `.hdf` / `.hdf5` → HDF5.

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
: If `True`, temperature grid is uniform in `fast_log2` space rather than `log10` space. Default `False`.

**format** : _str, optional_
: `"hdf5"`, `"txt"`, or `"auto"` (inferred from extension). Default `"auto"`.

**include_all** : _bool, optional_
: If `True`, include all reactions even if they contain `NaN` values or have a constant rate. Default `False`.

**verbose** : _bool, optional_
: Print tabulation progress. Default `False`.

**nT_max** : _int or None, optional_
: Largest grid the adaptive refinement may build. When the next refinement would exceed it, refinement stops with a warning even if `err_tol` has not been reached. `None` removes the cap. Default `16385`.

**Raises**

_ValueError_
: If *format* is not one of the supported strings, or if `format="auto"` and the file extension is not recognised.
