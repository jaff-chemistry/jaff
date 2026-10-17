---
tags:
    - Api
    - Network
---

# to_hdf5

`#!python to_hdf5(fname, label=None, T_min=None, T_max=None, nT=64, err_tol=0.01, rate_min=1e-30, rate_max=1e100, fast_log=False, include_all=False, verbose=False, nT_max=16385)`

Tabulates reaction rate coefficients on a temperature grid and writes to HDF5.

**Parameters**

**fname** : _str or Path_
: Output file path. If the file name doesn't end with a .hdf5 extension, the extension is forced to `.hdf5`.

**label** : _str or None, optional_
: Dataset label. Defaults to filename without extension.

**T_min** : _float or None, optional_
: Minimum temperature in Kelvin. Defaults to reaction.tmin.

**T_max** : _float or None, optional_
: Maximum temperature in Kelvin. Defaults to reaction.tmax.

**nT** : _int, optional_
: Initial number of temperature grid points. The adaptive refinement may increase this significantly. Default `64`.

**err_tol** : _float, optional_
: Maximum permitted relative interpolation error. `None` skips adaptive refinement and return exactly `nT` points. Default `0.01`.

**rate_min** : _float, optional_
: Small positive floor used in the relative-error denominator to avoid division by zero for near-zero rates. Default `1e-30`.

**rate_max** : _float, optional_
: Upper clamp applied to each rate before storing; prevents overflow in log-rate arithmetic. Default `1e100`.

**fast_log** : _bool, optional_
: Use fast logarithm approximation. Default `False`.

**include_all** : _bool, optional_
: If `True`, include all reactions in the output even if they contain `NaN` values or have a constant rate. Default `False`.

**verbose** : _bool, optional_
: Print progress. Default `False`.

**nT_max** : _int or None, optional_
: Largest grid the adaptive refinement may build. When the next refinement would exceed it, refinement stops with a warning even if `err_tol` has not been reached. `None` removes the cap. Default `16385`.
