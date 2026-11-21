# ABOUTME: .jaff network serialization (gzip JSON) and quokka-compatible
# ABOUTME: rate-coefficient table building/writing (txt and HDF5).
"""
Serialization and tabulation utilities for JAFF network files.

This module provides these public entry-points:

* ``to_jaff_file`` -- serialize a :class:`~jaff.Network` to a gzip-compressed
  JSON file with the ``jaff.network_json`` format marker.  Loading ``.jaff``
  files is not supported at the moment.
* ``get_table`` / ``write_data_table`` -- build a rate-coefficient lookup
  table over a temperature grid and optionally write it to a quokka-compatible
  text or HDF5 file.

The file format
---------------
A ``.jaff`` file is a gzip-compressed UTF-8 JSON document with the top-level
keys ``format``, ``schema_version``, ``jaff_version``, ``sympy_schema_version``,
``sympy_version``, ``label``, ``file_name``, ``species``, ``rate_symbols``,
``reactions``, ``dEdt_extra``, ``emissivity``, and ``absorption_coefficient``
(the last two ``null`` without radiation).  SymPy expressions are stored via the versioned compact encoding
in :mod:`jaff.common._sympy_json` (``SCHEMA_VERSION = 2``).
"""

from __future__ import annotations

import gzip
import json
import logging
from dataclasses import asdict
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable

import numpy as np
from sympy import (
    Basic,
    Float,
    Heaviside,
    Max,
    Min,
    Piecewise,
    __version__,
    expand_log,
    lambdify,
    log,
)
from sympy.core.function import AppliedUndef

from .. import __version__ as jaff_version
from ..common import SCHEMA_VERSION as SYMPY_SCHEMA
from ..common import fast_log2, inverse_fast_log2, is_jaff_file
from ..common import to_jsonable as sympy_to_jsonable
from ..drivers.hdf5 import HDF5
from ..types import HDF5Dict
from ._logger import JaffLogger, jaff_progress

if TYPE_CHECKING:
    from .. import Network, Reaction, Specie, Species
    from ..core.reaction import Reactions
else:
    Specie = "Specie"
    Species = "Species"
    Reaction = "Reaction"
    Network = "Network"


def to_jaff_file(filename: str | Path, net: "Network"):
    """
    Serialize a :class:`~jaff.Network` to a gzip-compressed ``.jaff`` JSON file.

    The output document uses the ``jaff.network_json`` format marker and
    ``schema_version = 1``.  SymPy expressions (rates, energy releases, …) are
    encoded with the compact :mod:`~jaff.common._sympy_json` serializer at its
    current ``SCHEMA_VERSION``.  Plain string rates (e.g. custom C/Fortran
    expressions) are wrapped in a ``{"kind": "string", "value": …}`` envelope.

    Parameters
    ----------
    filename : str or pathlib.Path
        Destination path.  If it is not already a ``.jaff`` (or ``.jaff.gz``)
        path, its suffix is replaced by ``.jaff`` (added if it has none).
    net : Network
        The reaction network to serialize.

    Raises
    ------
    ValueError
        If any reaction rate contains an undefined SymPy function (``AppliedUndef``),
        which cannot be round-tripped through the JSON schema.
    TypeError
        If a rate field is neither a :class:`sympy.Basic` expression, a plain
        :class:`str`, nor ``None``.
    """
    if isinstance(filename, str):
        filename = Path(filename)

    if not is_jaff_file(filename):
        filename = filename.with_suffix(".jaff")

    def has_undefined_functions(expr):
        """Return True if *expr* contains any undefined (``AppliedUndef``) SymPy functions."""
        if not isinstance(expr, Basic):
            return False

        if expr.atoms(AppliedUndef):
            return True

        return False

    def encode_maybe_sympy(value):
        """
        Encode a rate field to a JSON-compatible object.

        Plain strings are wrapped as ``{"kind": "string", "value": …}``.
        SymPy expressions are encoded with the compact sympy_json encoder.
        ``None`` passes through unchanged.
        """
        if isinstance(value, str):
            return {"kind": "string", "value": value}
        if isinstance(value, Basic):
            if has_undefined_functions(value):
                raise ValueError(
                    "Cannot serialize: expression contains undefined SymPy function(s)"
                )
            return sympy_to_jsonable(value, include_assumptions=False)
        if value is None:
            return None

        raise TypeError(f"Unsupported value type for serialization: {type(value)!r}")

    def jsonable(obj):
        """
        Recursively convert an object to a JSON-serializable form.

        NumPy arrays become lists; NumPy scalars become Python scalars; dicts,
        lists, and tuples are traversed recursively.  All other types are
        returned unchanged (assumed to be natively JSON-serializable).
        """
        if obj is None or isinstance(obj, (str, int, float, bool)):
            return obj
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.floating, np.integer)):
            return obj.item()
        if isinstance(obj, dict):
            return {str(k): jsonable(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [jsonable(v) for v in obj]
        return obj

    payload = {
        "format": "jaff.network_json",
        "schema_version": 1,
        "jaff_version": jaff_version,
        "sympy_schema_version": SYMPY_SCHEMA,
        "sympy_version": __version__,
        "label": net.label,
        "file_name": str(net.filepath),
        "species": [
            {
                "name": sp.name,
                "index": int(sp.index),
                "mass": float(sp.mass) if sp.mass is not None else None,
                "charge": int(sp.charge) if sp.charge is not None else None,
            }
            for sp in net.species
        ],
        "rate_symbols": [
            {
                "name": sym.name,  # type: ignore
                "assumptions": {
                    k: v
                    for k, v in (sym.assumptions0 or {}).items()
                    if isinstance(k, str) and isinstance(v, bool)
                },
            }
            for sym in sorted(
                {
                    s
                    for r in net.reactions
                    if isinstance(r.rate.coefficient, Basic)
                    for s in r.rate.coefficient.free_symbols
                }
                | set(net.thermodynamics.dEdt_extra.volumetric.free_symbols),
                key=lambda s: s.name,
            )
        ],
        "reactions": [
            {
                "reactants": [s.name for s in r.reactants],
                "products": [s.name for s in r.products],
                "rate": encode_maybe_sympy(r.rate.coefficient),
                "tmin": r.tmin,
                "tmax": r.tmax,
                "t_cutoff": r.config.t_cutoff,
                "rate_segments": [
                    {
                        "rate": encode_maybe_sympy(s.rate),
                        "tmin": s.tmin,
                        "tmax": s.tmax,
                    }
                    for s in r.rate_segments
                ],
                "dE": encode_maybe_sympy(r.dE),
                "custom_rad_rate": r.photo is not None and r.photo.custom_rad_rate,
                "original_string": r.source.line,
                "index": r.index,
                "type": r.type,
                "format": r.source.format,
                "pi_database": _photo_setting(r, "pi_database"),
                "xsecs": _encode_custom_xsec_props(r, encode_maybe_sympy),
                "spectral_yield": encode_maybe_sympy(r.rad.spectral_yield),
                "shielding": _encode_shielding(r),
            }
            for r in net.reactions
        ],
        "dEdt_extra": encode_maybe_sympy(net.thermodynamics.dEdt_extra.volumetric),
        "emissivity": encode_maybe_sympy(
            net.radiation.emissivity if net.radiation is not None else None
        ),
        "absorption_coefficient": encode_maybe_sympy(
            net.radiation.absorption_coefficient if net.radiation is not None else None
        ),
    }

    with gzip.open(filename, "wt", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)


def _photo_setting(reaction: Reaction, name: str) -> Any:
    """A photo reaction's resolved photo setting; ``None`` otherwise."""
    if reaction.photo is None:
        return None

    return getattr(reaction.photo.config, name)


def _encode_shielding(reaction: Reaction) -> dict[str, Any] | None:
    """A photo reaction's shielding as ``{"type", "props"}``, or ``None``."""
    shielding = _photo_setting(reaction, "shielding")
    if shielding is None:
        return None

    return {
        "type": shielding.type,
        "props": shielding.params.model_dump(exclude_unset=True),
    }


def _encode_custom_xsec_props(
    reaction: Reaction, encode_expr: Callable[[Any], Any]
) -> dict | None:
    """JSON form of a ``photoCrossSection`` cross section, else ``None``.

    Only custom cross sections are stored: database ones are looked up again
    when the ``.jaff`` file is loaded, but a custom one would otherwise need the
    ``.jfunc`` it came from.  *encode_expr* serializes the ``photodecay_expr``.
    """
    xsec = reaction.photo.xsec if reaction.photo is not None else None
    if xsec is None or xsec.database != "custom":
        return None

    return {
        "units": asdict(xsec.units),
        "process": asdict(xsec.process),
        "database": xsec.database,
        "photodecay_expr": encode_expr(xsec.photodecay_expr),
    }


def get_table(
    reactions: Reactions,
    logger: logging.Logger | None,
    T_min=None,
    T_max=None,
    nT=64,
    err_tol=0.01,
    rate_min=1e-30,
    rate_max=1e100,
    fast_log=False,
    verbose=False,
):
    """
    Build a temperature-indexed rate-coefficient lookup table for a reaction list.

    Before classification, ``av`` (visual extinction) is substituted with
    ``0.0`` and ``crate`` (cosmic-ray rate) with ``1.0``, so rates depending
    on them are tabulated at those values.  After that substitution, constant
    rates and rates depending **solely** on ``tgas`` (gas temperature) are
    tabulated; rates depending on any other symbol, or containing an undefined
    SymPy function, are filled with ``NaN``.

    The temperature grid is built adaptively: starting from ``nT`` points the
    algorithm doubles the grid by inserting midpoints until the maximum relative
    interpolation error across all reactions falls below ``err_tol``.  All
    arithmetic is performed in log-rate space to avoid catastrophic cancellation
    and overflow.

    When ``fast_log=True`` the grid is uniform in :func:`~jaff.common.fast_log2`
    space (rather than in ``log10(T)`` space) so that it can be directly indexed
    by the fast-log approximation used in quokka's runtime interpolation kernel.

    Parameters
    ----------
    reactions : Reactions
        Iterable of :class:`~jaff.core.Reaction` objects for which rates are to
        be tabulated.
    logger : logging.Logger or None
        Logger instance used for progress/diagnostic messages.  When ``None``
        a default :class:`~jaff.io._logger.JaffLogger` is constructed.
    T_min : float or None, optional
        Minimum temperature in K.  The smallest ``reaction.tmin`` when not
        supplied.
    T_max : float or None, optional
        Maximum temperature in K.  The largest ``reaction.tmax`` when not
        supplied.
    nT : int, optional
        Initial number of temperature grid points (default ``64``).  The
        adaptive refinement may increase this significantly.
    err_tol : float or None, optional
        Maximum permitted relative interpolation error (default ``0.01``).
        Pass ``None`` to skip adaptive refinement and return exactly ``nT``
        points.
    rate_min : float, optional
        Small positive floor used in the relative-error denominator to avoid
        division by zero for near-zero rates (default ``1e-30``).
    rate_max : float, optional
        Upper clamp applied to each rate before storing; prevents overflow in
        log-rate arithmetic (default ``1e100``).
    fast_log : bool, optional
        If ``True``, the temperature grid is uniform in
        :func:`~jaff.common.fast_log2` space instead of ``log10`` space
        (default ``False``).
    verbose : bool, optional
        If ``True``, log the current grid size, worst-case error, and the
        offending reaction and temperature after each refinement step
        (default ``False``).

    Returns
    -------
    temp : numpy.ndarray, shape (nTemp,)
        Temperature grid in K.
    coeff : numpy.ndarray, shape (nReactions, nTemp)
        Rate coefficients in cm^3 s^-1 (or the appropriate units for the
        reaction type).  Entries for reactions that cannot be tabulated are
        ``NaN``.

    Raises
    ------
    ValueError
        If ``T_min`` (``T_max``) is not supplied and no reaction has a finite
        ``tmin`` (``tmax``), e.g. *reactions* is empty; or if
        ``T_min >= T_max``.

    Notes
    -----
    The log-expand trick (``expand_log(log(r))``) is applied before
    :func:`sympy.lambdify` to keep intermediate values in a numerically safe
    range.  Each lambdified rate is evaluated with a single vectorised NumPy
    call over the grid; rates containing ``Piecewise``, ``Heaviside``, ``Min``
    or ``Max`` (and any rate whose vectorised result has the wrong shape) are
    evaluated point by point in a Python loop instead.
    """

    if logger is None:
        logger = JaffLogger().get_logger()

    if T_min is None:
        T_min = min((r.tmin for r in reactions if r.tmin is not None), default=None)
    if T_max is None:
        T_max = max((r.tmax for r in reactions if r.tmax is not None), default=None)
    if T_min is None or T_max is None or not np.isfinite([T_min, T_max]).all():
        raise ValueError(
            "could not determine T_min or T_max from "
            "reaction list; set T_min and T_max manually"
        )
    if T_min >= T_max:
        raise ValueError(f"T_min ({T_min}) must be below T_max ({T_max})")

    react_sympy = [r.get_sympy() for r in reactions]

    # Local import: a module-level one would cycle via jaff.core.network.
    from ..core import NetworkSymbols

    trivial_subs = {NetworkSymbols.av: Float(0.0), NetworkSymbols.crate: Float(1.0)}
    react_subst = [r.xreplace(trivial_subs) for r in react_sympy]

    tgas = NetworkSymbols.tgas
    react_func = []
    react_vectorizable = []
    branchy = (Piecewise, Heaviside, Min, Max)
    for r in jaff_progress.track(react_subst, description="Compiling rate functions"):
        free_symbols = r.free_symbols
        if len(free_symbols) == 0:
            react_func.append(np.log(float(r)))
            react_vectorizable.append(False)
        elif (
            (len(free_symbols) > 1)
            or (tgas not in free_symbols)
            or bool(r.atoms(AppliedUndef))
        ):
            react_func.append(None)
            react_vectorizable.append(False)
        else:
            # log-expand before lambdify to avoid overflow; exponentiate at the end
            logr = expand_log(log(r))
            react_func.append(lambdify(tgas, logr, "numpy"))
            react_vectorizable.append(not logr.has(*branchy))

    def eval_rate(i, f, grid):
        """Evaluate lambdified rate *f* (index *i*) over a temperature array.

        Uses a single vectorised call for analytic rates and falls back to a
        scalar Python loop for branching rates.  The vectorised result is
        verified against the scalar one the first time it is used so that a
        mis-broadcasting lambdified function silently demotes itself to scalar.
        """
        if not react_vectorizable[i]:
            return np.array([f(t) for t in grid])
        vec = np.asarray(f(grid), dtype=float)
        if vec.shape != grid.shape:
            react_vectorizable[i] = False
            return np.array([f(t) for t in grid])
        return vec

    nTemp = nT
    if not fast_log:
        temp = np.logspace(np.log10(T_min), np.log10(T_max), nTemp)
    else:
        log_temp_min = fast_log2(T_min)
        log_temp_max = fast_log2(T_max)
        log_temp = np.linspace(log_temp_min, log_temp_max, nTemp)
        temp = inverse_fast_log2(log_temp)
    log_rates = np.zeros((len(react_func), nTemp))
    for i, f in jaff_progress.track(
        list(enumerate(react_func)), description="Evaluating rates"
    ):
        if isinstance(f, float):
            log_rates[i, :] = f
        elif f is None:
            log_rates[i, :] = np.nan
        else:
            f_eval = eval_rate(i, f, temp)
            log_rates[i, :] = np.clip(f_eval, a_min=None, a_max=np.log(rate_max))

    # Adaptive growth of the table
    if err_tol is not None:
        with jaff_progress.indeterminate("Refining rate table"):
            while True:
                # Compute estimates at half-way points
                nTemp = 2 * nTemp - 1
                temp_grow = np.zeros(nTemp)
                temp_grow[::2] = temp
                if not fast_log:
                    temp_grow[1::2] = np.sqrt(temp[1:] * temp[:-1])
                else:
                    log_temp_lo = fast_log2(temp[:-1])
                    log_temp_hi = fast_log2(temp[1:])
                    temp_grow[1::2] = inverse_fast_log2(0.5 * (log_temp_lo + log_temp_hi))
                log_rates_grow = np.zeros((len(react_func), nTemp))
                log_rates_grow[:, ::2] = log_rates
                log_rates_approx = np.zeros((len(react_func), (nTemp - 1) // 2))
                for i, f in enumerate(react_func):
                    if isinstance(f, float):
                        log_rates_grow[i, 1::2] = f
                        log_rates_approx[i, :] = f
                    elif f is None:
                        log_rates_grow[i, 1::2] = np.nan
                        log_rates_approx[i, :] = np.nan
                    else:
                        f_eval = eval_rate(i, f, temp_grow[1::2])
                        log_rates_grow[i, 1::2] = np.clip(
                            f_eval, a_min=None, a_max=np.log(rate_max)
                        )
                        log_rates_approx[i, :] = 0.5 * (
                            log_rates_grow[i, :-1:2] + log_rates_grow[i, 2::2]
                        )

                temp = temp_grow
                log_rates = log_rates_grow

                rel_err = np.abs(
                    (np.exp(log_rates_approx) - np.exp(log_rates[:, 1::2]))
                    / (np.exp(log_rates[:, 1::2]) + rate_min)
                )
                max_err = np.nanmax(rel_err)

                if verbose:
                    idx_max = np.unravel_index(np.nanargmax(rel_err), rel_err.shape)
                    logger.info(
                        f"nTemp = {nTemp}, max_err = {max_err} in reaction "
                        f"{reactions[idx_max[0]].get_verbatim()} at T = {temp[idx_max[1]]}"
                    )

                if max_err < err_tol:
                    break

    return temp, np.exp(log_rates)


def write_data_table(
    reactions: Reactions,
    logger: logging.Logger | None,
    fname: str | Path,
    label: str | None = None,
    T_min=None,
    T_max=None,
    nT=64,
    err_tol=0.01,
    rate_min=1e-30,
    rate_max=1e100,
    fast_log=False,
    format="auto",
    include_all=False,
    verbose=False,
):
    """
    Build a rate-coefficient table and write it to a text or HDF5 file.

    Calls :func:`get_table` with the given sampling parameters and then
    serializes the result in a format understood by the quokka hydro code.
    When ``include_all=False`` (the default), reactions with any ``NaN``
    entries or a completely flat rate curve are excluded from the output.

    Supported output formats
    ------------------------
    ``txt``
        A plain-text file with a header block describing the reactions followed
        by the data in the quokka 1-D table layout:

        ```text
        1              # table dimensionality
        <nReact>       # number of outputs per entry
        2 or 3         # axis spacing (2 = log10, 3 = fast_log)
        <nTemp>        # number of temperature points
        <T_min> <T_max>
        <coeff row 0>
        ...
        ```

    ``hdf5``
        An HDF5 file structured as a quokka ``reaction_coeff`` group with the
        group attributes ``input_names``, ``input_units``, ``xlo``, ``xhi``,
        and ``spacing``.

    Parameters
    ----------
    reactions : Reactions
        Iterable of :class:`~jaff.core.Reaction` objects.
    logger : logging.Logger or None
        Logger instance.  ``None`` creates a default :class:`JaffLogger`.
    fname : str or pathlib.Path
        Destination path.  The output format is deduced from the file
        extension when ``format="auto"``.
    label : str or None, optional
        Human-readable network label written into the file header.  Defaults
        to the file stem when ``None``.
    T_min : float or None, optional
        Minimum temperature in K.  See :func:`get_table`.
    T_max : float or None, optional
        Maximum temperature in K.  See :func:`get_table`.
    nT : int, optional
        Initial number of temperature grid points (default ``64``).
    err_tol : float or None, optional
        Adaptive refinement tolerance (default ``0.01``).
    rate_min : float, optional
        Error-denominator floor (default ``1e-30``).
    rate_max : float, optional
        Rate upper clamp (default ``1e100``).
    fast_log : bool, optional
        Use :func:`~jaff.common.fast_log2` spacing for the temperature axis
        (default ``False``).
    format : {"auto", "txt", "hdf5"}, optional
        Output format.  ``"auto"`` (the default) infers the format from the
        file extension (``.txt`` → text, ``.hdf`` / ``.hdf5`` → HDF5).
    include_all : bool, optional
        If ``True``, include all reactions in the output even if they contain
        ``NaN`` values or have a constant rate.  Default ``False``.
    verbose : bool, optional
        Forward to :func:`get_table` for per-refinement-step logging.

    Raises
    ------
    ValueError
        If *format* is not one of the supported strings, or if ``format="auto"``
        and the file extension is not recognised.
    """
    if logger is None:
        logger = JaffLogger().get_logger()

    if isinstance(fname, str):
        fname = Path(fname)

    supported_formats = ["auto", "hdf5", "txt"]

    if format not in supported_formats:
        raise ValueError(
            f"Unknown output format {format}\n"
            f"Supported output formats are {', '.join(supported_formats)}"
        )

    if format in supported_formats and format != "auto":
        out_type = supported_formats[supported_formats.index(format)]
    elif format == "auto":
        # Map each recognised extension to its canonical output format.
        ext_to_type = {".hdf": "hdf5", ".hdf5": "hdf5", ".txt": "txt"}
        if fname.suffix in ext_to_type:
            out_type = ext_to_type[fname.suffix]
        else:
            raise ValueError(
                f"Cannot deduce output type from extension for file: {fname}\n"
                f"Supported extensions are .txt, .hdf5 and .hdf"
            )

    if label is None:
        label = fname.stem

    temp, coef = get_table(
        reactions=reactions,
        logger=logger,
        T_min=T_min,
        T_max=T_max,
        nT=nT,
        err_tol=err_tol,
        rate_min=rate_min,
        rate_max=rate_max,
        fast_log=fast_log,
        verbose=verbose,
    )

    if include_all:
        react_list = list(range(len(coef)))
    else:
        react_list = []
        for i, c in enumerate(coef):
            if np.sum(np.isnan(c)) > 0 or np.amax(c) - np.amin(c) == 0.0:
                continue
            react_list.append(i)
    coef = coef[react_list]

    rtype = []
    reactants = []
    products = []
    for i in react_list:
        if reactions[i].type == "unknown":
            rtype.append("2_body")
        else:
            rtype.append(reactions[i].type)
        reactants_ = {}
        for r in reactions[i].reactants:
            if r.name in reactants_.keys():
                reactants_[r.name] += 1
            else:
                reactants_[r.name] = 1
        reactants.append(reactants_)
        products_ = {}
        for p in reactions[i].products:
            if p.name in products_.keys():
                products_[p.name] += 1
            else:
                products_[p.name] = 1
        products.append(products_)

    def to_text():
        """Write the table in quokka plain-text 1-D lookup format."""
        with open(fname, "w", encoding="utf-8") as fp:
            fp.write("# JAFF auto-generated rate coefficient table\n")
            fp.write(f"# Network name: {label}\n")
            fp.write("# Reactions included\n")
            fp.write("#   (reactants) (products) (reaction type)\n")
            for rt, r, p in zip(rtype, reactants, products):
                fp.write(f"#   {r} {p} {rt}\n")

            # Write data in quokka table format
            fp.write("1\n")  # Table is 1d
            fp.write(f"{len(coef)}\n")  # N outputs per table entry
            if fast_log:
                fp.write("3\n")  # Table is uniform in fast_log
            else:
                fp.write("2\n")  # Table is uniform in log
            fp.write(f"{len(temp)}\n")  # Number of temperature entries
            fp.write(f"{temp[0]} {temp[-1]}\n")  # Min/max temperature

            # Now write the data
            for c in coef:
                for c_ in c:
                    fp.write(f"{c_} ")
                fp.write("\n")

    def to_hdf5():
        """Write the table as a quokka-compatible HDF5 ``reaction_coeff`` group."""
        output_names = []
        output_units = []
        for i, rt, r, p in zip(range(len(rtype)), rtype, reactants, products):
            output_names.append(f"{rt} rate coefficient: {r} --> {p}")
            output_units.append("cm^3 s^-1")

        hdfdict = HDF5Dict(
            {
                "reaction_coeff": {
                    "_attrs": {
                        "input_names": ["temperature"],
                        "input_units": ["K"],
                        "xlo": np.array([temp[0]]),
                        "xhi": np.array([temp[-1]]),
                        "spacing": ["fast_log"] if fast_log else ["log"],
                    },
                    "output_names": {
                        "_data": output_names,
                        "_dtype": "s",
                        "_kind": "linear",
                    },
                    "output_units": {
                        "_data": output_units,
                        "_dtype": "s",
                        "_kind": "linear",
                    },
                    "data": {
                        "_data": coef,
                        "_kind": "linear",
                    },
                }
            }
        )

        HDF5().from_dict(fname, hdfdict, mode="w")

    # Write output in appropriate format
    if out_type == "txt":
        to_text()
    elif out_type == "hdf5":
        to_hdf5()
