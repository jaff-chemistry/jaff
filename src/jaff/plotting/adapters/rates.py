# ABOUTME: Rate-coefficient adapter: reaction / SymPy expr / (x, y) -> series dict
# ABOUTME: Evaluates on a log-spaced temperature grid; unevaluable rates skipped

"""Rate coefficients as plotting series.

Duck-typed: a "reaction" is anything with a ``rate.coefficient`` (and optionally
``tmin`` / ``tmax`` / ``get_latex``), so this module never imports
:mod:`jaff.core`.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from sympy import Basic, Symbol, lambdify

from ...io import JaffLogger
from ._common import label_of as _label_of

#: Fallback temperature bounds (K) when neither caller nor reaction gives one.
DEFAULT_TMIN: float = 2.73
DEFAULT_TMAX: float = 1e6


def _first(*values: float | None) -> float:
    """First non-``None`` value as a float."""
    return float(next(v for v in values if v is not None))


def _usable(value: Any) -> float | None:
    """A reaction's own bound if it is finite and positive, else ``None``."""
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) and v > 0 else None


def _evaluate(
    expr: Basic,
    var: str,
    t0: float,
    t1: float,
    npoints: int,
    label: str,
    quiet: bool = False,
) -> dict | None:
    """Evaluate *expr* on a log grid over ``[t0, t1]``; ``None`` if impossible.

    Failures are logged as warnings, or at debug level when *quiet*.
    """
    if not 0 < t0 < t1:
        raise ValueError(
            f"Temperature bounds must satisfy 0 < tmin < tmax; got {t0}, {t1}"
        )
    logger = JaffLogger().get_logger()
    log = logger.debug if quiet else logger.warning
    if expr is None:
        log(f"Rate {label!r} is None. Skipping.")
        return None
    grid = np.logspace(np.log10(t0), np.log10(t1), npoints)
    try:
        f = lambdify(var, expr, "numpy")
        y = np.array([float(f(t)) for t in grid])
    except (NameError, TypeError, ValueError, AttributeError, KeyError) as exc:
        log(f"Cannot evaluate rate {label!r}: {exc}. Skipping.")
        return None
    return {"x": grid, "y": y, "label": label}


def rate_series(
    item: Any,
    *,
    tmin: float | None = None,
    tmax: float | None = None,
    var: str = "tgas",
    npoints: int = 200,
    label: str | None = None,
    quiet: bool = False,
) -> dict | None:
    """Turn one rate input into a series dict ``{"x", "y", "label"}``.

    Parameters
    ----------
    item
        A reaction-like object (has ``rate.coefficient``), a SymPy expression, or an
        ``(x, y)`` tuple of arrays.
    tmin, tmax : float or None
        Temperature range (K).  For reactions falls back to the reaction's
        own bounds (ignored when ``None``, non-finite or <= 0), then to
        2.73 K / 1e6 K.  Required for bare expressions.
    var : str
        Temperature symbol name substituted in reaction rates.
    npoints : int
        Number of log-spaced samples.
    label : str or None
        Legend label; defaults to the reaction's LaTeX equation (or ``str``).
    quiet : bool
        Log an unevaluable rate at debug level instead of as a warning (for
        callers that report a summary).

    Returns
    -------
    dict or None
        ``None`` if the rate cannot be evaluated numerically, e.g. photo or
        cosmic-ray rates with free symbols; this is logged as a warning, or
        at debug level when *quiet*.

    Raises
    ------
    ValueError
        A bare SymPy expression without ``tmin`` and ``tmax``, or (for a
        reaction or SymPy expression) temperature bounds that do not satisfy
        ``0 < tmin < tmax``.
    TypeError
        An unsupported *item* type.
    """
    if isinstance(item, tuple) and len(item) == 2:
        x, y = item
        return {
            "x": np.asarray(x, dtype=float),
            "y": np.asarray(y, dtype=float),
            "label": label or "",
        }

    if hasattr(getattr(item, "rate", None), "coefficient"):
        t0 = _first(tmin, _usable(getattr(item, "tmin", None)), DEFAULT_TMIN)
        t1 = _first(tmax, _usable(getattr(item, "tmax", None)), DEFAULT_TMAX)
        lab = label if label is not None else _label_of(item)
        return _evaluate(item.rate.coefficient, var, t0, t1, npoints, lab, quiet)

    if isinstance(item, Basic):
        if tmin is None or tmax is None:
            raise ValueError("tmin and tmax are required for bare SymPy expressions.")
        syms = list(item.free_symbols)
        sym = syms[0].name if len(syms) == 1 and isinstance(syms[0], Symbol) else var
        lab = label if label is not None else str(item)
        return _evaluate(item, sym, tmin, tmax, npoints, lab, quiet)

    raise TypeError(
        f"Cannot plot rate input of type {type(item).__name__!r}; expected a "
        "reaction-like object, a SymPy expression, or an (x, y) tuple."
    )
