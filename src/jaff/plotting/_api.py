"""High-level plotting functions for rates and cross sections.

These are the preferred entry points for plotting.  They accept a single item
or a list, build a tidy long frame, and delegate rendering to
:meth:`jaff.plotting.plotter.Plotter.render_series`.

The functions are *domain-agnostic* by duck typing: a "reaction" is anything
exposing the attributes used here (``rate`` / ``tmin`` / ``tmax`` /
``get_latex`` for rates; ``xsecs_dict`` / ``band_xsecs`` / ``get_latex`` for
cross sections).  Nothing in this module imports :mod:`jaff.core`, so plotting
stays a leaf dependency and can also plot bare SymPy expressions and arrays.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from sympy import Basic, lambdify

from ..io import JaffLogger
from . import _frames, _units
from ._frames import _sampled_xsecs  # noqa: F401  (re-exported)
from .plotter import Plotter

# Valid photo cross-section process keys.
_XSEC_PROCESSES: tuple[str, ...] = ("photo_absorption", "photodecay")

# Default temperature bounds for rate coefficients when none are supplied.
_DEFAULT_TMIN: float = 2.73
_DEFAULT_TMAX: float = 1e6


def _as_list(obj: Any) -> list[Any]:
    """Wrap a single item in a list; pass a list/tuple through as a list."""
    if isinstance(obj, (list, tuple)):
        return list(obj)
    return [obj]


def _rate_series(
    item: Any,
    tmin: float | None,
    tmax: float | None,
    var: str,
    npoints: int,
    label: str | None,
    logger: Any,
) -> tuple[np.ndarray, np.ndarray, str] | None:
    """Coerce one rate input to an ``(x, y, label)`` triple, or ``None``.

    Accepts (duck-typed, dispatched per item):

    - a reaction-like object (has ``rate``): evaluated over its temperature
      range on a log grid; label defaults to its LaTeX equation;
    - a SymPy expression: evaluated over ``[tmin, tmax]`` (required);
    - an ``(x, y)`` pair of arrays: used directly.
    """
    # (x, y) array pair -- used verbatim.
    if isinstance(item, tuple) and len(item) == 2:
        x, y = item
        return np.asarray(x, dtype=float), np.asarray(y, dtype=float), label or ""

    # Reaction-like: has a symbolic `rate` (and usually tmin/tmax/get_latex).
    if hasattr(item, "rate"):
        t0 = (
            tmin
            if tmin is not None
            else _coalesce(getattr(item, "tmin", None), _DEFAULT_TMIN)
        )
        t1 = (
            tmax
            if tmax is not None
            else _coalesce(getattr(item, "tmax", None), _DEFAULT_TMAX)
        )
        expr = item.rate
        lab = label if label is not None else _latex_or_str(item)
        return _eval_expr(expr, var, t0, t1, npoints, lab, logger)

    # Bare SymPy expression.
    if isinstance(item, Basic):
        if tmin is None or tmax is None:
            raise ValueError(
                "tmin and tmax are required when plotting bare SymPy expressions."
            )
        syms = list(item.free_symbols)
        sym = syms[0].name if len(syms) == 1 else var
        lab = label if label is not None else str(item)
        return _eval_expr(item, sym, tmin, tmax, npoints, lab, logger)

    raise TypeError(
        f"Cannot plot rate input of type {type(item).__name__!r}; expected a "
        "reaction-like object, a SymPy expression, or an (x, y) array pair."
    )


def _eval_expr(
    expr: Basic,
    var: str,
    t0: float,
    t1: float,
    npoints: int,
    label: str,
    logger: Any,
) -> tuple[np.ndarray, np.ndarray, str] | None:
    """Evaluate *expr* over a log-spaced grid in ``[t0, t1]``.

    Returns ``None`` (with a warning) if the expression cannot be evaluated
    numerically -- e.g. a photo-reaction rate that still carries the symbolic
    radiation-density variable.
    """
    grid = np.logspace(np.log10(t0), np.log10(t1), npoints)
    try:
        f = lambdify(var, expr, "numpy")
        y = np.array([float(f(t)) for t in grid])
    except (NameError, TypeError, ValueError) as exc:
        logger.warning(f"Cannot evaluate rate {label!r}: {exc}. Skipping.")
        return None
    return grid, y, label


def _coalesce(value: Any, default: float) -> float:
    """Return *value* as a float, or *default* when *value* is ``None``."""
    return default if value is None else float(value)


def _latex_or_str(item: Any) -> str:
    """LaTeX label for a reaction-like object, falling back to ``str``."""
    return item.get_latex() if hasattr(item, "get_latex") else str(item)


def plot_rates(
    rates: Any,
    *,
    tmin: float | None = None,
    tmax: float | None = None,
    var: str = "tgas",
    npoints: int = 100,
    labels: list[str] | None = None,
    palette: list[str] | None = None,
    xlabel: str = "Temperature (K)",
    ylabel: str = r"Rate coefficient $k$",
    xscale: str = "log",
    yscale: str = "log",
    shade: bool | float = False,
    title: str = "",
    fig: Any = None,
    ax: Any = None,
    grid: bool = True,
    show: bool = True,
    save: bool = False,
    filename: str = "rates.png",
) -> tuple[Any, Any] | None:
    """Plot one or more reaction rate coefficients on shared axes.

    Parameters
    ----------
    rates
        A single item or a list of: reaction-like objects (exposing ``rate``),
        SymPy expressions, or ``(x, y)`` array pairs.  All are drawn on the
        same axes with a legend.
    tmin, tmax : float or None, optional
        Temperature range (K).  For reactions, falls back to each reaction's
        own bounds, then to ``2.73`` / ``1e6``.  Required for bare SymPy
        expressions.
    var : str, optional
        Symbol name to substitute when evaluating reaction rates, by default
        ``"tgas"``.
    npoints : int, optional
        Number of log-spaced sample points, by default ``100``.
    labels : list[str] or None, optional
        Legend labels aligned to *rates*.  Defaults to each reaction's LaTeX
        equation (or ``str``).
    palette : list[str] or None, optional
        Colour cycle override.
    shade : bool or float, optional
        Shade under each curve.
    title, xlabel, ylabel, xscale, yscale, fig, ax, grid, show, save, filename
        Standard rendering controls.

    Returns
    -------
    tuple[matplotlib.figure.Figure, matplotlib.axes.Axes] or None
        ``None`` when no input could be evaluated.

    Notes
    -----
    Photo-reaction rates carry a symbolic radiation-density variable and cannot
    be evaluated as a function of temperature; such inputs are skipped with a
    warning.
    """
    logger = JaffLogger().get_logger()
    items = _as_list(rates)
    label_list = labels if labels is not None else [None] * len(items)
    if len(label_list) != len(items):
        raise ValueError("labels must have the same length as rates.")

    series = []
    for item, label in zip(items, label_list):
        result = _rate_series(item, tmin, tmax, var, npoints, label, logger)
        if result is not None:
            series.append(result)

    if not series:
        logger.info("plot_rates: nothing to plot.")
        return None

    df = _frames.long_frame(series, "x", "y", "series")
    return Plotter(palette=palette).render_series(
        df,
        x_col="x",
        y_col="y",
        group_col="series",
        xlabel=xlabel,
        ylabel=ylabel,
        log_x=(xscale == "log"),
        log_y=(yscale == "log"),
        dynamic_x=False,
        trim=False,
        shade=shade,
        title=title,
        fig=fig,
        ax=ax,
        grid=grid,
        show=show,
        save=save,
        filename=filename,
    )


def _normalize_processes(processes: str | list[str] | None) -> list[str]:
    """Normalise a process selection to a validated list of keys."""
    if processes is None or processes == "all":
        return list(_XSEC_PROCESSES)
    procs = [processes] if isinstance(processes, str) else list(processes)
    invalid = [p for p in procs if p not in _XSEC_PROCESSES]
    if invalid:
        raise KeyError(
            f"Invalid cross-section(s) {invalid}. Supported: {', '.join(_XSEC_PROCESSES)}"
        )
    return procs


def plot_xsecs(
    reactions: Any,
    *,
    processes: str | list[str] | None = "all",
    layout: str = "overlay",
    energy_unit: str = "eV",
    xsec_unit: str = "Mb",
    energy_log: bool = True,
    xsecs_log: bool = True,
    trim: bool = True,
    shade: bool | float = False,
    show_bands: bool = False,
    palette: list[str] | None = None,
    title: str | None = None,
    fig: Any = None,
    ax: Any = None,
    grid: bool = True,
    show: bool = True,
    save: bool = False,
    filename: str = "",
) -> tuple[Any, Any] | None:
    """Plot photo cross sections for one or more reactions on shared axes.

    Parameters
    ----------
    reactions
        A single reaction-like object (exposing ``xsecs_dict`` and
        ``band_xsecs``) or a list of them.  Their cross sections are overlaid.
    processes : str | list[str] | None, optional
        Which processes to draw: ``"all"`` (default), a single key, or a list.
        Valid keys: ``"photo_absorption"``, ``"photodecay"``.
    layout : str, optional
        ``"overlay"`` (default) or ``"subplots"`` (one panel per curve).
    energy_unit, xsec_unit : str, optional
        Axis units; defaults ``"eV"`` and ``"Mb"``.
    energy_log, xsecs_log : bool, optional
        Log-scale the energy / cross-section axis (default ``True``).
    trim : bool, optional
        Tighten the energy axis to the positive data span (default ``True``).
    shade : bool or float, optional
        Shade under each curve.
    show_bands : bool, optional
        Overlay band-averaged cross sections as bars for every reaction that
        has them.
    palette : list[str] or None, optional
        Colour cycle override.
    title, fig, ax, grid, show, save, filename
        Standard rendering controls.

    Returns
    -------
    tuple[matplotlib.figure.Figure, matplotlib.axes.Axes] or None
        ``None`` when no reaction has cross-section data to plot.

    Notes
    -----
    With a single reaction the legend uses the process names; with several it
    prefixes each with the reaction so curves stay distinguishable.
    """
    logger = JaffLogger().get_logger()
    reaction_list = _as_list(reactions)
    procs = _normalize_processes(processes)
    multi_reaction = len(reaction_list) > 1

    series: list[tuple[np.ndarray, np.ndarray, str]] = []
    band_frames = []
    for r in reaction_list:
        xsecs = getattr(r, "xsecs_dict", None)
        if xsecs is None:
            logger.info(f"No cross sections available for: {r}")
            continue
        xsecs = _sampled_xsecs(xsecs)
        energy = xsecs.get("photon_energy")
        if energy is None:
            continue
        available = [p for p in procs if xsecs.get(p) is not None]
        if not available:
            logger.info(f"No data for requested cross-section(s) {procs} in: {r}")
            continue

        x = np.asarray(_units.convert_energy(np.asarray(energy), "eV", energy_unit))
        for p in available:
            y = np.asarray(_units.convert_xsec(np.asarray(xsecs[p]), "cm2", xsec_unit))
            proc_label = _units.PROCESS_LABELS.get(p, p)
            label = f"{r} — {proc_label}" if multi_reaction else proc_label
            series.append((x, y, label))

        if show_bands and hasattr(r, "band_xsecs"):
            band_frames.append(_frames.band_frame(r.band_xsecs, energy_unit, xsec_unit))

    if not series:
        logger.info("plot_xsecs: nothing to plot.")
        return None

    df = _frames.long_frame(series, "energy", "xsec", "process")

    band_df = None
    if band_frames:
        import pandas as pd

        band_df = pd.concat(band_frames, ignore_index=True)
        if band_df.empty:
            band_df = None

    if not filename:
        if not multi_reaction:
            stem = "cross_sections"
            filename = f"{reaction_list[0]}_{stem}.png"
        else:
            filename = "cross_sections.png"

    if title is None:
        title = _latex_or_str(reaction_list[0]) if not multi_reaction else ""

    return Plotter(palette=palette).render_series(
        df,
        x_col="energy",
        y_col="xsec",
        group_col="process",
        xlabel=_units.energy_label(energy_unit),
        ylabel=_units.xsec_label(xsec_unit),
        log_x=energy_log,
        log_y=xsecs_log,
        dynamic_x=True,
        trim=trim,
        shade=shade,
        bands_df=band_df,
        layout=layout,
        title=title,
        fig=fig,
        ax=ax,
        grid=grid,
        show=show,
        save=save,
        filename=filename,
    )
