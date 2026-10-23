# ABOUTME: Public plotting entry points: plot(spec), plot_rates, plot_xsecs
# ABOUTME: Builders group inputs into panels, call adapters, render, save/show

"""Plotting entry points.

:func:`plot` renders a figure spec dict.  :func:`plot_rates` and
:func:`plot_xsecs` build that spec from reactions: each top-level item is one
panel, a nested list is overlaid on one axes, and without ``grid`` every
curve shares a single axes -- unless there are more curves than the style's
``legend_max``, in which case they are split in order over an automatic grid.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from sympy import Basic

from ..io import JaffLogger
from . import _units
from .adapters import rate_series, xsec_series
from .render import effective_style, min_points, plottable_mask, render
from .spec import FigureSpec
from .themes import resolve_style

_VECTOR = frozenset({".pdf", ".svg", ".svgz", ".eps", ".ps", ".pgf"})
_NON_INTERACTIVE = frozenset({"agg", "pdf", "ps", "svg", "pgf", "cairo", "template"})

Result = tuple[Figure, np.ndarray] | None


def plot(
    spec: dict | FigureSpec,
    *,
    save: str | Path | None = None,
    show: bool = True,
    close: bool = True,
) -> Result:
    """Render a figure spec; optionally save and show it.

    The spec is validated before anything is drawn, and a ``ValueError``
    names the offending path, e.g.
    ``axes[1].series[0]: x and y lengths differ (200 vs 199)``.

    Parameters
    ----------
    spec : dict or FigureSpec
        Figure specification: a plain nested dict (schema in the Notes; keys
        as the attributes of :class:`~jaff.plotting.FigureSpec`,
        :class:`~jaff.plotting.AxesSpec` and :class:`~jaff.plotting.SeriesSpec`)
        or an already validated :class:`~jaff.plotting.FigureSpec`.
    save : str or pathlib.Path or None
        Output path; the format follows the suffix.  Raster formats use the
        style's ``save["dpi"]``.
    show : bool
        Call :func:`matplotlib.pyplot.show` after rendering (skipped on a
        non-interactive backend such as Agg).
    close : bool
        Remove the figure from pyplot (:func:`matplotlib.pyplot.close`) after
        saving / showing, so repeated calls do not pile up open figures.  The
        returned figure can still be edited and saved with ``fig.savefig``.
        Pass ``False`` to keep it open in pyplot (e.g. for inline notebook
        display, or when you need the canvas itself).  Ignored after an
        interactive, non-blocking show (``plt.isinteractive()``, e.g. IPython
        with a Qt/Tk or ipympl backend), which would otherwise close the
        window at once.

    Returns
    -------
    tuple[Figure, numpy.ndarray] or None
        ``(fig, axes)`` with ``axes`` a 2-D array shaped like the grid, or
        ``None`` when the spec holds no series at all.

    Raises
    ------
    ValueError
        Invalid spec, style or palette.

    Notes
    -----
    Schema (values shown are the defaults unless noted):

    ```python
    {
        "grid": (1, 1),  # (rows, cols); axes fill row-major
        "style": None,  # name or partial dict; None: "paper"
        "palette": None,  # name, colour list or {label: colour}; None: "jaff"
        "figsize": None,  # None: grid x style panel size, grown for legends
        "sharex": False,
        "sharey": False,
        "legend": "auto",  # auto | axes | figure | none
        "legend_title": None,  # bold heading of the shared legend / the groups
        "title": None,  # figure title
        "axes": [  # required
            {
                "series": [...],  # required
                "xscale": "linear",  # or "log"
                "yscale": "linear",  # or "log"
                "xlim": None,
                "ylim": None,
                "xlabel": "",
                "ylabel": "",
                "title": "",  # panel name, drawn inside (upper right)
                "tag": None,  # extra inside label, under the title
                "legend": None,  # None: figure policy; "axes" | "none"
            },
        ],
    }
    ```

    Series:

    ```python
    {
        "x": array,  # required; 1-D
        "y": array,  # required; same length as x
        "label": "",
        "color": None,  # None: palette (label map, then cycle)
        "style": None,  # linestyle, marker ("o", "X", "D") or both ("-o")
        "role": "data",  # "reference": thick dark-grey line beneath the rest
        "lw": None,
        "zorder": None,
        "kind": "line",  # "bars": horizontal bars, needs "x_lo" and "x_hi"
    }
    ```

    Non-finite points are dropped, as are non-positive values on log axes; a
    series left with fewer than two points (bars: none) is skipped with a
    warning.

    Layout:

    - Tight grid: the gap between panels is ``style.gap`` inches in both
      directions; margins are sized to the outer tick and axis labels.  Rows
      and columns of the grid with no used panel take no space (8 panels in
      a 4×4 grid are laid out as 2×4); the returned ``axes`` keep the
      requested shape, with those cells hidden.  Hidden cells inside used
      rows stay (legend host).  Without ``figsize`` each panel box is
      ``style.panel_size`` and the figure is made tall enough for
      right-margin legends; with ``figsize`` the figure is kept as given and
      the panels take what is left.
    - X tick labels by place in the column (hidden cells don't count): the
      topmost panel puts its x tick labels (and x label) on the top spine;
      the bottommost panel (also one above a hidden cell, or the only one) at
      the bottom; a middle panel shows none when its x range equals the
      bottommost panel's, else draws them inside above the bottom spine, with
      its y range widened downward so no curve runs under them (with an
      explicit ``ylim`` they stay outside instead).  Ticks stay on all four
      sides.
    - Axis labels shown by every panel are drawn once: a common y label over
      two or more rows becomes ``fig.supylabel``, a common x label over two
      or more columns becomes ``fig.supxlabel``.  Otherwise a middle panel
      drops an x label its column's bottom panel repeats, and a y label is
      dropped where the panel to the left repeats it.
    - With exponent ticks (``"paper"``), log axes get at most ~4 labelled
      decades (every 1, 2, 3, 4, 5, 10, ... decades) and unlabelled minor
      ticks on decades (at most 4 between labels); an axis spanning less
      than 1.5 decades gets ticks at fractional exponents (``1.05``,
      ``1.1``, ...) and no minor ticks.  These are fixed positions chosen
      from the final limits: changing the limits of the returned axes does
      not re-tick them (nor re-space hollow markers).
    - With inside tick labels (``"paper"``), the x range is widened on the
      left by the widest y label (plus the tick length); a panel with an
      explicit ``xlim`` keeps its y tick labels outside instead.
    - Headroom for inside titles is applied once per shared y range
      (``sharey``), using the largest data span of the panels sharing it.
    - Legends: see :mod:`jaff.plotting`.  A legend with more entries than the
      style's ``legend_max`` is not drawn, with a warning naming the axes (or
      ``figure``); a shared legend over the cap that keys several panels
      falls back to per-panel groups.

    Validation:

    - With ``sharex`` / ``sharey``, all axes must use the same ``xscale`` /
      ``yscale``.
    - ``color`` must be a string; ``style`` must be a string or a tuple.
    - Text fields (``label``, ``title``, ``xlabel``, ...) given as ``None``
      count as empty.
    - Style and palette names are checked when the spec is validated; errors
      are prefixed ``figure.style`` / ``figure.palette``.

    Examples
    --------
    Three panels on a 2×2 grid; the shared legend fills the free
    bottom-right cell:

    >>> import numpy as np
    >>> from jaff.plotting import plot
    >>> nh = np.logspace(17, 22, 300)
    >>> def panel(name, ref, sim):
    ...     return {
    ...         "title": name,
    ...         "xscale": "log",
    ...         "yscale": "log",
    ...         "xlabel": r"$N_{\\rm H}$ (cm$^{-2}$)",
    ...         "ylabel": r"$n_i/n_{\\rm H}$",
    ...         "series": [
    ...             {"x": nh, "y": ref, "label": "Gong+2017", "role": "reference"},
    ...             {"x": nh, "y": sim, "label": "simulation"},
    ...             {"x": nh[::20], "y": sim[::20], "label": "RT off", "style": "o"},
    ...         ],
    ...     }
    >>> drop = 1 / (1 + nh / 1e20)
    >>> fig, axes = plot(
    ...     {
    ...         "grid": (2, 2),
    ...         "sharex": True,
    ...         "axes": [
    ...             panel("CO", 1e-4 * drop, 1.1e-4 * drop),
    ...             panel("C", 1e-5 * drop, 0.9e-5 * drop),
    ...             panel("O", 3e-4 * drop, 3e-4 * drop),
    ...         ],
    ...     },
    ...     save="abundances.pdf",
    ... )
    >>> axes.shape
    (2, 2)
    """
    fs = spec if isinstance(spec, FigureSpec) else FigureSpec.from_dict(spec)
    if not any(a.series for a in fs.axes):
        JaffLogger().get_logger().info("plot: nothing to plot.")
        return None
    fig, axes = render(fs)
    if save is not None:
        _save(fig, Path(save), fs)
    shown = show and matplotlib.get_backend().lower() not in _NON_INTERACTIVE
    if shown:
        plt.show()
    # A non-blocking show (interactive mode) leaves the window open: keep it.
    if close and not (shown and plt.isinteractive()):
        plt.close(fig)
    return fig, axes


def _save(fig: Figure, path: Path, spec: FigureSpec) -> None:
    """Save *fig* with the spec's style (scoped rcParams, tight bbox)."""
    style = effective_style(spec)  # scaled rc, as when rendering
    kw: dict[str, Any] = {
        "facecolor": style.surface,
        "bbox_inches": "tight",
        "pad_inches": style.save["pad_inches"],
    }
    if path.suffix.lower() not in _VECTOR:
        kw["dpi"] = style.save["dpi"]
    with plt.rc_context(style.rc):
        fig.savefig(path, **kw)


def _is_group(item: Any, pair_leaves: bool) -> bool:
    """Whether *item* is a container of leaves rather than a leaf itself.

    Any iterable except ``str``/``bytes``/``dict``/SymPy objects is a group.
    A tuple is a group too, unless *pair_leaves* is set and it is an
    ``(x, y)`` pair of 1-D numeric array-likes (a tuple of two such pairs is
    a group of pairs).
    """
    if isinstance(item, (str, bytes, dict, Basic)) or not hasattr(item, "__iter__"):
        return False
    if isinstance(item, tuple) and pair_leaves and len(item) == 2:
        try:
            dims = [np.asarray(part, dtype=float).ndim for part in item]
        except (TypeError, ValueError):
            return True
        return any(d > 1 for d in dims)
    return True


def _as_panels(
    items: Any, grid: tuple[int, int] | None, pair_leaves: bool = False
) -> list[list[Any]]:
    """Group inputs into panels: top-level item = panel, nested group = overlay.

    Without *grid* everything is flattened into a single panel.
    """
    top = list(items) if _is_group(items, pair_leaves) else [items]
    panels = [list(entry) if _is_group(entry, pair_leaves) else [entry] for entry in top]
    if grid is None:
        return [[leaf for panel in panels for leaf in panel]]
    return panels


def _per_panel(value: str | list[str], n: int, name: str) -> list[str]:
    """Broadcast a scale to *n* panels, or check a per-panel list's length."""
    if isinstance(value, str):
        return [value] * n
    if len(value) != n:
        raise ValueError(f"{name}: got {len(value)} values for {n} panels")
    return list(value)


def _auto_grid(n: int) -> tuple[int, int]:
    """Grid for *n* automatic panels: one row up to 3, else near-square."""
    if n <= 3:
        return (1, n)
    cols = math.ceil(math.sqrt(n))
    return (math.ceil(n / cols), cols)


def _pack(weights: list[int], cap: int) -> list[list[int]]:
    """Split indices, in order, into runs whose weights sum to at most *cap*.

    An item heavier than *cap* gets a run of its own.
    """
    runs: list[list[int]] = []
    total = cap + 1
    for i, w in enumerate(weights):
        if total + w > cap:
            runs.append([])
            total = 0
        runs[-1].append(i)
        total += w
    return runs


def _plottable(series: dict, xscale: str, yscale: str) -> bool:
    """Whether a series dict keeps enough points on the given scales to draw."""
    kind = series.get("kind", "line")
    try:
        keep = plottable_mask(
            kind,
            series["x"],
            series["y"],
            xscale,
            yscale,
            series.get("x_lo"),
            series.get("x_hi"),
        )
    except (TypeError, ValueError):
        return True  # malformed data: let spec validation report it
    return int(keep.sum()) >= min_points(kind)


def _report_unplottable(builder: str, n: int) -> None:
    """One summary warning for series dropped before panel assembly."""
    if n:
        JaffLogger().get_logger().warning(
            f"{builder}: {n} series have no plottable points on the chosen scales; "
            "skipped."
        )


def _figure_legend_fallback(builder: str, legend: str) -> str:
    """Per-panel legend groups for automatic panels.

    One shared legend would hold every curve of every panel (more than
    ``legend_max``), so ``"auto"`` becomes ``"axes"``; an explicit
    ``"figure"`` does too, with a warning.
    """
    if legend == "figure":
        JaffLogger().get_logger().warning(
            f"{builder}: legend='figure' cannot cover automatically split panels; "
            "using per-panel legends."
        )
    return "axes" if legend in ("auto", "figure") else legend


def _chunk_scale(value: str | list[str], name: str, n: int) -> list[str]:
    """Broadcast a scalar scale over *n* automatic panels (lists need a grid)."""
    if isinstance(value, str):
        return [value] * n
    raise ValueError(
        f"{name}: a per-panel list needs an explicit grid; the curves were split "
        f"over {n} automatic panels (legend_max), so pass a single scale or a grid"
    )


def _figure(axes: list[dict], grid: tuple[int, int] | None, **opts: Any) -> dict:
    """Assemble a figure spec dict, omitting ``None`` options."""
    if grid is not None and len(axes) > grid[0] * grid[1]:
        raise ValueError(
            f"grid {grid} has {grid[0] * grid[1]} cells but {len(axes)} panels have data"
        )
    spec = {"grid": grid or (1, 1), "axes": axes}
    spec.update({k: v for k, v in opts.items() if v is not None})
    return spec


def plot_rates(
    rates: Any,
    *,
    grid: tuple[int, int] | None = None,
    tmin: float | None = None,
    tmax: float | None = None,
    var: str = "tgas",
    npoints: int = 200,
    xscale: str | list[str] = "log",
    yscale: str | list[str] = "log",
    xlabel: str = r"$T$ (K)",
    ylabel: str = r"$\kappa$",
    labels: list[str] | None = None,
    style: str | dict | None = None,
    palette: str | list | dict | None = None,
    sharex: bool = False,
    sharey: bool = False,
    legend: str = "auto",
    title: str | None = None,
    save: str | Path | None = None,
    show: bool = True,
    close: bool = True,
) -> Result:
    """Plot rate coefficients against temperature.

    Parameters
    ----------
    rates : Reaction or sympy.Basic or tuple or list
        A reaction (anything with ``rate_coefficient``, optionally ``tmin`` /
        ``tmax`` / ``get_latex``), a SymPy expression (needs *tmin* and
        *tmax*), or an ``(x, y)`` tuple used verbatim; or a list of them.
        With *grid*, each top-level item is a panel and a nested list is
        overlaid on one panel.  Without *grid* everything is overlaid,
        unless there are more curves than the style's ``legend_max``: then
        they are split in order into panels of at most ``legend_max`` curves
        on an automatic grid (up to 3 panels in one row, else
        ``cols = ceil(sqrt(n))``, ``rows = ceil(n / cols)``).
    grid : tuple[int, int] or None
        ``(rows, cols)``; panels fill row-major.  An explicit grid is never
        split automatically.
    tmin : float or None
        Lower temperature bound (K).  A reaction falls back to its own
        ``tmin`` (when finite and positive), then to 2.73 K.
    tmax : float or None
        Upper temperature bound (K); falls back to the reaction's ``tmax``,
        then to 1e6 K.
    var : str
        Temperature symbol substituted in reaction rates.
    npoints : int
        Log-spaced samples per curve.
    xscale : str or list[str]
        ``"log"`` / ``"linear"``, for all panels or one entry per panel (a
        list needs *grid* once the curves are split automatically).
    yscale : str or list[str]
        As *xscale*, for the rate axis.
    xlabel : str
        Temperature axis label for every panel.
    ylabel : str
        Rate axis label for every panel; the default ``$\\kappa$`` (no
        units) reads as "log₁₀ κ" with the ``"paper"`` exponent ticks.  A
        label shared by all panels is drawn once for the figure (over two or
        more columns / rows); otherwise it is kept only on the outer panels.
    labels : list[str] or None
        Legend labels, one per leaf item in input order.  Default: the
        reactions' LaTeX equations.
    style : str or dict or None
        Style name or partial style dict (see :class:`FigureSpec`).
    palette : str or list or dict or None
        Palette name, colour list or ``{label: colour}`` map.
    sharex : bool
        Share the x range across panels.
    sharey : bool
        Share the y range across panels.
    legend : str
        Figure legend policy: ``"auto"``, ``"figure"``, ``"axes"`` or
        ``"none"``.
    title : str or None
        Figure title.
    save : str or pathlib.Path or None
        As in :func:`plot`.
    show : bool
        As in :func:`plot`.
    close : bool
        As in :func:`plot`.

    Returns
    -------
    tuple[Figure, numpy.ndarray] or None
        ``None`` when no rate could be evaluated.

    Raises
    ------
    ValueError
        If the temperature range does not satisfy ``0 < tmin < tmax``, a
        bare SymPy expression is given without *tmin* and *tmax*, *labels*
        or a per-panel scale list has the wrong length, or more panels have
        data than *grid* has cells.
    TypeError
        An unsupported item type.

    Notes
    -----
    Series with fewer than two plottable points on the chosen scales (e.g. a
    zero rate on a log axis) are dropped before panels are assembled, with
    one summary warning, so the automatic split never produces blank panels.
    With the split, the legend switches to per-panel groups in the right
    margin (``legend="axes"``); an explicit ``legend="figure"`` does too,
    with a warning.  Otherwise one shared legend keys every panel.

    A panel holding a single curve shows its equation as the panel title
    (inside the panel, or above it when too long).  Rates that cannot be
    evaluated against temperature alone (photo, cosmic-ray,
    interpolated-function rates) are skipped with one summary warning (the
    per-reaction reason is logged at debug level); a panel left empty is
    dropped.

    Examples
    --------
    >>> from jaff import Network
    >>> from jaff.plotting import plot_rates
    >>> net = Network("GOW")
    >>> r = list(net.reactions.with_type("two_body"))
    >>> fig, axes = plot_rates(r[:5], tmin=10, tmax=1e4)  # overlay
    >>> fig, axes = plot_rates(r)  # 15 evaluable curves: split 12 + 3
    >>> axes.shape
    (1, 2)
    >>> fig, axes = plot_rates(r, style={"legend_max": 20})  # one panel
    >>> fig, axes = plot_rates([r[0], r[1:4]], grid=(1, 2), yscale=["log", "linear"])
    """
    panels = _as_panels(rates, grid, pair_leaves=True)
    n_leaves = sum(len(p) for p in panels)
    if labels is not None and len(labels) != n_leaves:
        raise ValueError(f"labels: got {len(labels)} labels for {n_leaves} items")
    label_iter = iter(labels or [None] * n_leaves)
    xs = _per_panel(xscale, len(panels), "xscale")
    ys = _per_panel(yscale, len(panels), "yscale")
    logger = JaffLogger().get_logger()

    evaluated: list[tuple[int, list[dict]]] = []
    unplottable = unevaluable = 0
    for i, panel in enumerate(panels):
        if not panel:
            continue
        series = []
        for item in panel:
            s = rate_series(
                item,
                tmin=tmin,
                tmax=tmax,
                var=var,
                npoints=npoints,
                label=next(label_iter),
                quiet=True,
            )
            if s is None:
                unevaluable += 1
            else:
                series.append(s)
        kept = [s for s in series if _plottable(s, xs[i], ys[i])]
        unplottable += len(series) - len(kept)
        if not kept:
            # Causes are summarised once below; per-panel detail at debug level.
            logger.debug(f"plot_rates: panel {i} has nothing to plot; dropped.")
            continue
        evaluated.append((i, kept))
    if unevaluable:
        logger.warning(
            f"plot_rates: {unevaluable} rates cannot be evaluated vs temperature; "
            "skipped."
        )
    _report_unplottable("plot_rates", unplottable)
    if not evaluated:
        logger.info("plot_rates: nothing to plot.")
        return None
    cap = resolve_style(style).legend_max
    if grid is None and len(evaluated[0][1]) > cap:
        flat = evaluated[0][1]
        runs = _pack([1] * len(flat), cap)
        evaluated = [(k, [flat[j] for j in run]) for k, run in enumerate(runs)]
        xs = _chunk_scale(xscale, "xscale", len(runs))
        ys = _chunk_scale(yscale, "yscale", len(runs))
        grid = _auto_grid(len(runs))
        legend = _figure_legend_fallback("plot_rates", legend)

    axes = []
    for i, series in evaluated:
        panel_title = ""
        if len(series) == 1:
            # Lone curve: its (often long) label becomes the panel title.
            panel_title, series[0]["label"] = series[0]["label"], ""
        axes.append(
            {
                "series": series,
                "xscale": xs[i],
                "yscale": ys[i],
                "xlabel": xlabel,
                "ylabel": ylabel,
                "title": panel_title,
            }
        )
    spec = _figure(
        axes,
        grid,
        style=style,
        palette=palette,
        sharex=sharex,
        sharey=sharey,
        legend=legend,
        title=title,
    )
    return plot(spec, save=save, show=show, close=close)


def _prefixed(
    panel: list[Any], processes: Any, energy_unit: str, xsec_unit: str, bands: bool
) -> list[tuple[Any, list[dict]]]:
    """``(reaction, series)`` for each reaction with data, labels prefixed."""
    per_reaction = [
        (
            r,
            xsec_series(
                r,
                processes=processes,
                energy_unit=energy_unit,
                xsec_unit=xsec_unit,
                prefix=True,
                bands=bands,
            ),
        )
        for r in panel
    ]
    return [(r, ss) for r, ss in per_reaction if ss]


def _xsec_panel(
    contributing: list[tuple[Any, list[dict]]],
    processes: Any,
    energy_unit: str,
    xsec_unit: str,
    bands: bool,
) -> tuple[list[dict], str]:
    """Series for one panel plus its title (the equation of a lone reaction)."""
    if len(contributing) != 1:
        return [s for _, ss in contributing for s in ss], ""
    r, _ = contributing[0]
    series = xsec_series(
        r,
        processes=processes,
        energy_unit=energy_unit,
        xsec_unit=xsec_unit,
        prefix=False,
        bands=bands,
    )
    return series, r.get_latex() if hasattr(r, "get_latex") else str(r)


def plot_xsecs(
    reactions: Any,
    *,
    grid: tuple[int, int] | None = None,
    processes: str | list[str] | None = "all",
    energy_unit: str = "eV",
    xsec_unit: str = "Mb",
    xscale: str | list[str] = "log",
    yscale: str | list[str] = "log",
    bands: bool = False,
    style: str | dict | None = None,
    palette: str | list | dict | None = None,
    sharex: bool = False,
    sharey: bool = False,
    legend: str = "auto",
    title: str | None = None,
    save: str | Path | None = None,
    show: bool = True,
    close: bool = True,
) -> Result:
    """Plot photo cross sections against photon energy (or wavelength).

    Parameters
    ----------
    reactions : Reaction or list
        A reaction (anything exposing ``photo.xsec``, optionally
        ``band_xsecs`` / ``get_latex``) or a list (nesting / *grid* as in
        :func:`plot_rates`, including the automatic split past
        ``legend_max`` curves).  Automatic splitting keeps all curves of one
        reaction on one panel.
    grid : tuple[int, int] or None
        As in :func:`plot_rates`.
    processes : str or list[str] or None
        ``"all"`` (or ``None``), ``"photo_absorption"``, ``"photodecay"`` or
        a list.
    energy_unit : str
        Photon axis unit: ``"eV"``, ``"erg"``, or the wavelengths ``"nm"`` /
        ``"um"``.
    xsec_unit : str
        Cross-section unit: ``"Mb"``, ``"barn"`` or ``"cm^2"`` (``"cm2"``).
    xscale : str or list[str]
        As in :func:`plot_rates`.
    yscale : str or list[str]
        As in :func:`plot_rates`.
    bands : bool
        Overlay band-average cross sections (from ``Reaction.band_xsecs``)
        as horizontal bars; reactions without bands add nothing.  An open
        band stays open in any unit: in wavelength units the open ionising
        band extends to the short-wavelength edge of the axes.
    style : str or dict or None
        As in :func:`plot_rates`.
    palette : str or list or dict or None
        As in :func:`plot_rates`.
    sharex : bool
        As in :func:`plot_rates`.
    sharey : bool
        As in :func:`plot_rates`.
    legend : str
        As in :func:`plot_rates`.
    title : str or None
        As in :func:`plot_rates`.
    save : str or pathlib.Path or None
        As in :func:`plot`.
    show : bool
        As in :func:`plot`.
    close : bool
        As in :func:`plot`.

    Returns
    -------
    tuple[Figure, numpy.ndarray] or None
        ``None`` when no reaction has cross-section data.

    Raises
    ------
    KeyError
        Unknown process name.
    ValueError
        A per-panel scale list of the wrong length, more panels with data
        than *grid* has cells, or an unknown *energy_unit* / *xsec_unit*.

    Notes
    -----
    A panel holding one reaction is titled with its equation (inside the
    panel, or above it when too long) and its curves are labelled by
    process.  Several reactions on one panel are labelled
    ``equation — process``.  Unless the curves are split automatically, one
    shared legend keys all panels, with one colour per process across the
    figure; with the automatic split the legend switches to per-panel groups
    as in :func:`plot_rates`.

    Examples
    --------
    >>> from jaff import Network
    >>> from jaff.plotting import plot_xsecs
    >>> net = Network(
    ...     "h_photoionization", radiation_props={"bands": [1, 13.6, 100, "inf"]}
    ... )
    >>> photo = [p for p in net.reactions.photo_reactions() if p.photo.xsec]
    >>> fig, axes = plot_xsecs(photo[0], bands=True)
    >>> [t.get_text() for t in fig.legends[0].get_texts()]
    ['Photodecay', 'Band average']
    >>> fig, axes = plot_xsecs(
    ...     photo, processes="photodecay", energy_unit="nm", xscale="linear"
    ... )
    """
    panels = _as_panels(reactions, grid)
    xs = _per_panel(xscale, len(panels), "xscale")
    ys = _per_panel(yscale, len(panels), "yscale")
    logger = JaffLogger().get_logger()
    data: list[tuple[int, list[tuple[Any, list[dict]]]]] = []
    unplottable = 0
    for i, panel in enumerate(panels):
        if not panel:
            continue
        contributing = []
        for r, ss in _prefixed(panel, processes, energy_unit, xsec_unit, bands):
            kept = [s for s in ss if _plottable(s, xs[i], ys[i])]
            unplottable += len(ss) - len(kept)
            if kept:
                contributing.append((r, kept))
        data.append((i, contributing))
    _report_unplottable("plot_xsecs", unplottable)

    cap = resolve_style(style).legend_max
    if grid is None and data and sum(len(ss) for _, ss in data[0][1]) > cap:
        # Keep each reaction's curves together on one panel.
        flat = data[0][1]
        runs = _pack([len(ss) for _, ss in flat], cap)
        data = [(k, [flat[j] for j in run]) for k, run in enumerate(runs)]
        xs = _chunk_scale(xscale, "xscale", len(runs))
        ys = _chunk_scale(yscale, "yscale", len(runs))
        grid = _auto_grid(len(runs))
        legend = _figure_legend_fallback("plot_xsecs", legend)

    axes = []
    empty = 0
    for i, contributing in data:
        series, panel_title = _xsec_panel(
            contributing, processes, energy_unit, xsec_unit, bands
        )
        series = [s for s in series if _plottable(s, xs[i], ys[i])]
        if not series:
            logger.debug(f"plot_xsecs: panel {i} has nothing to plot; dropped.")
            empty += 1
            continue
        ax = {
            "series": series,
            "xscale": xs[i],
            "yscale": ys[i],
            "xlabel": _units.energy_label(energy_unit),
            "ylabel": _units.xsec_label(xsec_unit),
            "title": panel_title,
        }
        axes.append(ax)
    if empty:
        logger.warning(f"plot_xsecs: {empty} panels have no cross-section data; dropped.")
    if not axes:
        logger.info("plot_xsecs: nothing to plot.")
        return None
    spec = _figure(
        axes,
        grid,
        style=style,
        palette=palette,
        sharex=sharex,
        sharey=sharey,
        legend=legend,
        title=title,
    )
    return plot(spec, save=save, show=show, close=close)
