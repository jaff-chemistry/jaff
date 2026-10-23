# ABOUTME: Figure / axes / series specification dataclasses for jaff.plotting
# ABOUTME: Plain nested dicts are validated into frozen dataclasses by from_dict

"""Declarative figure specification.

A figure is described as data: a ``grid`` of axes, each holding a list of
series.  Users write plain nested dicts; :meth:`FigureSpec.from_dict`
validates them into frozen dataclasses so the renderer only ever sees
well-formed input.  Every validation failure raises :class:`ValueError` whose
message starts with the offending path (e.g. ``axes[1].series[0]``).
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
from typing import Any, NoReturn

import numpy as np

from .themes import resolve_palette, resolve_style

SCALES: tuple[str, ...] = ("linear", "log")
FIGURE_LEGENDS: tuple[str, ...] = ("auto", "axes", "figure", "none")
AXES_LEGENDS: tuple[str, ...] = ("axes", "none")
KINDS: tuple[str, ...] = ("line", "bars")
ROLES: tuple[str, ...] = ("data", "reference")

_SERIES_KEYS = frozenset(
    {"x", "y", "label", "color", "style", "lw", "zorder", "kind", "x_lo", "x_hi", "role"}
)
_AXES_KEYS = frozenset(
    {
        "series",
        "xscale",
        "yscale",
        "xlim",
        "ylim",
        "xlabel",
        "ylabel",
        "title",
        "tag",
        "legend",
    }
)
_FIGURE_KEYS = frozenset(
    {
        "grid",
        "axes",
        "style",
        "palette",
        "figsize",
        "sharex",
        "sharey",
        "legend",
        "legend_title",
        "title",
    }
)


def _fail(path: str, msg: str) -> NoReturn:
    """Raise a ``ValueError`` prefixed with the spec *path*."""
    raise ValueError(f"{path}: {msg}")


def _check_keys(d: Any, allowed: frozenset[str], path: str) -> None:
    """Require *d* to be a dict whose keys are all in *allowed*."""
    if not isinstance(d, dict):
        _fail(path, f"expected a dict, got {type(d).__name__}")
    unknown = sorted(set(d) - allowed)
    if unknown:
        _fail(path, f"unknown key(s) {unknown}; allowed: {sorted(allowed)}")


def _require(d: dict, key: str, path: str) -> Any:
    """Return ``d[key]`` or fail naming the missing key."""
    if key not in d:
        _fail(path, f"missing required key {key!r}")
    return d[key]


def _choice(value: Any, choices: tuple[str, ...], path: str) -> str:
    """Require *value* to be one of *choices*."""
    if value not in choices:
        _fail(path, f"{value!r} not in {list(choices)}")
    return value


def _array(value: Any, path: str) -> np.ndarray:
    """Copy *value* into a read-only 1-D float array."""
    try:
        arr = np.array(value, dtype=float)
    except (TypeError, ValueError) as exc:
        _fail(path, f"not numeric data ({exc})")
    if arr.ndim != 1:
        _fail(path, f"expected 1-D data, got shape {arr.shape}")
    arr.setflags(write=False)
    return arr


def _pair(value: Any, path: str) -> tuple[float, float] | None:
    """Coerce an optional 2-sequence of numbers (limits, figsize)."""
    if value is None:
        return None
    if isinstance(value, str):
        _fail(path, f"expected a pair of numbers, got {value!r}")
    try:
        a, b = value
        return (float(a), float(b))
    except (TypeError, ValueError):
        _fail(path, f"expected a pair of numbers, got {value!r}")


def _grid(value: Any, path: str) -> tuple[int, int]:
    """Coerce ``grid`` to two positive ints."""
    try:
        m, n = value
    except (TypeError, ValueError):
        _fail(path, f"expected (rows, cols), got {value!r}")
    for v in (m, n):
        if isinstance(v, bool) or not isinstance(v, (int, np.integer)) or v < 1:
            _fail(path, f"expected two positive ints, got {value!r}")
    return (int(m), int(n))


def _text(d: dict, key: str) -> str:
    """Return ``d[key]`` as text; missing or ``None`` gives ``""``."""
    value = d.get(key)
    return "" if value is None else str(value)


def _opt_float(value: Any, path: str) -> float | None:
    """Coerce an optional number."""
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        _fail(path, f"expected a number, got {value!r}")


@dataclass(frozen=True, eq=False)
class SeriesSpec:
    """One curve (or set of band bars) on an axes.

    Written as a plain dict in an axes' ``"series"`` list (keys as the
    attributes below; only ``x`` and ``y`` are required) and validated by
    :meth:`from_dict`.

    Attributes
    ----------
    x, y : numpy.ndarray
        Read-only 1-D data of equal length.
    label : str
        Legend / tag text; empty (or ``None`` in the dict) means unlabelled.
    color : str or None
        Explicit colour string; ``None`` resolves from the palette (its
        label map, then the cycle).
    style : str or tuple or None
        Linestyle (``"-"``, ``(0, (6, 3))``), marker (``"o"``, ``"X"``,
        ``"D"``) or both as a fmt string (``"-o"``); ``None`` resolves from
        the style cycle.
    lw : float or None
        Line (or marker-edge) width; ``None`` uses the style default.
    zorder : float or None
        Draw order.
    kind : str
        ``"line"`` (default) or ``"bars"`` (horizontal bars from ``x_lo`` to
        ``x_hi`` at height ``y``).
    x_lo, x_hi : numpy.ndarray or None
        Bar edges (``kind == "bars"`` only, and required there); ``y`` is
        the bar height.
    role : str
        ``"data"`` (default) or ``"reference"`` (thick dark-grey line drawn
        beneath the rest, outside the colour cycle).

    Notes
    -----
    Non-finite points are dropped when rendering, as are non-positive values
    on log axes; a series left with fewer than two points (bars: none) is
    skipped with a warning.

    Examples
    --------
    >>> import numpy as np
    >>> from jaff.plotting import SeriesSpec
    >>> x = np.logspace(1, 4, 50)
    >>> s = SeriesSpec.from_dict({"x": x, "y": 1e-10 * x**0.5, "style": "--D"})
    >>> s.kind, s.role, s.x.flags.writeable
    ('line', 'data', False)
    """

    x: np.ndarray
    y: np.ndarray
    label: str = ""
    color: str | None = None
    style: str | tuple | None = None
    lw: float | None = None
    zorder: float | None = None
    kind: str = "line"
    x_lo: np.ndarray | None = None
    x_hi: np.ndarray | None = None
    role: str = "data"

    @classmethod
    def from_dict(cls, d: dict, path: str = "series") -> SeriesSpec:
        """Validate a series dict.

        Raises
        ------
        ValueError
            On unknown keys, missing/ragged data, or bad ``kind``/edges.
        """
        _check_keys(d, _SERIES_KEYS, path)
        x = _array(_require(d, "x", path), f"{path}.x")
        y = _array(_require(d, "y", path), f"{path}.y")
        if x.size != y.size:
            _fail(path, f"x and y lengths differ ({x.size} vs {y.size})")
        kind = _choice(d.get("kind", "line"), KINDS, f"{path}.kind")
        x_lo = x_hi = None
        if kind == "bars":
            if "x_lo" not in d or "x_hi" not in d:
                _fail(path, "kind 'bars' requires 'x_lo' and 'x_hi'")
            x_lo = _array(d["x_lo"], f"{path}.x_lo")
            x_hi = _array(d["x_hi"], f"{path}.x_hi")
            if x_lo.size != y.size or x_hi.size != y.size:
                _fail(path, "x_lo and x_hi must have the same length as y")
        elif "x_lo" in d or "x_hi" in d:
            _fail(path, "'x_lo'/'x_hi' are only valid for kind 'bars'")
        color, style = d.get("color"), d.get("style")
        if color is not None and not isinstance(color, str):
            _fail(f"{path}.color", f"expected a colour string, got {color!r}")
        if style is not None and not isinstance(style, (str, tuple)):
            _fail(f"{path}.style", f"expected a string or tuple, got {style!r}")
        return cls(
            x=x,
            y=y,
            label=_text(d, "label"),
            color=color,
            style=style,
            lw=_opt_float(d.get("lw"), f"{path}.lw"),
            zorder=_opt_float(d.get("zorder"), f"{path}.zorder"),
            kind=kind,
            x_lo=x_lo,
            x_hi=x_hi,
            role=_choice(d.get("role", "data"), ROLES, f"{path}.role"),
        )


@dataclass(frozen=True, eq=False)
class AxesSpec:
    """One panel: its series plus scale, limits, labels, tag and legend mode.

    ``title`` names the panel (drawn inside, upper right, or above when too
    wide); ``tag`` is an extra inside label under it (``None`` / ``""``: none).
    ``legend=None`` inherits the figure policy.

    Written as a plain dict in the figure's ``"axes"`` list (keys as the
    attributes below; only ``series`` is required) and validated by
    :meth:`from_dict`.

    Attributes
    ----------
    series : tuple[SeriesSpec, ...]
        The panel's series (a list of series dicts in the spec; may be
        empty).
    xscale, yscale : str
        ``"linear"`` (default) or ``"log"``.
    xlim, ylim : tuple[float, float] or None
        Axis limits; ``None`` autoscales.
    xlabel, ylabel : str
        Axis labels (``None`` in the dict counts as empty).
    title : str
        Panel name, drawn inside the panel (upper right) with headroom above
        the data, or above the panel when wider than ~85 % of it.
    tag : str or None
        Extra inside label under the title.
    legend : str or None
        ``None`` follows the figure policy; ``"axes"`` gives the panel its
        own legend group, ``"none"`` none.
    """

    series: tuple[SeriesSpec, ...] = ()
    xscale: str = "linear"
    yscale: str = "linear"
    xlim: tuple[float, float] | None = None
    ylim: tuple[float, float] | None = None
    xlabel: str = ""
    ylabel: str = ""
    title: str = ""
    tag: str | None = None
    legend: str | None = None

    @classmethod
    def from_dict(cls, d: dict, path: str = "axes") -> AxesSpec:
        """Validate an axes dict (series validated recursively).

        Raises
        ------
        ValueError
            On unknown keys or invalid values; the message names the path.
        """
        _check_keys(d, _AXES_KEYS, path)
        raw = _require(d, "series", path)
        if not isinstance(raw, (list, tuple)):
            _fail(f"{path}.series", f"expected a list, got {type(raw).__name__}")
        legend = d.get("legend")
        if legend is not None:
            _choice(legend, AXES_LEGENDS, f"{path}.legend")
        tag = d.get("tag")
        return cls(
            series=tuple(
                SeriesSpec.from_dict(s, f"{path}.series[{i}]") for i, s in enumerate(raw)
            ),
            xscale=_choice(d.get("xscale", "linear"), SCALES, f"{path}.xscale"),
            yscale=_choice(d.get("yscale", "linear"), SCALES, f"{path}.yscale"),
            xlim=_pair(d.get("xlim"), f"{path}.xlim"),
            ylim=_pair(d.get("ylim"), f"{path}.ylim"),
            xlabel=_text(d, "xlabel"),
            ylabel=_text(d, "ylabel"),
            title=_text(d, "title"),
            tag=None if tag is None else str(tag),
            legend=legend,
        )


@dataclass(frozen=True, eq=False)
class FigureSpec:
    """A whole figure: grid shape, row-major axes, theme and legend policy.

    ``legend_title`` (optional) heads the shared legend, or the stack of
    per-panel legend groups with ``legend="axes"``.

    ``style`` / ``palette`` are validated in :meth:`from_dict` but stored
    unresolved (a deep copy of the name, dict or list); the renderer resolves
    them with :func:`jaff.plotting.themes.resolve_style` /
    :func:`~jaff.plotting.themes.resolve_palette`.

    Users normally write the figure as a plain nested dict (keys as the
    attributes below; only ``axes`` is required) and pass it to
    :func:`~jaff.plotting.plot`, which validates it with :meth:`from_dict`.

    Attributes
    ----------
    axes : tuple[AxesSpec, ...]
        The panels (a list of axes dicts in the spec), filling the grid row
        by row; there may be fewer than grid cells, and the leftover cells
        are hidden.
    grid : tuple[int, int]
        ``(rows, cols)``, default ``(1, 1)``.
    style : str or dict or None
        Style name or partial style dict (see
        :func:`~jaff.plotting.register_style`); ``None`` is ``"paper"``.
    palette : str or list or dict or None
        Palette name, colour list, or ``{label: colour}`` map; ``None`` is
        ``"jaff"``.
    figsize : tuple[float, float] or None
        Figure size in inches.  ``None``: grid × the style's ``panel_size``,
        grown for right-margin legends, with sizes scaled for wide grids;
        given, the figure keeps it and the panels take what is left.
    sharex, sharey : bool
        Share x / y ranges across panels; all axes must then use the same
        ``xscale`` / ``yscale``.
    legend : str
        ``"auto"`` (default) or ``"figure"``: one shared legend;
        ``"axes"``: one legend group per panel; ``"none"``.
    title : str or None
        Figure title.
    legend_title : str or None
        Bold heading of the shared legend, or of the stack of legend groups.

    Examples
    --------
    >>> from jaff.plotting import FigureSpec
    >>> spec = FigureSpec.from_dict(
    ...     {"grid": (1, 2), "axes": [{"series": [{"x": [1, 2], "y": [3, 4]}]}]}
    ... )
    >>> spec.grid, len(spec.axes), spec.legend
    ((1, 2), 1, 'auto')

    Invalid input raises a ``ValueError`` naming the offending path:

    >>> FigureSpec.from_dict({"axes": [{"series": [{"x": [1, 2], "y": [3]}]}]})
    Traceback (most recent call last):
    ...
    ValueError: axes[0].series[0]: x and y lengths differ (2 vs 1)
    """

    axes: tuple[AxesSpec, ...]
    grid: tuple[int, int] = (1, 1)
    style: str | dict | None = None
    palette: str | list | dict | None = None
    figsize: tuple[float, float] | None = None
    sharex: bool = False
    sharey: bool = False
    legend: str = "auto"
    title: str | None = None
    legend_title: str | None = None

    @classmethod
    def from_dict(cls, d: dict) -> FigureSpec:
        """Validate a figure dict.

        Raises
        ------
        ValueError
            On any invalid key or value; the message names the path.
        """
        _check_keys(d, _FIGURE_KEYS, "figure")
        raw = _require(d, "axes", "figure")
        if not isinstance(raw, (list, tuple)):
            _fail("figure.axes", f"expected a list, got {type(raw).__name__}")
        grid = _grid(d.get("grid", (1, 1)), "figure.grid")
        if len(raw) > grid[0] * grid[1]:
            _fail("figure.axes", f"{len(raw)} axes do not fit grid {grid}")
        title = d.get("title")
        legend_title = d.get("legend_title")
        if legend_title is not None and not isinstance(legend_title, str):
            _fail("figure.legend_title", f"expected a string, got {legend_title!r}")
        for key in ("sharex", "sharey"):
            if not isinstance(d.get(key, False), bool):
                _fail(f"figure.{key}", f"expected a bool, got {d[key]!r}")
        try:
            resolve_style(d.get("style"))
        except ValueError as exc:
            raise ValueError(f"figure.style: {exc}") from exc
        try:
            resolve_palette(d.get("palette"))
        except ValueError as exc:
            raise ValueError(f"figure.palette: {exc}") from exc
        axes = tuple(AxesSpec.from_dict(a, f"axes[{i}]") for i, a in enumerate(raw))
        for key, attr in (("sharex", "xscale"), ("sharey", "yscale")):
            scales = sorted({getattr(a, attr) for a in axes})
            if d.get(key, False) and len(scales) > 1:
                _fail(
                    f"figure.{key}", f"all axes must use the same {attr} (got {scales})"
                )
        return cls(
            axes=axes,
            grid=grid,
            style=copy.deepcopy(d.get("style")),
            palette=copy.deepcopy(d.get("palette")),
            figsize=_pair(d.get("figsize"), "figure.figsize"),
            sharex=d.get("sharex", False),
            sharey=d.get("sharey", False),
            legend=_choice(d.get("legend", "auto"), FIGURE_LEGENDS, "figure.legend"),
            title=None if title is None else str(title),
            legend_title=legend_title,
        )
