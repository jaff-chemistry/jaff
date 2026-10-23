# ABOUTME: Plot themes: Style (look + line/marker cycle) and Palette (colours)
# ABOUTME: Built-in registries, name/dict resolution, user registration

"""Plot themes.

A theme has two independent halves:

- :class:`Style` -- everything but colour: rcParams, panel size, axes look,
  tick geometry, the line-style/marker ``cycle``, marker placement, legend
  and tag formatting, save options.
- :class:`Palette` -- colours: a ``cycle`` plus an optional label -> colour
  ``mapping``.

Both are chosen by built-in name or a plain user dict.  Series without an
explicit colour/style are assigned colour first, then style: series ``i``
gets ``palette.cycle[i % nc]`` and ``style.cycle[i // nc]``.  A series with
``role="reference"`` is drawn thick, dark grey and underneath, outside both
cycles.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, fields, replace
from typing import Any

from matplotlib.colors import is_color_like
from matplotlib.lines import Line2D

_LINESTYLES = frozenset({"-", "--", "-.", ":", "solid", "dashed", "dashdot", "dotted"})
_NO_MARKER = frozenset({"", " ", "None", "none"})


@dataclass(frozen=True)
class Style:
    """Everything about a figure's look except colour.

    Attributes
    ----------
    rc : dict
        rcParams applied (scoped) while rendering and saving.
    panel_size : tuple[float, float]
        Inches per grid cell ``(width, height)``.
    surface, ink : str
        Background and spine/tick colour.
    ticks : dict
        ``major_length``, ``minor_length``, ``major_width``, ``minor_width``.
    exp_ticks : bool
        Label log axes by exponent and prefix the axis label with
        ``log10``.
    inside_ticks : bool
        Draw y tick labels inside the panel, right of the left spine.
    gap : float
        Gap between neighbouring panels, in inches (rows and columns
        alike); ``0`` makes panels touch.
    cycle : tuple
        Style sequence (linestyles, then markers).
    lw : float
        Default line / marker-edge width.
    marker : dict
        ``size`` (pt), ``gap_pt`` (on-screen arc-length spacing of a
        marker-only series) and ``line_gap_pt`` (sparser spacing of the
        markers on a line+marker series, e.g. ``"-o"``).
    legend : dict
        Keyword arguments for ``Axes.legend`` / ``Figure.legend``.
    solo_color : str
        Colour of an unlabelled series alone on its panel (no explicit
        colour); a colour already used by the shared legend is skipped.
    legend_max : int
        Most entries a legend may hold; a larger legend is not drawn, and
        the builders split more curves than this over several panels.
    tag : dict
        ``fontsize`` and box ``alpha`` of inside panel titles and tags.
    save : dict
        ``dpi`` (raster formats) and ``pad_inches``.
    ref_cols : int
        Columns the style is sized for; a figure with more used columns has
        every visual size scaled by ``columns / ref_cols`` (:func:`scaled`).
    scale : float
        Derived, not user-settable: the factor :func:`scaled` applied.
    """

    rc: dict
    panel_size: tuple[float, float]
    surface: str
    ink: str
    ticks: dict
    exp_ticks: bool
    inside_ticks: bool
    gap: float
    cycle: tuple
    lw: float
    marker: dict
    legend: dict
    legend_max: int
    solo_color: str
    tag: dict
    save: dict
    ref_cols: int
    scale: float = 1.0


@dataclass(frozen=True)
class Palette:
    """Colours: a ``cycle`` plus a label -> colour ``mapping`` (may be empty)."""

    cycle: tuple[str, ...]
    mapping: dict


_DEFAULT_CYCLE: tuple = (
    "-",
    (0, (6, 3)),
    (0, (1, 2)),
    (0, (3, 1, 1, 1)),
    "o",
    "X",
    "P",
    "D",
    "s",
    "^",
)

#: Look of a ``role="reference"`` series: colour, width factor, draw order.
REFERENCE_COLOR: str = "#4d4d4d"
REFERENCE_LW_FACTOR: float = 1.9
REFERENCE_ZORDER: float = 1.9

_BASE_RC: dict[str, Any] = {
    "font.family": "DejaVu Sans",
    "mathtext.fontset": "stix",
    "svg.fonttype": "none",
    "lines.solid_capstyle": "round",
    "lines.dash_capstyle": "round",
    "figure.facecolor": "#ffffff",
}

_LEGEND: dict[str, Any] = {
    "frameon": False,
    "handlelength": 2.2,
    "labelspacing": 0.35,
    "title_fontproperties": {"weight": "bold"},
}

PAPER = Style(
    rc={
        **_BASE_RC,
        "font.size": 9,
        "axes.labelsize": 11,
        "axes.labelpad": 3,
        "axes.linewidth": 0.6,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "xtick.major.pad": 2.5,
        "ytick.major.pad": 2.5,
        "legend.fontsize": 8.5,
    },
    panel_size=(1.9, 1.6),
    surface="#ffffff",
    ink="#000000",
    ticks={
        "major_length": 3,
        "minor_length": 1.5,
        "major_width": 0.6,
        "minor_width": 0.5,
    },
    exp_ticks=True,
    inside_ticks=True,
    gap=0.08,
    cycle=_DEFAULT_CYCLE,
    lw=0.9,
    marker={"size": 3.0, "gap_pt": 8.0, "line_gap_pt": 24.0},
    legend=copy.deepcopy(_LEGEND),
    legend_max=12,
    solo_color="#00798c",
    tag={"fontsize": 9.5, "alpha": 0.85},
    save={"dpi": 400, "pad_inches": 0.03},
    ref_cols=4,
)

TALK = Style(
    rc={
        **_BASE_RC,
        "font.size": 20,
        "axes.titlesize": 23,
        "axes.labelsize": 23,
        "axes.labelpad": 12,
        "axes.linewidth": 0.8,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
        "legend.fontsize": 17,
    },
    panel_size=(5.0, 4.0),
    surface="#ffffff",
    ink="#000000",
    ticks={
        "major_length": 8,
        "minor_length": 4,
        "major_width": 0.8,
        "minor_width": 0.8,
    },
    exp_ticks=False,
    inside_ticks=False,
    gap=0.25,
    cycle=_DEFAULT_CYCLE,
    lw=2.1,
    marker={"size": 9.0, "gap_pt": 28.0, "line_gap_pt": 70.0},
    legend=copy.deepcopy(_LEGEND),
    legend_max=12,
    solo_color="#00798c",
    tag={"fontsize": 20, "alpha": 0.85},
    save={"dpi": 210, "pad_inches": 0.08},
    ref_cols=3,
)

#: Registered styles by name (built-ins plus :func:`register_style` entries).
STYLES: dict[str, Style] = {"paper": PAPER, "talk": TALK}
"""Registered styles by name (built-ins plus :func:`register_style` entries).

Built-ins: ``"paper"`` (default; print size, inside y tick labels, exponent
log axes) and ``"talk"`` (large fonts and ticks, thick lines).
"""

#: Registered palettes by name (built-ins plus :func:`register_palette`).
PALETTES: dict[str, Palette] = {
    "okabe-ito": Palette(
        (
            "#0072b2",
            "#d55e00",
            "#009e73",
            "#cc79a7",
            "#e69f00",
            "#56b4e9",
            "#f0e442",
            "#000000",
        ),
        {},
    ),
    "sunset": Palette(("#1d3557", "#6d597a", "#b56576", "#eaac8b"), {}),
    "vivid": Palette(("#d1495b", "#2e86ab", "#3c9d4e", "#9b5de5", "#4d4d4d"), {}),
    "jaff": Palette(
        (
            "#0072b2",
            "#d55e00",
            "#009e73",
            "#cc79a7",
            "#e69f00",
            "#56b4e9",
            "#000000",
            "#e6194b",
            "#4363d8",
            "#469990",
            "#9a6324",
            "#808000",
            "#000075",
            "#f032e6",
        ),
        {},
    ),
    "logo": Palette(("#8b6cff", "#e05fb0", "#ff6a5a", "#ffc24b"), {}),
    "gow": Palette(
        (
            "#e6194b",
            "#3cb44b",
            "#f58231",
            "#4363d8",
            "#911eb4",
            "#42d4f4",
            "#f032e6",
            "#9a6324",
            "#808000",
            "#008080",
            "#000075",
            "#e6a817",
            "#bfef45",
            "#fabed4",
            "#469990",
        ),
        {},
    ),
    "muted": Palette(
        (
            "#4878d0",
            "#ee854a",
            "#6acc64",
            "#d65f5f",
            "#956cb4",
            "#8c613c",
            "#dc7ec0",
            "#797979",
            "#d5bb67",
            "#82c6e2",
        ),
        {},
    ),
    "deep": Palette(
        (
            "#4C72B0",
            "#DD8452",
            "#55A868",
            "#C44E52",
            "#8172B3",
            "#937860",
            "#DA8BC3",
            "#8C8C8C",
            "#CCB974",
            "#64B5CD",
        ),
        {},
    ),
}
"""Registered palettes by name (built-ins plus :func:`register_palette`).

Built-ins: ``"jaff"`` (default), ``"okabe-ito"``, ``"gow"``, ``"logo"``,
``"sunset"``, ``"vivid"``, ``"muted"`` and ``"deep"``.
"""

DEFAULT_STYLE: str = "paper"
DEFAULT_PALETTE: str = "jaff"

_BUILTIN_STYLES = frozenset(STYLES)
_BUILTIN_PALETTES = frozenset(PALETTES)
_MERGED_FIELDS = frozenset({"rc", "ticks", "marker", "legend", "tag", "save"})
_STYLE_FIELDS = frozenset(f.name for f in fields(Style)) - {"scale"}
#: rcParams whose (numeric) values are sizes in points, scaled by :func:`scaled`.
_RC_SIZE_SUFFIXES = ("size", "width", "pad")


def _panel_size(value: Any) -> tuple[float, float]:
    """Coerce *value* to a pair of positive floats."""
    try:
        w, h = (float(v) for v in value)
    except (TypeError, ValueError):
        raise ValueError(f"style.panel_size: expected (width, height), got {value!r}")
    if w <= 0 or h <= 0:
        raise ValueError(f"style.panel_size: must be positive, got {value!r}")
    return (w, h)


def _style_from_dict(d: dict) -> Style:
    """Layer a partial style dict over its ``base`` (default ``"paper"``)."""
    base = resolve_style(d.get("base", DEFAULT_STYLE))
    unknown = sorted(set(d) - _STYLE_FIELDS - {"base"})
    if unknown:
        raise ValueError(
            f"style: unknown field(s) {unknown}; allowed: {sorted(_STYLE_FIELDS)}"
        )
    changes: dict[str, Any] = {}
    for key, value in d.items():
        if key == "base":
            continue
        if key in _MERGED_FIELDS:
            if not isinstance(value, dict):
                raise ValueError(f"style.{key}: expected a dict, got {value!r}")
            changes[key] = {**getattr(base, key), **value}
        elif key == "cycle":
            if isinstance(value, str) or not isinstance(value, (list, tuple)):
                raise ValueError(f"style.cycle: expected a sequence, got {value!r}")
            if not value:
                raise ValueError("style.cycle: must not be empty")
            changes[key] = tuple(value)
        elif key == "panel_size":
            changes[key] = _panel_size(value)
        elif key == "gap":
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or value < 0
            ):
                raise ValueError(f"style.gap: expected inches >= 0, got {value!r}")
            changes[key] = float(value)
        elif key == "solo_color":
            if not isinstance(value, str) or not is_color_like(value):
                raise ValueError(f"style.solo_color: expected a colour, got {value!r}")
            changes[key] = value
        elif key == "ref_cols":
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(
                    f"style.ref_cols: expected a positive int, got {value!r}"
                )
            changes[key] = value
        elif key == "legend_max":
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(
                    f"style.legend_max: expected a positive int, got {value!r}"
                )
            changes[key] = value
        else:
            changes[key] = value
    return replace(base, **changes)


def scaled(style: Style, k: float) -> Style:
    """A copy of *style* with every visual size multiplied by *k*.

    Scales point sizes in ``rc`` (fonts, line widths, pads), ``lw``, tick
    geometry, marker size and arc gaps, and the panel ``gap``.  Panel titles
    and tags (``tag``) keep their size, and ``panel_size`` stays, so a wide figure shrunk to the page reads like a
    ``ref_cols``-column one.  ``k == 1`` returns *style* itself; registered
    styles are never modified.
    """
    if k == 1:
        return style
    rc = {
        key: value * k
        if key.endswith(_RC_SIZE_SUFFIXES)
        and isinstance(value, (int, float))
        and not isinstance(value, bool)
        else value
        for key, value in style.rc.items()
    }
    marker = {
        key: value * k if key in ("size", "gap_pt", "line_gap_pt") else value
        for key, value in style.marker.items()
    }
    return replace(
        style,
        rc=rc,
        lw=style.lw * k,
        gap=style.gap * k,
        ticks={key: value * k for key, value in style.ticks.items()},
        marker=marker,
        scale=style.scale * k,
    )


def resolve_style(value: str | dict | Style | None) -> Style:
    """Return the :class:`Style` for a name, partial dict, instance or ``None``.

    Raises
    ------
    ValueError
        For an unknown name or unknown dict field.
    """
    if value is None:
        return STYLES[DEFAULT_STYLE]
    if isinstance(value, Style):
        return value
    if isinstance(value, dict):
        return _style_from_dict(value)
    if not isinstance(value, str):
        raise ValueError(
            f"style: expected a name, dict or Style, got {type(value).__name__}"
        )
    if value not in STYLES:
        raise ValueError(f"style: unknown style {value!r}; available: {sorted(STYLES)}")
    return STYLES[value]


def resolve_palette(value: str | list | tuple | dict | Palette | None) -> Palette:
    """Return the :class:`Palette` for a name, colour list, label dict or ``None``.

    A dict maps labels to colours; unmapped labels use the default cycle.

    Raises
    ------
    ValueError
        For an unknown name or an empty colour list.
    """
    if value is None:
        return PALETTES[DEFAULT_PALETTE]
    if isinstance(value, Palette):
        return value
    if isinstance(value, dict):
        return Palette(PALETTES[DEFAULT_PALETTE].cycle, dict(value))
    if isinstance(value, (list, tuple)):
        if not value:
            raise ValueError("palette: colour list is empty")
        return Palette(tuple(value), {})
    if not isinstance(value, str):
        raise ValueError(
            f"palette: expected a name, list, dict or Palette, got {type(value).__name__}"
        )
    if value not in PALETTES:
        raise ValueError(
            f"palette: unknown palette {value!r}; available: {sorted(PALETTES)}"
        )
    return PALETTES[value]


def register_style(name: str, style: dict) -> None:
    """Register a named style built from a partial dict (see :func:`resolve_style`).

    The new name can then be passed as ``style=`` to :func:`~jaff.plotting.plot`,
    :func:`~jaff.plotting.plot_rates` and :func:`~jaff.plotting.plot_xsecs`.  The
    same partial dict can also be passed as ``style=`` directly, without
    registering it.

    Parameters
    ----------
    name : str
        Name of the new style; registering an existing non-built-in name
        replaces it.
    style : dict
        Partial override of a base style.  ``"base"`` names the style to
        start from (default ``"paper"``).  The nested fields ``rc``,
        ``ticks``, ``marker``, ``legend``, ``tag`` and ``save`` are merged
        key-wise into the base's; every other field replaces the base value.
        Fields:

        - ``rc`` -- rcParams applied (scoped) while rendering and saving.
        - ``panel_size`` -- inches per panel box ``(width, height)``.
        - ``gap`` -- inches between panels, rows and columns alike; ``0``
          makes panels touch.
        - ``surface``, ``ink`` -- background and spine/tick colour.
        - ``ticks`` -- ``major_length``, ``minor_length``, ``major_width``,
          ``minor_width``.
        - ``exp_ticks`` -- label log axes by exponent, with at most ~4
          labelled decades, minor ticks on decades (at most 4 between
          labels, no 2-9 sub-ticks) and fractional exponents (``1.05``,
          ``1.1``, ``−1.9``) on an axis spanning less than 1.5 decades.
        - ``inside_ticks`` -- draw y tick labels inside the panel.
        - ``cycle`` -- line-style / marker sequence (non-empty).
        - ``lw`` -- default line / marker-edge width.
        - ``marker`` -- ``size`` (pt), ``gap_pt`` (spacing of a marker-only
          series) and ``line_gap_pt`` (sparser spacing on a line+marker
          series such as ``"-o"``).
        - ``legend`` -- keyword arguments for the matplotlib legends.
        - ``legend_max`` -- positive int (12 in both built-ins): a legend
          with more entries is not drawn, and :func:`~jaff.plotting.plot_rates`
          / :func:`~jaff.plotting.plot_xsecs` without a ``grid`` split more
          curves than this over an automatic grid.
        - ``solo_color`` -- colour of an unlabelled series alone on its
          panel (``"#00798c"``, teal, in both built-ins); when the shared
          legend already uses it, the next unused palette colour is taken.
        - ``tag`` -- ``fontsize`` and box ``alpha`` of inside panel titles
          and tags.
        - ``save`` -- ``dpi`` (raster formats) and ``pad_inches``.
        - ``ref_cols`` -- positive int, the column count the style is sized
          for (``"paper"`` 4, ``"talk"`` 3).  A figure with more used columns
          (after empty ones collapse) has every visual size (fonts, line
          widths, ticks, markers, gaps, pads) scaled by
          ``columns / ref_cols``, except panel titles and tags, while
          ``panel_size`` stays: a 6-column paper figure shrunk to the page
          reads like a 4-column one.  An explicit ``figsize`` turns the
          scaling off.

    Raises
    ------
    ValueError
        If *name* is a built-in style, or *style* has an unknown field or an
        invalid value.

    Notes
    -----
    Built-in styles (``"paper"``, ``"talk"``) cannot be overwritten, and
    registered styles are never modified by rendering.

    Examples
    --------
    >>> from jaff.plotting import STYLES, register_style
    >>> register_style("poster", {"base": "talk", "rc": {"font.size": 26}})
    >>> sorted(STYLES)
    ['paper', 'poster', 'talk']

    A one-off override needs no registration, e.g. fused panels:

    >>> from jaff.plotting import plot
    >>> spec = {"axes": [{"series": [{"x": [1, 2, 3], "y": [1, 4, 9]}]}]}
    >>> fig, axes = plot(spec | {"style": {"base": "talk", "gap": 0}})
    """
    if name in _BUILTIN_STYLES:
        raise ValueError(f"style: cannot overwrite built-in style {name!r}")
    STYLES[name] = _style_from_dict(style)


def register_palette(name: str, colours: list | tuple | dict) -> None:
    """Register a named palette from a colour list or label -> colour dict.

    The new name can then be passed as ``palette=`` to
    :func:`~jaff.plotting.plot`, :func:`~jaff.plotting.plot_rates` and
    :func:`~jaff.plotting.plot_xsecs`.  A list or dict can also be passed as
    ``palette=`` directly, without registering it.

    Parameters
    ----------
    name : str
        Name of the new palette; registering an existing non-built-in name
        replaces it.
    colours : list, tuple or dict
        A colour cycle, or a dict mapping series **labels** to colours
        (other labels use the default ``"jaff"`` cycle).

    Raises
    ------
    ValueError
        If *name* is a built-in palette or *colours* is empty.

    Examples
    --------
    >>> from jaff.plotting import register_palette
    >>> register_palette("sources", {"Gong+2017": "#4d4d4d", "sim": "#d55e00"})
    >>> register_palette("duo", ["#0072b2", "#d55e00"])

    Use it by name; the label-mapped colour wins over the cycle:

    >>> from jaff.plotting import plot
    >>> spec = {"axes": [{"series": [{"x": [1, 2], "y": [1, 2], "label": "sim"}]}]}
    >>> fig, axes = plot(spec | {"palette": "sources"})
    >>> axes[0, 0].lines[0].get_color()
    '#d55e00'
    """
    if name in _BUILTIN_PALETTES:
        raise ValueError(f"palette: cannot overwrite built-in palette {name!r}")
    PALETTES[name] = resolve_palette(
        list(colours) if isinstance(colours, tuple) else colours
    )


def is_marker(style: Any) -> bool:
    """``True`` when *style* names a marker rather than a linestyle."""
    return (
        isinstance(style, str)
        and style not in _LINESTYLES
        and style not in _NO_MARKER
        and style in Line2D.markers
    )


_LINE_PREFIXES = ("--", "-.", "-", ":")


def split_style(style: Any) -> tuple[Any, str | None]:
    """Split a series style into ``(linestyle, marker)``.

    ``"-o"`` / ``"--D"`` (matplotlib fmt strings) give both; a plain
    linestyle or dash tuple gives ``(style, None)``; a marker gives
    ``(None, marker)``.
    """
    if not isinstance(style, str):
        return style, None
    if style in _LINESTYLES:
        return style, None
    if is_marker(style):
        return None, style
    for prefix in _LINE_PREFIXES:
        if style.startswith(prefix) and is_marker(style[len(prefix) :]):
            return prefix, style[len(prefix) :]
    return style, None


def label_looks(labels: list[str], style: Style, palette: Palette) -> dict[str, tuple]:
    """One ``(colour, style)`` per distinct label, in first-seen order.

    Used when one shared legend keys several panels: a label must look the
    same everywhere.  Label *k* gets ``palette.cycle[k % nc]`` (or its
    ``mapping`` colour) and ``style.cycle[min(k // nc, len - 1)]``.
    """
    nc = len(palette.cycle)
    if nc == 0:
        raise ValueError("palette: colour cycle is empty")
    looks: dict[str, tuple] = {}
    for label in labels:
        if label and label not in looks:
            k = len(looks)
            color = palette.mapping.get(label, palette.cycle[k % nc])
            looks[label] = (color, style.cycle[min(k // nc, len(style.cycle) - 1)])
    return looks


def assign(
    series: list,
    style: Style,
    palette: Palette,
    by_label: dict[str, tuple] | None = None,
) -> list[tuple[str, Any]]:
    """Resolve ``(colour, style)`` for each series on one axes.

    Explicit ``color`` / ``style`` win; then *by_label* (figure-wide looks of a
    shared legend, see :func:`label_looks`); then the palette ``mapping`` by
    label (colour only); then the cycles, colour first: the k-th series still
    needing a colour gets ``palette.cycle[k % nc]``, the j-th still needing a
    style gets ``style.cycle[min(j // nc, len - 1)]``.  With *by_label* the
    cycles start after the looks it hands out, so a series outside the shared
    legend (e.g. an unlabelled one) never looks like a legend entry.
    Reference series sit outside all of this (dark grey, solid).
    """
    nc = len(palette.cycle)
    if nc == 0:
        raise ValueError("palette: colour cycle is empty")
    by_label = by_label or {}
    ci = si = len(by_label)
    out: list[tuple[str, Any]] = []
    for s in series:
        if getattr(s, "role", "data") == "reference":
            # Outside both cycles: grey, solid unless told otherwise.
            out.append((s.color or REFERENCE_COLOR, s.style or "-"))
            continue
        fixed_color, fixed_style = by_label.get(s.label, (None, None))
        color = s.color or fixed_color or palette.mapping.get(s.label)
        if color is None:
            color = palette.cycle[ci % nc]
            ci += 1
        ls = s.style
        if s.kind == "bars":
            ls = "-" if ls is None else ls  # bars ignore style; keep the cycle intact
        elif ls is None and fixed_style is not None:
            ls = fixed_style
        elif ls is None:
            ls = style.cycle[min(si // nc, len(style.cycle) - 1)]
            si += 1
        out.append((color, ls))
    return out
