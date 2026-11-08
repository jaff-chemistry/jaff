# ABOUTME: Values legend of rate-coloured network figures: number formatting, label
# ABOUTME: templates, groups from net.analysis and the boxed legend under the bar
"""The boxed legend of the analysis values a rate-coloured figure used."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Any

from matplotlib.axes import Axes
from matplotlib.patches import FancyBboxPatch
from sympy import Expr, Float, Symbol, latex

from ...core.analysis._evaluate import _function
from ._layout import (
    EQUATION_FONT,
    HEADER_FONT,
    MARGIN,
    VALUE_GAP,
    VALUE_PAD,
    VALUE_ROW,
    VALUE_SPACE,
    VALUE_TITLE,
    ValueBox,
    ValueGroup,
)
from ._themes import NetworkTheme

if TYPE_CHECKING:
    from ...core.analysis import Analysis

#: Title of each table's group, in drawing order.
TITLES = {
    "parameters": "PARAMETERS",
    "densities": r"NUMBER DENSITIES (cm$^{-3}$)",
    "functions": "FUNCTIONS",
}
#: Significant figures of a value.
FIGURES = 3


def value_tex(value: float) -> str:
    """*value* as mathtext without dollars: 3 significant figures.

    Plain in ``[1e-2, 1e3)``, else ``m\\times10^{n}`` (``10^{n}`` when ``m`` is 1).
    """
    if value == 0:
        return "0"
    sign = "-" if value < 0 else ""
    size = abs(value)
    exponent = math.floor(math.log10(size))
    mantissa = float(f"{size / 10**exponent:.{FIGURES - 1}f}")
    if mantissa >= 10:
        mantissa, exponent = mantissa / 10, exponent + 1
    if -2 <= exponent < 3:
        return f"{sign}{float(f'{size:.{FIGURES}g}'):g}"
    if mantissa == 1:
        return f"{sign}10^{{{exponent}}}"

    return rf"{sign}{mantissa:g}\times10^{{{exponent}}}"


def entry_text(label: str, tex: str) -> str:
    """A legend entry: *label* with the value *tex* at its ``{}``.

    Inside maths the value goes in bare, outside it in dollars; a label without
    ``{}`` reads ``label = value``, in one maths span when the label is one.
    """
    if "{}" not in label:
        if label.startswith("$") and label.endswith("$") and label.count("$") == 2:
            return f"{label[:-1]} = {tex}$"
        return f"{label} = ${tex}$"
    at = label.index("{}")
    inside = label[:at].replace(r"\$", "").count("$") % 2 == 1

    return label[:at] + (tex if inside else f"${tex}$") + label[at + 2 :]


def _function_tex(name: str, value: Any, parameters: dict[str, float]) -> str:
    """A function value: its number for the parameters, else its LaTeX."""
    expr: Expr = _function(name, value)
    numeric = expr.xreplace(
        {
            s: Float(parameters[s.name])
            for s in expr.free_symbols
            if isinstance(s, Symbol) and s.name in parameters
        }
    )
    if numeric.is_number:
        return value_tex(float(numeric.evalf()))

    return latex(expr)


def value_groups(analysis: Analysis) -> tuple[ValueGroup, ...]:
    """The labelled entries of *analysis*, one group per table that has any."""
    if not analysis.enabled:
        return ()
    params, dens, funcs = analysis.tables()
    numbers = params.values()
    groups = []
    for key, table in (("parameters", params), ("densities", dens), ("functions", funcs)):
        entries = []
        for name, entry in table.get().items():
            if entry.label is None or entry.value is None:
                continue
            if key == "functions":
                tex = _function_tex(name, entry.value, numbers)
            else:
                tex = value_tex(entry.value)
            entries.append(entry_text(entry.label, tex))
        if entries:
            groups.append(ValueGroup(TITLES[key], tuple(entries)))

    return tuple(groups)


def draw_value_box(
    ax: Axes,
    size: tuple[float, float],
    theme: NetworkTheme,
    box: ValueBox,
    groups: tuple[ValueGroup, ...],
    bottom: float,
) -> None:
    """The boxed values legend, centred, its lower edge at *bottom* (y down)."""
    width = size[0]
    x0 = max(MARGIN, (width - box.width) / 2)
    top = bottom - box.height + VALUE_SPACE
    ax.add_patch(
        FancyBboxPatch(
            (x0, top),
            box.width,
            bottom - top,
            boxstyle="round,pad=0,rounding_size=0.08",
            facecolor=theme.surface,
            edgecolor=theme.edge,
            lw=0.8,
            zorder=2,
            gid="values-box",
        )
    )
    y = top + VALUE_PAD
    left = x0 + VALUE_PAD
    for group in groups:
        ax.text(
            left,
            y + VALUE_TITLE / 2,
            group.title,
            ha="left",
            va="center",
            fontsize=HEADER_FONT,
            fontweight="bold",
            color=theme.muted,
            zorder=3,
        )
        y += VALUE_TITLE
        for k, entry in enumerate(group.entries):
            row, column = divmod(k, box.columns)
            ax.text(
                left + column * box.column_width,
                y + (row + 0.5) * VALUE_ROW,
                entry,
                ha="left",
                va="center",
                fontsize=EQUATION_FONT,
                color=theme.text,
                zorder=3,
            )
        y += math.ceil(len(group.entries) / box.columns) * VALUE_ROW + VALUE_GAP
