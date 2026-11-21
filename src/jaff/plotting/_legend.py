# ABOUTME: Legend policy: one shared figure legend, or per-panel groups (opt-in)
# ABOUTME: Creates the legends; the layout module positions the margin ones

"""Legend placement.

Figure policy ``"auto"`` (and ``"figure"``): ONE shared, de-duplicated
legend for the whole figure, in the bottom-right hidden grid cell when the
grid is unfinished, else in the right margin.  ``"axes"``: one legend group
per panel, collected in the right margin and headed by the panel's title, tag
or letter (``"(a)"``, ``"(b)"``, ...; an untitled panel is tagged with its
letter inside).  ``"none"``: no legend.  An axes-level ``legend`` (``"axes"``
or ``"none"``) overrides the figure policy for that panel.

A legend with more de-duplicated entries than ``style.legend_max`` is not
drawn (a warning names the axes, or ``figure``); a shared legend over the cap
that keys several panels falls back to per-panel groups instead
(:func:`fallback`).  An optional ``legend_title`` heads the shared legend, or
the stack of groups.  Margin legends are created here and positioned by
:mod:`jaff.plotting._layout`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from matplotlib.artist import Artist
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.legend import Legend
from matplotlib.text import Text

from ..io import JaffLogger
from .spec import AxesSpec
from .themes import Style

#: ``(label, legend handle)`` pairs for one axes, in draw order.
Entries = list[tuple[str, Artist]]

#: Per-panel legend mode resolved from the figure policy and axes overrides.
SHARED, GROUP, NONE = "shared", "group", "none"


@dataclass
class Placement:
    """The legends drawn, for the layout to position.

    Attributes
    ----------
    shared : Legend or None
        The shared figure legend.
    host : int or None
        Flat grid index of the hidden cell holding ``shared``; ``None`` when
        it belongs in the right margin.
    groups : list[Legend]
        Per-panel legends (``"axes"`` mode), stacked in the right margin.
    heading : Text or None
        Bold ``legend_title`` over the stack of groups.
    """

    shared: Legend | None = None
    host: int | None = None
    groups: list[Legend] = field(default_factory=list)
    heading: Text | None = None

    def margin(self) -> list[Legend | Text]:
        """Artists that live in the right margin, in stacking order."""
        lead = [self.shared] if self.shared is not None and self.host is None else []
        head = [self.heading] if self.heading is not None else []
        return head + lead + self.groups


def letter(index: int) -> str:
    """Panel letter ``"(a)"``, ``"(b)"``, ... for a 0-based panel index."""
    out, n = "", index
    while True:
        out = chr(ord("a") + n % 26) + out
        n = n // 26 - 1
        if n < 0:
            return f"({out})"


def modes(specs: tuple[AxesSpec, ...], policy: str) -> list[str]:
    """Legend mode per panel: ``shared``, ``group`` or ``none``."""
    default = {"auto": SHARED, "figure": SHARED, "axes": GROUP}.get(policy, NONE)
    override = {"axes": GROUP, "none": NONE}
    return [override.get(s.legend, default) if s.legend else default for s in specs]


def fallback(
    specs: tuple[AxesSpec, ...], panel_modes: list[str], style: Style
) -> list[str]:
    """Per-panel groups when a shared legend over several panels is too long.

    Counts the distinct labels of the shared panels; above ``legend_max``
    (and with more than one shared panel) every shared panel switches to its
    own group, with one warning.
    """
    shared = [i for i, m in enumerate(panel_modes) if m == SHARED]
    labels = {s.label for i in shared for s in specs[i].series if s.label}
    if len(shared) < 2 or len(labels) <= style.legend_max:
        return panel_modes
    JaffLogger().get_logger().warning(
        f"figure: shared legend with {len(labels)} entries exceeds legend_max="
        f"{style.legend_max}; using per-panel legends."
    )
    return [GROUP if m == SHARED else m for m in panel_modes]


def heading(spec: AxesSpec, index: int) -> str:
    """Group heading: the panel title, else its tag, else its letter."""
    return spec.title or spec.tag or letter(index)


def add_tag(ax: Axes, text: str, style: Style, drop_pt: float = 0.0) -> None:
    """Write *text* inside the upper-right corner of *ax* on a soft box.

    *drop_pt* lowers it by that many points (to sit under an inside title),
    independent of the panel's size.
    """
    ax.annotate(
        text,
        xy=(0.95, 0.94),
        xycoords="axes fraction",
        xytext=(0.0, -drop_pt),
        textcoords="offset points",
        ha="right",
        va="top",
        fontsize=style.tag["fontsize"],
        bbox=box(style),
        zorder=10,
    )


def box(style: Style) -> dict:
    """The soft surface box behind inside panel text."""
    return dict(
        boxstyle="round,pad=0.15", fc=style.surface, ec="none", alpha=style.tag["alpha"]
    )


def _dedupe(entries: list[Entries]) -> dict[str, Artist]:
    """Map label to handle, the first handle winning."""
    seen: dict[str, Artist] = {}
    for e in entries:
        for label, handle in e:
            seen.setdefault(label, handle)
    return seen


def _over_cap(seen: dict[str, Artist], where: str, style: Style) -> bool:
    """Whether a legend for *seen* exceeds ``style.legend_max`` (warns if so)."""
    if len(seen) <= style.legend_max:
        return False
    JaffLogger().get_logger().warning(
        f"{where}: legend with {len(seen)} entries exceeds legend_max="
        f"{style.legend_max}; not drawn."
    )
    return True


def _kw(style: Style) -> dict[str, Any]:
    """Legend keywords: the style's, upper-left anchored, left-aligned title."""
    return {
        **style.legend,
        "loc": "upper left",
        "borderaxespad": 0.0,
        "alignment": "left",
    }


def apply(
    fig: Figure,
    used: list[Axes],
    free: list[Axes],
    specs: tuple[AxesSpec, ...],
    entries: list[Entries],
    panel_modes: list[str],
    style: Style,
    title: str | None = None,
) -> Placement:
    """Create the shared legend and per-panel group legends.

    *title* heads the shared legend, else the stack of groups.  Margin
    legends are anchored provisionally; :func:`jaff.plotting._layout.apply`
    moves them into place.
    """
    placed = Placement()
    shared = _dedupe([e for e, m in zip(entries, panel_modes) if m == SHARED])
    if shared and not _over_cap(shared, "figure", style):
        handles, labels = list(shared.values()), list(shared)
        kw: dict[str, Any] = {**_kw(style), "title": title}
        if free:
            placed.host = len(used) + len(free) - 1
            placed.shared = free[-1].legend(
                handles, labels, bbox_to_anchor=(0.0, 1.0), **kw
            )
        else:
            placed.shared = fig.legend(handles, labels, bbox_to_anchor=(1.0, 1.0), **kw)
    for i, (ax, spec, ent, mode) in enumerate(zip(used, specs, entries, panel_modes)):
        if mode != GROUP or not ent:
            continue
        seen = _dedupe([ent])
        if _over_cap(seen, f"axes[{i}]", style):
            continue
        placed.groups.append(
            ax.legend(
                list(seen.values()),
                list(seen),
                title=heading(spec, i),
                bbox_to_anchor=(1.0, 1.0),
                bbox_transform=fig.transFigure,
                **_kw(style),
            )
        )
    if title and placed.groups and placed.shared is None:
        placed.heading = fig.text(
            0.0,
            1.0,
            title,
            ha="left",
            va="top",
            fontweight="bold",
            fontsize=placed.groups[0].get_title().get_fontsize(),
        )
    return placed
