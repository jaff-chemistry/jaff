# ABOUTME: Colour themes of the network figures (net.visualize): NetworkTheme,
# ABOUTME: the NETWORK_THEMES registry and reaction-type colour assignment
"""Colour themes for :class:`~jaff.plotting.NetworkVisualizer` figures."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

OKABE_ITO = ("#0072b2", "#d55e00", "#009e73", "#cc79a7", "#e69f00", "#56b4e9", "#f0e442")
SLATE = ("#7b8794", "#3e4c59")


@dataclass(frozen=True)
class NetworkTheme:
    """Colours of a network figure.

    Attributes
    ----------
    surface : str
        Figure background and species-pill fill.
    band : str
        Background band behind each atom-count column.
    ink : str
        Species outlines and names.
    text : str
        Equations and legend text.
    muted : str
        Column headers.
    halo : str
        Outline behind equations, so edges do not cut through them.
    edge : str
        Reactant edges (and product edges when ``typed_edges`` is off).
    typed_edges : bool
        Draw product edges in their reaction-type colour.
    typed_arrows : bool
        Draw product arrowheads in their reaction-type colour.
    types : tuple of str
        Reaction-type colours, given out in order of first appearance.
    other : str
        Colour of ``unknown`` reactions and of types beyond ``types``.
    consumed : tuple of str
        Cell fill for a consumed species: count 1, count 2 or more.
    produced : tuple of str
        Cell fill for a produced species: count 1, count 2 or more.
    produced_outline : bool
        Draw produced cells as an outline (in ``outline``) instead of a fill.
    outline : str
        Outline colour of produced cells when ``produced_outline`` is on.
    """

    surface: str
    band: str
    ink: str
    text: str
    muted: str
    halo: str
    edge: str
    typed_edges: bool
    typed_arrows: bool
    types: tuple[str, ...]
    other: str
    consumed: tuple[str, str]
    produced: tuple[str, str]
    produced_outline: bool
    outline: str


_LIGHT = dict(
    surface="#ffffff",
    band="#f5f6f8",
    ink="#1d2025",
    text="#33373d",
    muted="#8a8f98",
    halo="#ffffff",
    edge="#9aa0a8",
    other="#8a8f98",
    outline="#3e4c59",
)

NETWORK_THEMES: dict[str, NetworkTheme] = {
    "vivid": NetworkTheme(
        **_LIGHT,
        typed_edges=True,
        typed_arrows=True,
        types=(
            "#2e86ab",
            "#d1495b",
            "#3c9d4e",
            "#9b5de5",
            "#f18f01",
            "#00a6a6",
            "#c2185b",
            "#6b8e23",
            "#1d3557",
            "#b5651d",
        ),
        consumed=SLATE,
        produced=("#ffffff", "#ffffff"),
        produced_outline=True,
    ),
    "okabe": NetworkTheme(
        **_LIGHT,
        typed_edges=True,
        typed_arrows=True,
        types=OKABE_ITO,
        consumed=("#f4a582", "#d6604d"),
        produced=("#92c5de", "#4393c3"),
        produced_outline=False,
    ),
    "neutral": NetworkTheme(
        **_LIGHT,
        typed_edges=False,
        typed_arrows=True,
        types=OKABE_ITO,
        consumed=SLATE,
        produced=("#ffffff", "#ffffff"),
        produced_outline=True,
    ),
    "muted": NetworkTheme(
        **_LIGHT,
        typed_edges=True,
        typed_arrows=True,
        types=(
            "#5b7fa6",
            "#b5835a",
            "#6f9a7d",
            "#9c7fb0",
            "#c9a227",
            "#7a9eb1",
            "#b07a8c",
        ),
        consumed=("#c3b5d9", "#7a5ca6"),
        produced=("#f6c27a", "#e08a1e"),
        produced_outline=False,
    ),
    "print": NetworkTheme(
        **_LIGHT,
        typed_edges=False,
        typed_arrows=False,
        types=("#1d2025", "#6b7280", "#a7adb5"),
        consumed=SLATE,
        produced=("#ffffff", "#ffffff"),
        produced_outline=True,
    ),
    "dark": NetworkTheme(
        surface="#1e2228",
        band="#272c33",
        ink="#e6edf3",
        text="#c9d1d9",
        muted="#8b949e",
        halo="#1e2228",
        edge="#6e7681",
        typed_edges=True,
        typed_arrows=True,
        types=(
            "#56b4e9",
            "#e69f00",
            "#2bc49b",
            "#e78ac3",
            "#f0e442",
            "#ff7f50",
            "#a6d854",
        ),
        other="#8b949e",
        consumed=("#8b949e", "#c9d1d9"),
        produced=("#1e2228", "#1e2228"),
        produced_outline=True,
        outline="#c9d1d9",
    ),
}
"""Registered themes by name; ``"vivid"`` is the default."""


def resolve_theme(theme: str | NetworkTheme) -> NetworkTheme:
    """The theme named *theme*, or *theme* itself.

    Raises
    ------
    ValueError
        If *theme* is a name that is not registered.
    """
    if isinstance(theme, NetworkTheme):
        return theme
    if theme not in NETWORK_THEMES:
        raise ValueError(
            f"Unknown network theme {theme!r}; choose one of {sorted(NETWORK_THEMES)}"
        )

    return NETWORK_THEMES[theme]


def _ordered(types: Sequence[str]) -> list[str]:
    """Distinct types in order of first appearance, ``unknown`` left out."""
    return [t for t in dict.fromkeys(types) if t != "unknown"]


def type_colours(theme: NetworkTheme, types: Sequence[str]) -> dict[str, str]:
    """Colour of each reaction type in *types* (one entry per distinct type)."""
    colours = {
        t: theme.types[k] if k < len(theme.types) else theme.other
        for k, t in enumerate(_ordered(types))
    }
    if "unknown" in types:
        colours["unknown"] = theme.other

    return colours


def legend_types(theme: NetworkTheme, types: Sequence[str]) -> list[tuple[str, str]]:
    """Legend entries ``(label, colour)``: one per coloured type, then ``other``."""
    ordered = _ordered(types)
    named = ordered[: len(theme.types)]
    entries = [(t.replace("_", " "), theme.types[k]) for k, t in enumerate(named)]
    if "unknown" in types or len(ordered) > len(named):
        entries.append(("other", theme.other))

    return entries
