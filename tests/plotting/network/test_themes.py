# ABOUTME: Tests for network-visualisation colour themes: registry, resolution and
# ABOUTME: reaction-type colour assignment (first appearance, unknown, overflow)

from dataclasses import fields

import pytest
from matplotlib.colors import is_color_like

from jaff.plotting.network._themes import (
    NETWORK_THEMES,
    NetworkTheme,
    legend_types,
    resolve_theme,
    type_colours,
)

NAMES = ["vivid", "okabe", "neutral", "muted", "print", "dark"]


def test_registered_themes():
    assert list(NETWORK_THEMES) == NAMES


@pytest.mark.parametrize("name", NAMES)
def test_every_colour_is_valid(name):
    theme = NETWORK_THEMES[name]
    for f in fields(theme):
        value = getattr(theme, f.name)
        colours = value if isinstance(value, tuple) else (value,)
        for colour in colours:
            if isinstance(colour, str):
                assert is_color_like(colour), (name, f.name, colour)


def test_resolve_by_name_and_object():
    assert resolve_theme("vivid") is NETWORK_THEMES["vivid"]
    custom = NETWORK_THEMES["print"]
    assert resolve_theme(custom) is custom


def test_unknown_theme_lists_choices():
    with pytest.raises(ValueError, match="vivid"):
        resolve_theme("nope")


def test_vivid_has_ten_type_colours():
    assert len(NETWORK_THEMES["vivid"].types) == 10


def test_type_colours_by_first_appearance():
    theme = NETWORK_THEMES["vivid"]
    colours = type_colours(theme, ["photo", "two_body", "photo", "unknown"])

    assert colours == {
        "photo": theme.types[0],
        "two_body": theme.types[1],
        "unknown": theme.other,
    }


def test_overflow_types_share_other():
    theme = NetworkTheme(**{**NETWORK_THEMES["vivid"].__dict__, "types": ("#000000",)})
    colours = type_colours(theme, ["a", "b", "c"])

    assert colours == {"a": "#000000", "b": theme.other, "c": theme.other}


def test_legend_types_collapse_other():
    theme = NetworkTheme(**{**NETWORK_THEMES["vivid"].__dict__, "types": ("#000000",)})

    assert legend_types(theme, ["a", "unknown", "b", "c"]) == [
        ("a", "#000000"),
        ("other", theme.other),
    ]


def test_legend_labels_are_readable():
    theme = NETWORK_THEMES["vivid"]

    assert legend_types(theme, ["two_body"])[0][0] == "two body"


def _contrast(a: str, b: str) -> float:
    """WCAG contrast ratio of two colours."""
    from matplotlib.colors import to_rgb

    def luminance(colour: str) -> float:
        r, g, b = (
            c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
            for c in to_rgb(colour)
        )
        return 0.2126 * r + 0.7152 * g + 0.0722 * b

    high, low = sorted((luminance(a), luminance(b)), reverse=True)
    return (high + 0.05) / (low + 0.05)


@pytest.mark.parametrize("name", ["vivid", "neutral", "print"])
def test_slate_cell_counts_are_easy_to_read(name):
    from jaff.plotting.network._common import text_on

    theme = NETWORK_THEMES[name]
    for fill in theme.consumed:
        assert _contrast(fill, text_on(fill, theme)) >= 7, fill


@pytest.mark.parametrize("name", NAMES)
def test_other_marker_numbers_are_easy_to_read(name):
    from jaff.plotting.network._common import text_on

    theme = NETWORK_THEMES[name]

    assert _contrast(theme.other, text_on(theme.other, theme)) >= 7


@pytest.mark.parametrize("name", NAMES)
def test_text_on_picks_the_higher_contrast(name):
    from jaff.plotting.network._common import text_on

    theme = NETWORK_THEMES[name]
    for fill in [*theme.types, theme.other, *theme.consumed]:
        chosen = text_on(fill, theme)
        other = theme.ink if chosen == theme.surface else theme.surface
        assert _contrast(fill, chosen) >= _contrast(fill, other), fill
