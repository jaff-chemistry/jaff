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
