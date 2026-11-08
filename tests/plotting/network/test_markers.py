# ABOUTME: Tests for reaction-type marker shapes in net.visualize: assignment by
# ABOUTME: first appearance, "other" shape, drawing, legend handles and the views

import matplotlib.pyplot as plt
import pytest
from matplotlib.lines import Line2D

from jaff.plotting.network import NETWORK_THEMES, NetworkVisualizer
from jaff.plotting.network._markers import (
    MPL_MARKERS,
    OTHER_SHAPE,
    SHAPES,
    draw_marker,
    type_shapes,
)

TYPED = [
    "H2 -> H + H  []  PHOTO, 11.2",
    "H -> H+ + e-  []  1e-17*crate",
    "H + H + H -> H2 + H  []  1e-30",
    "H+ + e- -> H  []  1e-12",
]


# --------------------------------------------------------------------------- #
# Assignment                                                                   #
# --------------------------------------------------------------------------- #


def test_seven_distinct_shapes_and_an_other():
    assert len(SHAPES) == len(set(SHAPES)) == 7
    assert OTHER_SHAPE not in SHAPES
    assert set(MPL_MARKERS) == {*SHAPES, OTHER_SHAPE}


def test_shapes_by_first_appearance():
    shapes = type_shapes(["photo", "two_body", "photo", "cosmic_ray"])

    assert shapes == {"photo": SHAPES[0], "two_body": SHAPES[1], "cosmic_ray": SHAPES[2]}


def test_unknown_gets_the_other_shape():
    assert type_shapes(["unknown", "photo"]) == {
        "unknown": OTHER_SHAPE,
        "photo": SHAPES[0],
    }


def test_types_beyond_the_shapes_get_the_other_shape():
    types = [f"t{i}" for i in range(9)]
    shapes = type_shapes(types)

    assert [shapes[t] for t in types[:7]] == list(SHAPES)
    assert shapes["t7"] == shapes["t8"] == OTHER_SHAPE


# --------------------------------------------------------------------------- #
# Drawing                                                                      #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("shape", [*SHAPES, OTHER_SHAPE])
def test_index_sits_inside_every_shape(shape):
    fig = plt.figure(figsize=(2, 2))
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(0, 2)
    ax.set_ylim(2, 0)
    theme = NETWORK_THEMES["vivid"]
    patch = draw_marker(ax, 1.0, 1.0, shape, "#0072b2", theme, 12, gid="reaction:12")
    fig.canvas.draw()
    box = patch.get_window_extent()
    (text,) = ax.texts
    inner = text.get_window_extent()
    plt.close(fig)

    assert patch.get_gid() == "reaction:12"
    assert patch.get_label() == shape
    assert text.get_text() == "12"
    cx, cy = (inner.x0 + inner.x1) / 2, (inner.y0 + inner.y1) / 2
    assert box.x0 < cx < box.x1 and box.y0 < cy < box.y1


# --------------------------------------------------------------------------- #
# Views                                                                        #
# --------------------------------------------------------------------------- #


def _markers(ax):
    return {
        p.get_gid(): p for p in ax.patches if (p.get_gid() or "").startswith("reaction:")
    }


def _legend_markers(ax):
    return {
        h.get_label(): h.get_marker()
        for h in ax.get_legend().legend_handles
        if isinstance(h, Line2D) and h.get_linestyle() in ("None", "none", "")
    }


def test_graph_markers_follow_the_type(make_network):
    net = make_network(TYPED)
    _, ax = NetworkVisualizer(net).graph(show=False)
    shapes = type_shapes([r.type for r in net.reactions])

    for i, r in enumerate(net.reactions):
        assert _markers(ax)[f"reaction:{i}"].get_label() == shapes[r.type]


def test_graph_legend_shows_each_type_shape(make_network):
    net = make_network(TYPED)
    _, ax = NetworkVisualizer(net).graph(show=False)
    shapes = type_shapes([r.type for r in net.reactions])
    legend = _legend_markers(ax)

    for t, shape in shapes.items():
        label = "other" if shape == OTHER_SHAPE else t.replace("_", " ")
        assert legend[label] == MPL_MARKERS[shape]


def test_species_view_keeps_the_network_shapes(make_network):
    net = make_network(TYPED)
    _, ax = NetworkVisualizer(net).species("H", show=False)
    shapes = type_shapes([r.type for r in net.reactions])

    for gid, patch in _markers(ax).items():
        i = int(gid.split(":")[1])
        assert patch.get_label() == shapes[net.reactions[i].type]


@pytest.mark.parametrize("shape", [*SHAPES, OTHER_SHAPE])
def test_shapes_stay_clear_of_the_equation_above(shape):
    from jaff.plotting.network._layout import MARKER

    fig = plt.figure(figsize=(2, 2))
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(0, 2)
    ax.set_ylim(2, 0)
    patch = draw_marker(ax, 1.0, 1.0, shape, "#0072b2", NETWORK_THEMES["vivid"], 3)
    top = (
        ax.transData.inverted()
        .transform(patch.get_window_extent().corners())
        .min(axis=0)[1]
    )
    plt.close(fig)

    # The equation's baseline sits MARKER / 2 + 0.04 in above the centre.
    assert 1.0 - top <= MARKER / 2 + 0.04
