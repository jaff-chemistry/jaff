# ABOUTME: Tests for rate-coloured reactions in net.visualize: colour_by k or R from
# ABOUTME: net.analysis, log positions on a colormap, the colour bar, the errors

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib import colormaps
from matplotlib.colors import Colormap, to_hex
from matplotlib.patches import FancyArrowPatch

from jaff.plotting.network import NETWORK_THEMES, NetworkVisualizer
from jaff.plotting.network._rates import (
    DECADES,
    decade_ticks,
    rate_colormap,
    rate_positions,
)

# H3+ formed by 0, destroyed by 1-3; reaction 4 (CR) does not involve H3+.
LINES = [
    "H2 + H2+ -> H + H3+ [10,1000] 2e-9",
    "H3+ + e- -> H + H2 [10,1000] 2e-8*(tgas/300)**(-0.5)",
    "H3+ + O -> H2 + OH+ [10,1000] 8e-10",
    "H3+ + CO -> H2 + HCO+ [10,1000] 2e-11",
    "H2 -> H2+ + e- [] 1.2*crate",
]
VALUES = {"tgas": 300.0, "crate": 1e-17}
K = "rate_coefficient"


@pytest.fixture
def plain(make_network):
    return make_network(LINES)


@pytest.fixture
def net(plain):
    plain.analysis.enable(parameters=VALUES)

    return plain


# --------------------------------------------------------------------------- #
# Positions on the colormap                                                    #
# --------------------------------------------------------------------------- #


def test_largest_is_one_and_smallest_is_zero():
    pos, _ = rate_positions(np.array([1e-9, 1e-12, 1e-15]))

    assert pos == pytest.approx([1.0, 0.5, 0.0])


def test_positions_are_logarithmic():
    pos, _ = rate_positions(np.array([1e-8, 1e-9, 1e-10, 1e-12]))

    assert pos == pytest.approx([1.0, 0.75, 0.5, 0.0])


def test_values_below_the_floor_are_zero():
    pos, (low, top) = rate_positions(np.array([1.0, 10.0**-DECADES, 1e-30]))

    assert pos == pytest.approx([1.0, 0.0, 0.0])
    assert (low, top) == pytest.approx((-DECADES, 0.0))


def test_zero_negative_and_non_finite_are_zero():
    pos, _ = rate_positions(np.array([1e-9, 3e-11, 0.0, -1.0, np.nan, np.inf]))

    assert pos == pytest.approx([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])


def test_equal_values_are_one():
    pos, _ = rate_positions(np.array([3e-10, 3e-10]))

    assert pos == pytest.approx([1.0, 1.0])


def test_no_positive_value_has_no_range():
    pos, span = rate_positions(np.array([0.0, np.nan]))

    assert pos == pytest.approx([0.0, 0.0])
    assert span is None


# --------------------------------------------------------------------------- #
# Ticks and colormaps                                                          #
# --------------------------------------------------------------------------- #


def test_ticks_are_decades_in_range_at_most_six():
    ticks = decade_ticks(-17.8, -5.8)

    assert 1 <= len(ticks) <= 6
    assert all(isinstance(t, int) and -17.8 <= t <= -5.8 for t in ticks)
    assert ticks == sorted(ticks)
    assert ticks[-1] == -6


def test_narrow_range_without_a_decade_has_no_ticks():
    assert decade_ticks(-9.4, -9.1) == []


@pytest.mark.parametrize("name", list(NETWORK_THEMES))
def test_every_theme_has_a_rate_colormap(name):
    assert isinstance(rate_colormap(NETWORK_THEMES[name]), Colormap)


def test_vivid_uses_trimmed_ylgnbu():
    cmap = rate_colormap(NETWORK_THEMES["vivid"])
    base = colormaps["YlGnBu"]

    assert to_hex(cmap(0.0)) == to_hex(base(0.3))
    assert to_hex(cmap(1.0)) == to_hex(base(1.0))


def test_a_cmap_override_is_used_as_given():
    cmap = rate_colormap(NETWORK_THEMES["vivid"], "magma")

    assert to_hex(cmap(0.0)) == to_hex(colormaps["magma"](0.0))


# --------------------------------------------------------------------------- #
# Figures                                                                      #
# --------------------------------------------------------------------------- #


def _patches(ax, prefix):
    return {p.get_gid(): p for p in ax.patches if (p.get_gid() or "").startswith(prefix)}


def _expected_colours(net, indices, cmap, rates=False):
    exprs = [
        net.reactions[i].rate.expanded if rates else net.reactions[i].rate.coefficient
        for i in indices
    ]
    k = net.analysis.evaluate(exprs)
    pos, _ = rate_positions(k)

    return {i: to_hex(cmap(p)) for i, p in zip(indices, pos)}


def test_without_colour_by_the_type_colours_stay(net):
    _, ax = NetworkVisualizer(net).graph(show=False)

    assert len(ax.images) == 0
    edges = _patches(ax, "edge:reactant:")
    assert all(p.get_linewidth() == pytest.approx(0.8) for p in edges.values())


def test_graph_colours_markers_and_product_arrows(net):
    _, ax = NetworkVisualizer(net).graph(show=False, colour_by=K)
    cmap = rate_colormap(NETWORK_THEMES["vivid"])
    expected = _expected_colours(net, range(len(net.reactions)), cmap)
    theme = NETWORK_THEMES["vivid"]

    for gid, marker in _patches(ax, "reaction:").items():
        i = int(gid.split(":")[1])
        assert to_hex(marker.get_facecolor()) == expected[i]
        assert to_hex(marker.get_edgecolor()) == theme.ink
    for gid, edge in _patches(ax, "edge:product:").items():
        assert to_hex(edge.get_edgecolor()) == expected[int(gid.split(":")[2])]


def test_reactant_lines_stay_grey_and_thin(net):
    _, ax = NetworkVisualizer(net).graph(show=False, colour_by=K)
    theme = NETWORK_THEMES["vivid"]

    for edge in _patches(ax, "edge:reactant:").values():
        assert to_hex(edge.get_edgecolor()) == theme.edge
        assert edge.get_linewidth() == pytest.approx(0.8)


def test_arrowheads_take_the_rate_colour(net):
    _, ax = NetworkVisualizer(net).graph(show=False, colour_by=K)
    lines = {to_hex(p.get_edgecolor()) for p in _patches(ax, "edge:product:").values()}
    heads = {
        to_hex(p.get_facecolor()) for p in ax.patches if isinstance(p, FancyArrowPatch)
    }

    assert heads <= lines


def test_species_view_scales_over_its_own_reactions(plain):
    # crate is only used by reaction 4, which the H3+ view does not show.
    plain.analysis.enable(parameters={"tgas": 300.0})
    _, ax = NetworkVisualizer(plain).species("H3+", show=False, colour_by=K)
    cmap = rate_colormap(NETWORK_THEMES["vivid"])
    markers = _patches(ax, "reaction:")

    # Shown: 2e-9, 2e-8, 8e-10, 2e-11.
    assert to_hex(markers["reaction:1"].get_facecolor()) == to_hex(cmap(1.0))
    assert to_hex(markers["reaction:3"].get_facecolor()) == to_hex(cmap(0.0))


def test_rate_cmap_overrides_the_theme(net):
    _, ax = NetworkVisualizer(net).graph(show=False, colour_by=K, rate_cmap="magma")
    expected = _expected_colours(net, range(len(net.reactions)), colormaps["magma"])

    for gid, marker in _patches(ax, "reaction:").items():
        assert to_hex(marker.get_facecolor()) == expected[int(gid.split(":")[1])]


def test_type_legend_shows_shapes_only(net):
    _, ax = NetworkVisualizer(net).graph(show=False, colour_by=K)
    theme = NETWORK_THEMES["vivid"]
    markers = [
        h for h in ax.get_legend().legend_handles if h.get_linestyle() in ("None", "none")
    ]

    assert markers
    assert all(to_hex(h.get_markerfacecolor()) == theme.surface for h in markers)


# --------------------------------------------------------------------------- #
# Reaction rates                                                               #
# --------------------------------------------------------------------------- #


def test_colour_by_rate_uses_the_reaction_rates(net):
    dens = np.linspace(1.0, 2.0, len(net.species))
    net.analysis.enable(densities=dens)
    _, ax = NetworkVisualizer(net).graph(show=False, colour_by="rate")
    expected = _expected_colours(
        net, range(len(net.reactions)), rate_colormap(NETWORK_THEMES["vivid"]), True
    )

    for gid, marker in _patches(ax, "reaction:").items():
        assert to_hex(marker.get_facecolor()) == expected[int(gid.split(":")[1])]


def test_colour_bar_title_follows_the_quantity(net):
    net.analysis.enable(densities=np.ones(len(net.species)))
    _, k_ax = NetworkVisualizer(net).graph(show=False, colour_by=K)
    _, r_ax = NetworkVisualizer(net).graph(show=False, colour_by="rate")

    assert "RATE COEFFICIENT $k$" in [t.get_text() for t in k_ax.texts]
    assert any(
        t.get_text().startswith(r"REACTION RATE $\mathcal{R}$") for t in r_ax.texts
    )


def test_colour_by_rate_needs_densities(net):
    with pytest.raises(ValueError, match="no density"):
        NetworkVisualizer(net).graph(show=False, colour_by="rate")


# --------------------------------------------------------------------------- #
# Colour bar                                                                   #
# --------------------------------------------------------------------------- #


def test_colour_bar_is_a_rounded_gradient_with_decade_labels(net):
    _, ax = NetworkVisualizer(net).graph(show=False, colour_by=K)

    (image,) = ax.images
    assert image.get_clip_path() is not None
    texts = [t.get_text() for t in ax.texts]
    assert "RATE COEFFICIENT $k$" in texts
    ticks = [t for t in texts if t.startswith("$10^{")]
    assert 1 <= len(ticks) <= 6


def test_colour_bar_lies_above_the_legend_inside_the_figure(net):
    fig, ax = NetworkVisualizer(net).graph(show=False, close=False, colour_by=K)
    fig.canvas.draw()
    bar = ax.images[0].get_window_extent()
    legend = ax.get_legend().get_window_extent()
    plt.close(fig)

    assert bar.y0 > legend.y1  # display y points up
    assert 0 < bar.x0 and bar.x1 < fig.bbox.width


def test_tight_crop_keeps_the_colour_bar(net):
    fig, ax = NetworkVisualizer(net).graph(
        show=False, close=False, tight=True, colour_by=K
    )
    fig.canvas.draw()
    labels = [t for t in ax.texts if t.get_text().startswith("$10^{")]
    boxes = [t.get_window_extent() for t in labels]
    boxes.append(ax.images[0].get_window_extent())
    plt.close(fig)

    for box in boxes:
        assert box.y0 >= -1 and box.x0 >= -1 and box.x1 <= fig.bbox.width + 1


# --------------------------------------------------------------------------- #
# Errors                                                                       #
# --------------------------------------------------------------------------- #


def test_missing_value_raises(plain):
    plain.analysis.enable(parameters={"tgas": 300.0})

    with pytest.raises(ValueError, match="crate"):
        NetworkVisualizer(plain).graph(show=False, colour_by=K)


def test_colour_by_needs_analysis_mode(plain):
    with pytest.raises(ValueError, match="analysis"):
        NetworkVisualizer(plain).graph(show=False, colour_by=K)


def test_unknown_colour_by_raises(net):
    with pytest.raises(ValueError, match="colour_by"):
        NetworkVisualizer(net).graph(show=False, colour_by="speed")


def test_rate_cmap_without_colour_by_raises(net):
    with pytest.raises(ValueError, match="colour_by"):
        NetworkVisualizer(net).graph(show=False, rate_cmap="magma")


def test_no_figure_left_open(net):
    NetworkVisualizer(net).graph(show=False, colour_by=K)

    assert plt.get_fignums() == []
