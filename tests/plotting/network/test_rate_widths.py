# ABOUTME: Tests for rate-scaled arrows in net.visualize: log-normalised widths
# ABOUTME: between a min and max, the 12-decade floor and the legend width key

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.legend import Legend

from jaff.plotting.network import NetworkVisualizer
from jaff.plotting.network._widths import DECADES, edge_widths, width_key

LO, HI = 0.4, 3.5

# H3+ formed by 0, destroyed by 1-3; reaction 4 (CR) does not involve H3+.
LINES = [
    "H2 + H2+ -> H + H3+ [10,1000] 2e-9",
    "H3+ + e- -> H + H2 [10,1000] 2e-8*(tgas/300)**(-0.5)",
    "H3+ + O -> H2 + OH+ [10,1000] 8e-10",
    "H3+ + CO -> H2 + HCO+ [10,1000] 2e-11",
    "H2 -> H2+ + e- [] 1.2*crate",
]


@pytest.fixture
def net(make_network):
    return make_network(LINES)


# --------------------------------------------------------------------------- #
# Normalisation                                                                #
# --------------------------------------------------------------------------- #


def test_largest_gets_max_and_smallest_gets_min():
    widths = edge_widths(np.array([1e-9, 1e-12, 1e-15]), LO, HI)

    assert widths == pytest.approx([HI, (LO + HI) / 2, LO])


def test_widths_are_logarithmic():
    widths = edge_widths(np.array([1e-8, 1e-9, 1e-10, 1e-12]), LO, HI)
    step = (HI - LO) / 4

    assert widths == pytest.approx([HI, HI - step, HI - 2 * step, LO])


def test_values_below_the_floor_get_the_min_width():
    values = np.array([1.0, 10.0**-DECADES, 1e-30])
    widths = edge_widths(values, LO, HI)

    assert widths == pytest.approx([HI, LO, LO])


def test_zero_negative_and_non_finite_get_the_min_width():
    widths = edge_widths(np.array([1e-9, 0.0, -1.0, np.nan, np.inf]), LO, HI)

    assert widths == pytest.approx([HI, LO, LO, LO, LO])


def test_equal_values_get_the_max_width():
    assert edge_widths(np.array([3e-10, 3e-10]), LO, HI) == pytest.approx([HI, HI])


def test_no_positive_value_gives_min_widths():
    assert edge_widths(np.array([0.0, np.nan]), LO, HI) == pytest.approx([LO, LO])


# --------------------------------------------------------------------------- #
# Legend key                                                                   #
# --------------------------------------------------------------------------- #


def test_key_has_at_most_three_entries():
    key = width_key(np.array([1e-8, 1e-11, 1e-15]), LO, HI)

    assert 1 <= len(key) <= 3


def test_key_values_are_powers_of_ten_in_range_widest_first():
    key = width_key(np.array([2e-8, 3e-15]), LO, HI)
    values = [v for v, _ in key]
    widths = [w for _, w in key]

    assert all(np.log10(v) == pytest.approx(round(np.log10(v))) for v in values)
    assert all(3e-15 <= v <= 2e-8 for v in values)
    assert widths == sorted(widths, reverse=True)
    assert widths == pytest.approx(
        list(edge_widths(np.array([2e-8, 3e-15, *values]), LO, HI)[2:])
    )


def test_key_for_one_value():
    key = width_key(np.array([5e-10]), LO, HI)

    assert len(key) == 1 and key[0][1] == pytest.approx(HI)


# --------------------------------------------------------------------------- #
# Figures                                                                      #
# --------------------------------------------------------------------------- #


def _edges(ax):
    """Edge line widths by (role, reaction index)."""
    out = {}
    for p in ax.patches:
        gid = p.get_gid() or ""
        if gid.startswith("edge:"):
            _, role, index = gid.split(":")
            out.setdefault((role, int(index)), set()).add(p.get_linewidth())

    return out


def _one(widths: set) -> float:
    """The single width of a reaction's edges."""
    (width,) = widths

    return width


def _legends(ax):
    """Every legend on *ax*: the type legend and, when scaled, the width key."""
    return [a for a in ax.get_children() if isinstance(a, Legend)]


def _key(ax):
    """Labels of the width-key entries, in drawing order."""
    return [
        t.get_text()
        for legend in _legends(ax)
        for t in legend.get_texts()
        if "10^" in t.get_text()
    ]


def test_graph_without_values_keeps_the_fixed_widths(net):
    _, ax = NetworkVisualizer(net).graph(show=False)

    widths = _edges(ax)
    assert all(w == {0.8} for (role, _), w in widths.items() if role == "reactant")
    assert all(w == {1.1} for (role, _), w in widths.items() if role == "product")
    assert _key(ax) == []


def test_graph_scales_the_product_arrows(net):
    values = {"tgas": 300.0, "crate": 1e-17}
    _, ax = NetworkVisualizer(net).graph(show=False, rate_values=values)
    widths = _edges(ax)
    k = net.reactions.evaluate_rate_coefficients(values)
    expected = edge_widths(k, LO, HI)

    for (role, i), w in widths.items():
        if role == "product":
            assert _one(w) == pytest.approx(expected[i])
    assert 1 <= len(_key(ax)) <= 3


def test_reactant_lines_keep_the_fixed_width(net):
    values = {"tgas": 300.0, "crate": 1e-17}
    _, ax = NetworkVisualizer(net).graph(show=False, rate_values=values)

    reactant = [w for (role, _), w in _edges(ax).items() if role == "reactant"]
    assert reactant and all(w == {0.8} for w in reactant)


def test_width_key_is_its_own_row_widest_first(net):
    values = {"tgas": 300.0, "crate": 1e-17}
    _, ax = NetworkVisualizer(net).graph(show=False, rate_values=values)
    key = [g for g in _legends(ax) if any("10^" in t.get_text() for t in g.get_texts())]

    (legend,) = key
    lines = legend.legend_handles
    assert all("10^" in t.get_text() for t in legend.get_texts())
    assert legend._ncols == len(lines)
    widths = [h.get_linewidth() for h in lines]
    assert widths == sorted(widths, reverse=True)


def test_width_key_lies_inside_the_figure_below_the_types(net):
    values = {"tgas": 300.0, "crate": 1e-17}
    fig, ax = NetworkVisualizer(net).graph(show=False, close=False, rate_values=values)
    fig.canvas.draw()
    boxes = {
        "10^" in g.get_texts()[0].get_text(): g.get_window_extent() for g in _legends(ax)
    }
    plt.close(fig)

    key, types = boxes[True], boxes[False]
    assert key.y1 <= types.y0 + 1  # display y points up: key below the types
    assert key.y0 >= 0


def test_species_view_normalises_over_its_own_reactions(net):
    # crate is only used by reaction 4, which the H3+ view does not show.
    _, ax = NetworkVisualizer(net).species("H3+", show=False, rate_values={"tgas": 300.0})
    widths = _edges(ax)

    # Shown: 2e-9, 2e-8, 8e-10, 2e-11 -> three decades between min and max.
    assert _one(widths[("product", 1)]) == pytest.approx(HI)
    assert _one(widths[("product", 3)]) == pytest.approx(LO)
    assert _one(widths[("product", 0)]) == pytest.approx(LO + (HI - LO) * 2 / 3)


def test_arrow_width_sets_the_limits(net):
    values = {"tgas": 300.0, "crate": 1e-17}
    _, ax = NetworkVisualizer(net).graph(
        show=False, rate_values=values, arrow_width=(0.2, 5.0)
    )
    all_widths = {w for ws in _edges(ax).values() for w in ws}

    assert max(all_widths) == pytest.approx(5.0)
    assert min(all_widths) >= 0.2


def test_missing_value_raises(net):
    with pytest.raises(ValueError, match="crate"):
        NetworkVisualizer(net).graph(show=False, rate_values={"tgas": 300.0})


def test_rate_functions_without_values_raises(net):
    with pytest.raises(ValueError, match="rate_values"):
        NetworkVisualizer(net).graph(show=False, rate_functions={"f": 1.0})


def test_bad_arrow_width_raises(net):
    with pytest.raises(ValueError, match="arrow_width"):
        NetworkVisualizer(net).graph(
            show=False, rate_values={"tgas": 1.0, "crate": 1.0}, arrow_width=(3.0, 1.0)
        )


def test_no_figure_left_open(net):
    NetworkVisualizer(net).graph(show=False, rate_values={"tgas": 300.0, "crate": 1e-17})

    assert plt.get_fignums() == []


def test_tight_crop_keeps_the_width_key(net):
    values = {"tgas": 300.0, "crate": 1e-17}
    fig, ax = NetworkVisualizer(net).graph(
        show=False, close=False, tight=True, rate_values=values
    )
    fig.canvas.draw()
    figure_box = fig.bbox
    boxes = [g.get_window_extent() for g in _legends(ax)]
    plt.close(fig)

    assert len(boxes) == 2
    for box in boxes:
        assert box.y0 >= figure_box.y0 - 1 and box.y1 <= figure_box.y1 + 1
