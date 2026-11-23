# ABOUTME: Tests for jaff.plotting._legend — shared legend vs per-panel groups
# ABOUTME: Covers the auto policy, forced policies, overrides, caps and layout

import numpy as np

from jaff.plotting._legend import letter
from jaff.plotting import render
from jaff.plotting.spec import FigureSpec

X = np.linspace(1.0, 2.0, 20)


def _s(label):
    return {"x": X, "y": X, "label": label}


def _render(axes, **kw):
    return render(FigureSpec.from_dict({"axes": axes, **kw}))


def _texts(ax):
    return [t.get_text() for t in ax.texts]


def _legend_texts(legend):
    return [t.get_text() for t in legend.get_texts()]


def _drawn(fig):
    fig.canvas.draw()
    return fig.canvas.get_renderer()


def _overlap(a, b, tol=0.5):
    return (
        min(a.x1, b.x1) - max(a.x0, b.x0) > tol
        and min(a.y1, b.y1) - max(a.y0, b.y0) > tol
    )


# --- auto: one shared legend -------------------------------------------------


def test_single_panel_legend_in_right_margin():
    fig, axes = _render([{"series": [_s("a"), _s("b")]}])
    (legend,) = fig.legends
    assert _legend_texts(legend) == ["a", "b"] and axes[0, 0].get_legend() is None
    r = _drawn(fig)
    assert legend.get_window_extent(r).x0 >= axes[0, 0].get_window_extent(r).x1


def test_lone_series_gets_legend_not_tag():
    fig, axes = _render([{"series": [_s("CO")]}])
    assert _texts(axes[0, 0]) == [] and _legend_texts(fig.legends[0]) == ["CO"]


def test_unlabelled_series_no_legend_no_tag():
    fig, axes = _render([{"series": [_s(""), _s("")]}])
    assert not fig.legends and axes[0, 0].get_legend() is None
    assert _texts(axes[0, 0]) == []


def test_shared_legend_in_bottom_right_free_cell():
    panel = {"series": [_s("ref"), _s("sim")]}
    _, axes = _render([panel, panel, panel], grid=(2, 2))
    used, free = axes.ravel()[:3], axes.ravel()[3]
    assert all(a.get_legend() is None for a in used)
    assert _legend_texts(free.get_legend()) == ["ref", "sim"]


def test_shared_legend_clear_of_labels_hanging_into_its_cell():
    panel = {"series": [_s("ref"), _s("sim")], "xlabel": "x"}
    fig, axes = _render([panel, panel, panel], grid=(2, 2), sharex=True)
    r = _drawn(fig)
    legend = axes[1, 1].get_legend().get_window_extent(r)
    above = axes[0, 1]
    hanging = [t.get_window_extent(r) for t in above.get_xticklabels() if t.get_text()]
    assert hanging and all(legend.y1 <= t.y0 + 0.5 for t in hanging)


def test_full_grid_shared_legend_in_right_margin():
    panel = {"series": [_s("ref"), _s("sim")]}
    fig, axes = _render([panel, panel], grid=(1, 2))
    (legend,) = fig.legends
    r = _drawn(fig)
    assert legend.get_window_extent(r).x0 >= axes[0, 1].get_window_extent(r).x1
    assert legend.get_window_extent(r).x1 <= fig.bbox.x1 + 0.5


def test_shared_legend_dedupes_across_panels():
    fig, _ = _render([{"series": [_s("a")]}, {"series": [_s("a"), _s("b")]}], grid=(1, 2))
    (legend,) = fig.legends
    assert _legend_texts(legend) == ["a", "b"]


def test_policy_figure_same_as_auto():
    fig, _ = _render(
        [{"series": [_s("a")]}, {"series": [_s("a"), _s("b")]}],
        grid=(1, 2),
        legend="figure",
    )
    (legend,) = fig.legends
    assert _legend_texts(legend) == ["a", "b"]


def test_policy_none():
    fig, axes = _render([{"series": [_s("a"), _s("b")]}], legend="none")
    assert not fig.legends and axes[0, 0].get_legend() is None


def test_axes_level_none_excluded_from_shared():
    fig, _ = _render(
        [{"series": [_s("a")]}, {"series": [_s("b")], "legend": "none"}], grid=(1, 2)
    )
    assert _legend_texts(fig.legends[0]) == ["a"]


# --- tags ---------------------------------------------------------------------


def test_explicit_tag_drawn_inside():
    panel = {"series": [_s("ref"), _s("sim")], "tag": "CO"}
    _, axes = _render([panel, panel], grid=(1, 2))
    assert all(_texts(a) == ["CO"] for a in axes.ravel())


def test_empty_tag_draws_nothing():
    _, axes = _render([{"series": [_s("CO")], "tag": ""}])
    assert _texts(axes[0, 0]) == []


def test_tag_below_inside_title():
    fig, axes = _render([{"series": [_s("a")], "title": "CO", "tag": "n=1"}])
    r = _drawn(fig)
    ax = axes[0, 0]
    (tag,) = ax.texts
    assert tag.get_window_extent(r).y1 <= ax.title.get_window_extent(r).y0 + 0.5


# --- axes: per-panel groups in the right margin ------------------------------


def test_policy_axes_group_per_panel():
    _, axes = _render(
        [{"series": [_s("a")], "title": "CO"}, {"series": [_s("b"), _s("c")]}],
        grid=(1, 2),
        legend="axes",
    )
    left, right = axes.ravel()
    assert _legend_texts(left.get_legend()) == ["a"]
    assert left.get_legend().get_title().get_text() == "CO"
    assert right.get_legend().get_title().get_text() == "(b)"
    assert _texts(right) == ["(b)"] and _texts(left) == []


def test_axes_level_override_to_group():
    fig, axes = _render(
        [{"series": [_s("a")]}, {"series": [_s("b")], "legend": "axes"}], grid=(1, 2)
    )
    assert _legend_texts(fig.legends[0]) == ["a"]
    assert _legend_texts(axes[0, 1].get_legend()) == ["b"]


def test_group_legend_dedupes_labels():
    _, axes = _render([{"series": [_s("a"), _s("a")]}], legend="axes")
    assert _legend_texts(axes[0, 0].get_legend()) == ["a"]


def _long(n, panel):
    return [
        _s(f"panel {panel}: a fairly long reaction label number {k}") for k in range(n)
    ]


def test_groups_stack_without_overlap_and_wrap():
    fig, axes = _render(
        [{"series": _long(12, p)} for p in range(4)], grid=(2, 2), legend="axes"
    )
    r = _drawn(fig)
    boxes = [a.get_legend().get_window_extent(r) for a in axes.ravel()]
    panels = [a.get_window_extent(r) for a in axes.ravel()]
    for lb in boxes:
        assert not any(_overlap(lb, pb) for pb in panels)
        assert lb.x1 <= fig.bbox.x1 + 0.5 and lb.y0 >= -0.5
    for i, a in enumerate(boxes):
        assert not any(_overlap(a, b) for b in boxes[i + 1 :])


def test_letters():
    assert [letter(i) for i in (0, 1, 25, 26)] == ["(a)", "(b)", "(z)", "(aa)"]


# --- legend_max cap ------------------------------------------------------------


def _many(n, prefix="s"):
    return [_s(f"{prefix}{i}") for i in range(n)]


def test_shared_legend_over_cap_not_drawn(jaff_warnings):
    fig, _ = _render([{"series": _many(13)}])
    assert not fig.legends
    assert any("figure" in m and "legend" in m for m in jaff_warnings)


def test_shared_legend_at_cap_drawn(jaff_warnings):
    fig, _ = _render([{"series": _many(12)}])
    assert len(fig.legends[0].get_texts()) == 12
    assert not jaff_warnings


def test_cap_counts_deduplicated_entries():
    fig, _ = _render([{"series": _many(10) + _many(10)}])
    assert len(fig.legends[0].get_texts()) == 10


def test_cap_follows_style():
    fig, _ = _render([{"series": _many(13)}], style={"legend_max": 20})
    assert len(fig.legends[0].get_texts()) == 13


def test_group_over_cap_not_drawn(jaff_warnings):
    _, axes = _render([{"series": _many(13)}], legend="axes")
    assert axes[0, 0].get_legend() is None
    assert any("axes[0]" in m and "legend" in m for m in jaff_warnings)


def test_figure_legend_over_cap_across_panels_falls_back(jaff_warnings):
    # review item 9: per-panel groups instead of no legend at all
    fig, axes = _render(
        [{"series": _many(7, "a")}, {"series": _many(7, "b")}], grid=(1, 2)
    )
    assert not fig.legends
    assert [len(a.get_legend().get_texts()) for a in axes.ravel()] == [7, 7]
    assert any("figure" in m and "per-panel" in m for m in jaff_warnings)


# --- figure size ---------------------------------------------------------------


def test_explicit_figsize_respected():
    fig, _ = _render([{"series": _long(12, 0)}], figsize=(3, 2))
    assert tuple(fig.get_size_inches()) == (3, 2)


def test_figure_grows_for_margin_legend():
    with_legend, _ = _render([{"series": _long(4, 0)}])
    without, _ = _render([{"series": _long(4, 0)}], legend="none")
    assert with_legend.get_size_inches()[0] > without.get_size_inches()[0] + 1
    assert with_legend.get_size_inches()[1] == without.get_size_inches()[1]
