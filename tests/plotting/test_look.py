# ABOUTME: Tests for the spec §14 look: tight grid, inside ticks/titles, legends
# ABOUTME: Reference series, line+marker styles, palettes and decade tick rule

import matplotlib.colors as mcolors
import numpy as np
import pytest

from jaff.plotting import themes
from jaff.plotting.render import render
from jaff.plotting.spec import FigureSpec, SeriesSpec
from jaff.plotting.themes import resolve_style, split_style
from tests._paths import REPO_ROOT

X = np.logspace(17, 22, 120)
Y = 1e-6 * (1 + np.tanh(np.log10(X) - 19)) + 1e-12


def _render(axes, **kw):
    return render(FigureSpec.from_dict({"axes": axes, **kw}))


def _panel(label="a", **kw):
    return {
        "series": [{"x": X, "y": Y, "label": label}],
        "xscale": "log",
        "yscale": "log",
        **kw,
    }


def _inches(fig, box):
    return [v / fig.dpi for v in (box.x0, box.y0, box.x1, box.y1)]


def _draw(fig):
    fig.canvas.draw()
    return fig.canvas.get_renderer()


# --- 1. tight uniform grid ---------------------------------------------------


def test_equal_absolute_gaps():
    fig, axes = _render([_panel()] * 4, grid=(2, 2), sharex=True)
    _draw(fig)
    a, b, c = (
        _inches(fig, ax.get_position(original=False).transformed(fig.transFigure))
        for ax in (axes[0, 0], axes[0, 1], axes[1, 0])
    )
    col_gap, row_gap = b[0] - a[2], a[1] - c[3]
    assert col_gap == pytest.approx(row_gap, abs=0.01)
    assert col_gap == pytest.approx(resolve_style("paper").gap, abs=0.01)


def test_axes_have_panel_size():
    fig, axes = _render([_panel()] * 2, grid=(1, 2))
    _draw(fig)
    w, h = resolve_style("paper").panel_size
    for ax in axes.ravel():
        box = ax.get_window_extent()
        assert box.width / fig.dpi == pytest.approx(w, rel=0.02)
        assert box.height / fig.dpi == pytest.approx(h, rel=0.02)


def test_explicit_figsize_kept():
    fig, _ = _render([_panel()] * 2, grid=(1, 2), figsize=(5, 2.5))
    assert tuple(fig.get_size_inches()) == pytest.approx((5, 2.5))


def test_x_tick_label_sides_2x2():
    # §15 (supersedes §14.1 "lowest panel only"): top of a column -> top spine
    fig, axes = _render([_panel()] * 3, grid=(2, 2), sharex=True)
    _draw(fig)
    assert _xlabels(axes[0, 0], top=True) and not _xlabels(axes[0, 0])
    assert _xlabels(axes[1, 0]) and _xlabels(axes[0, 1])


def test_gap_zero_panels_touch():
    fig, axes = _render(
        [_panel()] * 2, grid=(1, 2), style={"base": "talk", "gap": 0}, sharey=True
    )
    _draw(fig)
    left, right = (ax.get_window_extent() for ax in axes.ravel())
    assert right.x0 - left.x1 == pytest.approx(0, abs=1)


@pytest.mark.parametrize("bad", [-0.1, "wide", True])
def test_gap_validation(bad):
    with pytest.raises(ValueError, match="style.gap"):
        resolve_style({"gap": bad})


def test_margins_tight_around_labels():
    fig, axes = _render([_panel(xlabel="N", ylabel="n")] * 2, grid=(1, 2))
    r = _draw(fig)
    tight = fig.get_tightbbox(r)
    w, h = fig.get_size_inches()
    assert tight.x0 >= -0.01 and tight.y0 >= -0.01
    assert tight.x1 <= w + 0.01 and tight.y1 <= h + 0.01
    assert tight.width >= w - 0.15 and tight.height >= h - 0.15


# --- 2. inside y tick labels -------------------------------------------------


def _ylabels(ax):
    return [t for t in ax.get_yticklabels() if t.get_text() and t.get_visible()]


def test_inside_y_tick_labels_paper():
    fig, axes = _render([_panel()])
    r = _draw(fig)
    ax = axes[0, 0]
    box = ax.get_window_extent(r)
    labels = _ylabels(ax)
    assert labels
    for t in labels:
        e = t.get_window_extent(r)
        assert e.x0 >= box.x0 and e.y0 >= box.y0 - 0.5 and e.y1 <= box.y1 + 0.5


def test_curves_start_right_of_inside_labels():
    fig, axes = _render([_panel()])
    r = _draw(fig)
    ax = axes[0, 0]
    widest = max(t.get_window_extent(r).x1 for t in _ylabels(ax))
    first = ax.transData.transform((X[0], Y[0]))[0]
    assert first >= widest


def test_talk_y_tick_labels_outside():
    fig, axes = _render([_panel()], style="talk")
    r = _draw(fig)
    ax = axes[0, 0]
    assert all(
        t.get_window_extent(r).x1 <= ax.get_window_extent(r).x0 + 0.5
        for t in _ylabels(ax)
    )


# --- 3. inside titles + headroom ---------------------------------------------


def test_short_title_inside_upper_right():
    fig, axes = _render([_panel(title="CO")])
    r = _draw(fig)
    ax = axes[0, 0]
    box, t = ax.get_window_extent(r), ax.title.get_window_extent(r)
    assert ax.get_title() == "CO"
    assert t.y1 <= box.y1 and t.x1 <= box.x1 and t.x0 > box.x0 + box.width / 2


def test_title_clear_of_data():
    fig, axes = _render([_panel(title="CO")])
    r = _draw(fig)
    ax = axes[0, 0]
    top = ax.transData.transform((X[-1], Y.max()))[1]
    assert top < ax.title.get_window_extent(r).y0


def test_long_title_goes_above():
    long = "a very long panel title that cannot fit inside the axes at all"
    fig, axes = _render([_panel(title=long)])
    r = _draw(fig)
    ax = axes[0, 0]
    assert ax.title.get_window_extent(r).y0 >= ax.get_window_extent(r).y1


# --- 4. legend policy --------------------------------------------------------


def _all_legend_texts(fig, axes):
    legends = list(fig.legends) + [a.get_legend() for a in axes.ravel() if a.get_legend()]
    return [[t.get_text() for t in lg.get_texts()] for lg in legends]


def test_auto_single_shared_legend_deduped():
    panels = [
        {"series": [{"x": X, "y": Y, "label": "ref"}, {"x": X, "y": Y, "label": "sim"}]},
        {"series": [{"x": X, "y": Y, "label": "ref"}, {"x": X, "y": Y, "label": "new"}]},
    ] * 2
    fig, axes = _render(panels, grid=(2, 2))
    assert _all_legend_texts(fig, axes) == [["ref", "sim", "new"]]


def test_auto_legend_in_free_bottom_right_cell():
    fig, axes = _render([_panel("a"), _panel("b"), _panel("c")], grid=(2, 2))
    host = axes[1, 1]
    assert host.get_legend() is not None
    assert [t.get_text() for t in host.get_legend().get_texts()] == ["a", "b", "c"]


def test_single_panel_legend_in_right_margin():
    fig, axes = _render(
        [
            {
                "series": [
                    {"x": X, "y": Y, "label": "a"},
                    {"x": X, "y": 2 * Y, "label": "b"},
                ],
                "xscale": "log",
                "yscale": "log",
            }
        ]
    )
    r = _draw(fig)
    (legend,) = fig.legends
    assert legend.get_window_extent(r).x0 >= axes[0, 0].get_window_extent(r).x1
    assert legend.get_window_extent(r).x1 <= fig.bbox.x1 + 0.5


def test_lone_label_gets_no_auto_tag():
    _, axes = _render([_panel("CO")])
    assert len(axes[0, 0].texts) == 0


def test_axes_mode_groups_in_right_margin():
    panels = [
        {
            "series": [{"x": X, "y": Y, "label": f"s{p}{k}"} for k in range(3)],
            "xscale": "log",
            "yscale": "log",
            **({"title": "CO"} if p == 0 else {}),
        }
        for p in range(3)
    ]
    fig, axes = _render(panels, grid=(1, 3), legend="axes")
    r = _draw(fig)
    right = max(a.get_window_extent(r).x1 for a in axes.ravel())
    legends = [a.get_legend() for a in axes.ravel()]
    assert all(lg is not None for lg in legends)
    assert [lg.get_title().get_text() for lg in legends] == ["CO", "(b)", "(c)"]
    boxes = [lg.get_window_extent(r) for lg in legends]
    assert all(b.x0 >= right for b in boxes)
    for i, a in enumerate(boxes):
        for b in boxes[i + 1 :]:
            assert (
                min(a.y1, b.y1) - max(a.y0, b.y0) <= 0.5
                or min(a.x1, b.x1) - max(a.x0, b.x0) <= 0.5
            )
    assert [t.get_text() for t in axes[0, 1].texts] == ["(b)"]
    assert all(b.x1 <= fig.bbox.x1 + 0.5 for b in boxes)


# --- 5. reference series -----------------------------------------------------


def test_reference_role_look():
    ref = {"x": X, "y": Y, "label": "Gong+2017", "role": "reference"}
    sim = {"x": X, "y": 2 * Y, "label": "sim"}
    _, axes = _render([{"series": [ref, sim], "xscale": "log", "yscale": "log"}])
    r_line, s_line = [ln for ln in axes[0, 0].lines if ln.get_visible()]
    style = resolve_style("paper")
    assert mcolors.to_hex(r_line.get_color()) == "#4d4d4d"
    assert r_line.get_linewidth() == pytest.approx(1.9 * style.lw)
    assert r_line.get_zorder() < s_line.get_zorder()
    assert mcolors.to_hex(s_line.get_color()) == themes.PALETTES["jaff"].cycle[0]


def test_role_validated():
    with pytest.raises(ValueError, match="role"):
        SeriesSpec.from_dict({"x": [1, 2], "y": [1, 2], "role": "ref"})


# --- 6. markers --------------------------------------------------------------


def test_marker_cycle():
    assert resolve_style("paper").cycle[-6:] == ("o", "X", "P", "D", "s", "^")


@pytest.mark.parametrize(
    ("style", "parts"),
    [
        ("-o", ("-", "o")),
        ("--D", ("--", "D")),
        ("o", (None, "o")),
        ("-", ("-", None)),
        ((0, (6, 3)), ((0, (6, 3)), None)),
    ],
)
def test_split_style(style, parts):
    assert split_style(style) == parts


def test_line_plus_marker_series():
    s = {"x": X, "y": Y, "style": "-o", "color": "#123456"}
    fig, axes = _render([{"series": [s], "xscale": "log", "yscale": "log"}])
    lines = [ln for ln in axes[0, 0].lines if ln.get_visible()]
    solid = [ln for ln in lines if ln.get_linestyle() == "-"]
    marks = [ln for ln in lines if ln.get_marker() == "o"]
    assert solid and marks
    m = marks[0]
    assert mcolors.to_hex(m.get_markerfacecolor()) == "#ffffff"
    assert 2 <= len(m.get_xdata()) < 20


# --- 7. palettes -------------------------------------------------------------


def test_default_palette_is_jaff():
    jaff = themes.PALETTES["jaff"].cycle
    assert themes.DEFAULT_PALETTE == "jaff" and len(jaff) == 14
    assert jaff[:3] == ("#0072b2", "#d55e00", "#009e73") and "#f0e442" not in jaff
    assert themes.PALETTES["logo"].cycle == ("#8b6cff", "#e05fb0", "#ff6a5a", "#ffc24b")
    assert len(themes.PALETTES["gow"].cycle) == 15


# --- 8. decade ticks ---------------------------------------------------------


def test_wide_log_axis_at_most_four_labels():
    y = np.logspace(-20, 0, 50)
    fig, axes = _render(
        [{"series": [{"x": X[:50], "y": y}], "xscale": "log", "yscale": "log"}]
    )
    _draw(fig)
    ax = axes[0, 0]
    labelled = [t for t in ax.get_yticklabels() if t.get_text()]
    assert 2 <= len(labelled) <= 4
    lo, hi = sorted(ax.get_ylim())
    minors = [v for v in ax.yaxis.get_minorticklocs() if lo <= v <= hi]
    exps = np.log10(minors)
    assert np.allclose(exps, np.round(exps))  # one minor per decade, no 2..9


def test_three_decades_labelled_every_decade():
    # inside labels off: their clear strip would widen the x range
    x = np.logspace(0, 3, 30)
    fig, axes = _render(
        [{"series": [{"x": x, "y": x}], "xscale": "log"}], style={"inside_ticks": False}
    )
    _draw(fig)
    labels = [t.get_text() for t in axes[0, 0].get_xticklabels() if t.get_text()]
    assert labels == ["0", "1", "2", "3"]


@pytest.mark.parametrize(
    ("lo", "hi", "major", "every"),
    [
        (1e-3, 1e1, [-2, 0], 1),
        (1e-21, 10, [-20, -10, 0], 2),  # step 10: minors every 2 decades
        (1, 1e3, [0, 1, 2, 3], 1),
        (1e-115, 1e-5, [-100, -75, -50, -25], 5),  # step 25: every 5
    ],
)
def test_decade_ticks_rule(lo, hi, major, every):
    # minors thinned for large steps (<= 4 between majors), never under a decade
    from jaff.plotting._ticks import decade_ticks

    labelled, minor = decade_ticks(lo, hi)
    assert labelled == major
    k0, k1 = round(np.log10(lo)), round(np.log10(hi))
    assert minor == [k for k in range(k0, k1 + 1) if k % every == 0]


@pytest.mark.parametrize(
    ("step", "every"),
    [(1, 1), (2, 1), (5, 1), (10, 2), (20, 5), (25, 5), (50, 10), (100, 25)],
)
def test_minor_spacing_by_step(step, every):
    from jaff.plotting._ticks import minor_every

    assert minor_every(step) == every


def test_seventy_decade_axis_sparse_minors():
    y = np.logspace(-80, -10, 50)
    fig, axes = _render(
        [{"series": [{"x": XS, "y": y}], "xscale": "log", "yscale": "log"}]
    )
    _draw(fig)
    ax = axes[0, 0]
    lo, hi = sorted(ax.get_ylim())
    majors = sorted(np.log10([v for v in ax.yaxis.get_majorticklocs() if lo <= v <= hi]))
    minors = np.log10([v for v in ax.yaxis.get_minorticklocs() if lo <= v <= hi])
    for a, b in zip(majors, majors[1:]):
        assert sum(a < m < b for m in minors) <= 4


def test_four_decade_axis_minor_every_decade():
    x = np.logspace(0, 4, 30)
    fig, axes = _render(
        [{"series": [{"x": x, "y": x}], "xscale": "log"}], style={"inside_ticks": False}
    )
    _draw(fig)
    ax = axes[0, 0]
    lo, hi = sorted(ax.get_xlim())
    locs = [*ax.xaxis.get_majorticklocs(), *ax.xaxis.get_minorticklocs()]
    decades = sorted({round(np.log10(v)) for v in locs if lo <= v <= hi})
    assert decades == [0, 1, 2, 3, 4]  # a tick (major or minor) on every decade


# --- shared legend needs one look per label ----------------------------------


def test_shared_legend_colour_consistent_across_panels():
    p = {"x": X, "y": Y, "label": "p"}
    q = {"x": X, "y": 2 * Y, "label": "q"}
    fig, axes = _render([{"series": [p, q]}, {"series": [q, p]}], grid=(1, 2))
    left, right = (
        [mcolors.to_hex(ln.get_color()) for ln in ax.lines] for ax in axes.ravel()
    )
    assert left == right[::-1] and left[0] != left[1]
    (legend,) = fig.legends
    assert [mcolors.to_hex(h.get_color()) for h in legend.legend_handles] == left


def test_shared_legend_distinct_colours_for_distinct_labels():
    fig, axes = _render(
        [
            {"series": [{"x": X, "y": Y, "label": f"a{k}"} for k in range(2)]},
            {"series": [{"x": X, "y": Y, "label": f"b{k}"} for k in range(2)]},
        ],
        grid=(1, 2),
    )
    colours = [mcolors.to_hex(ln.get_color()) for ax in axes.ravel() for ln in ax.lines]
    assert len(set(colours)) == 4


def test_touching_panels_trim_x_label_at_shared_spine():
    fig, axes = _render(
        [_panel()] * 2, grid=(1, 2), style={"base": "talk", "gap": 0}, sharey=True
    )
    r = _draw(fig)
    left, right = axes.ravel()
    spine = left.get_window_extent(r).x1
    labels = [t for t in left.get_xticklabels() if t.get_text()]
    assert all(t.get_window_extent(r).x1 <= spine + 0.5 for t in labels)


# --- review fixes (§14 round) -------------------------------------------------

XS = np.logspace(0, 2, 50)


def _titled(title, n=1, **kw):
    series = [{"x": XS, "y": XS * 1e-3 * (k + 1), "label": f"s{k}"} for k in range(n)]
    return {"xscale": "log", "yscale": "log", "title": title, "series": series, **kw}


def test_headroom_not_compounded_by_sharey():
    _, alone = _render([_titled("A")])
    _, shared = _render(
        [_titled("A"), _titled("B"), _titled("C")], grid=(1, 3), sharey=True
    )
    assert shared[0, 0].get_ylim()[1] == pytest.approx(
        alone[0, 0].get_ylim()[1], rel=1e-6
    )


def test_long_outside_titles_do_not_overlap():
    long = "A very very long panel title indeed"
    fig, axes = _render([_titled(long), _titled(long)], grid=(1, 2))
    r = _draw(fig)
    left, right = (ax.title.get_window_extent(r) for ax in axes.ravel())
    assert left.x1 <= right.x0 + 0.5
    assert left.x0 >= -0.5 and right.x1 <= fig.bbox.x1 + 0.5


def test_inside_labels_clear_of_inward_ticks():
    fig, axes = _render([_titled("A")])
    r = _draw(fig)
    ax = axes[0, 0]
    style = resolve_style("paper")
    tick_px = style.ticks["major_length"] * fig.dpi / 72
    x0 = ax.get_window_extent(r).x0
    assert all(t.get_window_extent(r).x0 >= x0 + tick_px for t in _ylabels(ax))


def test_unlabelled_series_never_share_a_legend_colour():
    lone = {"title": "single", "series": [{"x": XS, "y": XS}]}
    pair = {
        "series": [{"x": XS, "y": XS, "label": "A"}, {"x": XS, "y": 2 * XS, "label": "B"}]
    }
    fig, axes = _render([lone, pair], grid=(1, 2))
    (legend,) = fig.legends
    legend_colours = {mcolors.to_hex(h.get_color()) for h in legend.legend_handles}
    assert mcolors.to_hex(axes[0, 0].lines[0].get_color()) not in legend_colours


def test_figure_tall_enough_for_margin_legend():
    fig, _ = _render([_titled("A", n=12)])
    r = _draw(fig)
    (legend,) = fig.legends
    assert legend.get_window_extent(r).y0 >= -0.5
    assert fig.get_tightbbox(r).y0 >= -0.01


def test_legend_title_on_shared_legend():
    fig, _ = _render([_titled("A", n=2)], legend_title="Source")
    (legend,) = fig.legends
    title = legend.get_title()
    assert title.get_text() == "Source" and title.get_fontweight() in ("bold", 700)


def test_legend_title_heads_axes_groups():
    fig, axes = _render(
        [_titled("A", n=2), _titled("B", n=2)],
        grid=(1, 2),
        legend="axes",
        legend_title="Source",
    )
    r = _draw(fig)
    (head,) = [t for t in fig.texts if t.get_text() == "Source"]
    first = axes[0, 0].get_legend().get_window_extent(r)
    box = head.get_window_extent(r)
    assert box.y0 >= first.y1 - 0.5 and box.x0 >= axes[0, 1].get_window_extent(r).x1


def test_legend_title_validated():
    with pytest.raises(ValueError, match="legend_title"):
        FigureSpec.from_dict({"axes": [], "legend_title": 3})


def test_tag_sits_just_under_title_whatever_the_figsize():
    spec = {"title": "Title", "tag": "tag", "series": [{"x": XS, "y": XS, "label": "a"}]}
    fig, axes = _render([spec], figsize=(4, 6))
    r = _draw(fig)
    title = axes[0, 0].title.get_window_extent(r)
    (tag,) = axes[0, 0].texts
    gap_in = (title.y0 - tag.get_window_extent(r).y1) / fig.dpi
    assert 0 <= gap_in < 0.1


def test_explicit_xlim_puts_y_labels_outside():
    fig, axes = _render([_titled("A", xlim=(1, 100))])
    r = _draw(fig)
    ax = axes[0, 0]
    assert ax.get_xlim() == pytest.approx((1, 100))
    labels = _ylabels(ax)
    assert labels and all(
        t.get_window_extent(r).x1 <= ax.get_window_extent(r).x0 + 0.5 for t in labels
    )


def test_explicit_grid_over_cap_falls_back_to_groups(jaff_warnings):
    panels = [
        {"series": [{"x": XS, "y": XS * k, "label": f"p{i}-{k}"} for k in range(1, 8)]}
        for i in range(2)
    ]
    fig, axes = _render(panels, grid=(1, 2))
    assert not fig.legends
    assert all(len(ax.get_legend().get_texts()) == 7 for ax in axes.ravel())
    assert sum("legend_max" in m for m in jaff_warnings) == 1


@pytest.mark.parametrize(
    "kw",
    [
        {"style": {"panel_size": (1.9, 0.7)}},
        {"figsize": (2.2, 6.0)},
        {"style": {"panel_size": (1.9, 3.5)}},
    ],
)
def test_strip_matches_final_labels(kw):
    fig, axes = _render([_titled("A")], **kw)
    r = _draw(fig)
    ax = axes[0, 0]
    widest = max(t.get_window_extent(r).x1 for t in _ylabels(ax))
    assert ax.transData.transform((XS[0], 1e-3))[0] >= widest


def test_decorations_measured_at_final_geometry():
    xl = np.linspace(-5, 5, 50)
    panels = [
        {"title": "lin", "series": [{"x": xl, "y": -1000 * xl**2, "label": "q"}]},
        {"title": "xlim", "xlim": (-5, 5), "series": [{"x": xl, "y": xl, "label": "l"}]},
    ]
    fig, axes = _render(panels, grid=(1, 2))
    r = _draw(fig)
    left, right = (ax.get_window_extent(r) for ax in axes.ravel())
    lo, hi = axes[0, 1].get_ylim()
    shown = [
        t
        for t in axes[0, 1].get_yticklabels()
        if t.get_text() and lo <= t.get_position()[1] <= hi
    ]
    assert shown
    for t in shown:
        assert t.get_window_extent(r).x0 >= left.x1 - 0.5
    from jaff.plotting._layout import _decor

    final_left = _decor(axes[0, 1], r).left
    gap_in = (right.x0 - left.x1) / fig.dpi
    assert gap_in == pytest.approx(final_left + resolve_style("paper").gap, abs=0.01)


# --- §15: x-axis placement by grid row ---------------------------------------


def _xlabels(ax, top=False):
    ticks = ax.xaxis.get_major_ticks()
    lo, hi = sorted(ax.get_xlim())
    texts = [tk.label2 if top else tk.label1 for tk in ticks if lo <= tk.get_loc() <= hi]
    return [t for t in texts if t.get_visible() and t.get_text()]


def _col(n_rows, xs=None, **kw):
    xs = xs or [XS] * n_rows
    return [
        {
            "xscale": "log",
            "yscale": "log",
            "xlabel": "x",
            "series": [{"x": x, "y": x, "label": "a"}],
        }
        for x in xs
    ], kw


def test_column_top_middle_bottom_with_shared_x():
    panels, _ = _col(3)
    fig, axes = _render(panels, grid=(3, 1), sharex=True)
    _draw(fig)
    top, mid, bot = axes.ravel()
    assert _xlabels(top, top=True) and not _xlabels(top)
    assert not _xlabels(mid) and not _xlabels(mid, top=True)
    assert _xlabels(bot) and not _xlabels(bot, top=True)
    assert top.xaxis.get_label_position() == "top" and top.get_xlabel()
    assert bot.xaxis.get_label_position() == "bottom" and bot.get_xlabel()
    assert mid.get_xlabel() == ""


def test_middle_panel_with_other_range_gets_inside_labels():
    panels, _ = _col(3, xs=[XS, np.logspace(1, 4, 50), XS])
    fig, axes = _render(panels, grid=(3, 1))
    r = _draw(fig)
    mid = axes[1, 0]
    box = mid.get_window_extent(r)
    labels = _xlabels(mid)
    assert labels
    for t in labels:
        e = t.get_window_extent(r)
        assert e.y0 >= box.y0 and e.x0 >= box.x0 - 0.5 and e.x1 <= box.x1 + 0.5
    lowest = mid.transData.transform((1e1, 1e1))[1]
    assert lowest >= max(t.get_window_extent(r).y1 for t in labels)


def test_middle_inside_labels_clear_of_neighbours():
    panels, _ = _col(3, xs=[XS, np.logspace(1, 4, 50), XS])
    fig, axes = _render(panels, grid=(3, 1))
    _draw(fig)
    top, mid, bot = (a.get_window_extent() for a in axes.ravel())
    gap = resolve_style("paper").gap * fig.dpi
    assert mid.y0 - bot.y1 == pytest.approx(gap, abs=1)


def test_middle_other_range_with_ylim_keeps_labels_outside():
    panels, _ = _col(3, xs=[XS, np.logspace(1, 4, 50), XS])
    panels[1]["ylim"] = (1, 1e4)
    fig, axes = _render(panels, grid=(3, 1))
    r = _draw(fig)
    mid = axes[1, 0]
    labels = _xlabels(mid)
    assert labels and all(
        t.get_window_extent(r).y1 <= mid.get_window_extent(r).y0 + 0.5 for t in labels
    )
    assert mid.get_ylim() == pytest.approx((1, 1e4))


def test_single_panel_column_labels_at_bottom():
    fig, axes = _render([_panel(), _panel(), _panel()], grid=(2, 2), sharex=True)
    _draw(fig)
    right = axes[0, 1]
    assert _xlabels(right) and not _xlabels(right, top=True)
    assert _xlabels(axes[0, 0], top=True) and _xlabels(axes[1, 0])


def test_top_row_labels_inside_figure():
    fig, axes = _render([_panel(xlabel="N")] * 4, grid=(2, 2), sharex=True, title="T")
    r = _draw(fig)
    tight = fig.get_tightbbox(r)
    assert tight.y1 <= fig.get_size_inches()[1] + 0.01
    labels = _xlabels(axes[0, 0], top=True)
    assert (
        labels
        and max(t.get_window_extent(r).y1 for t in labels)
        < fig._suptitle.get_window_extent(r).y0
    )


def test_touching_top_labels_trimmed_at_spine():
    fig, axes = _render(
        [_panel()] * 4,
        grid=(2, 2),
        style={"base": "talk", "gap": 0},
        sharey=True,
        sharex=True,
    )
    r = _draw(fig)
    spine = axes[0, 0].get_window_extent(r).x1
    assert all(
        t.get_window_extent(r).x1 <= spine + 0.5 for t in _xlabels(axes[0, 0], top=True)
    )


# --- §15 review: outside titles vs top-spine labels, mixed roles ------------

LONG = "A very very long panel title indeed"


def _boxes_overlap(a, b, tol=0.5):
    return (
        min(a.x1, b.x1) - max(a.x0, b.x0) > tol
        and min(a.y1, b.y1) - max(a.y0, b.y0) > tol
    )


def _title_clear_of_top_labels(fig, axes):
    r = _draw(fig)
    tight = fig.get_tightbbox(r)
    w, h = fig.get_size_inches()
    assert tight.y1 <= h + 0.01 and tight.x0 >= -0.01 and tight.x1 <= w + 0.01
    for ax in axes.ravel():
        if not ax.axison or not ax.get_title():
            continue
        title = ax.title.get_window_extent(r)
        tops = [t.get_window_extent(r) for t in _xlabels(ax, top=True)]
        if ax.xaxis.get_label_position() == "top" and ax.get_xlabel():
            tops.append(ax.xaxis.label.get_window_extent(r))
        assert not any(_boxes_overlap(title, b) for b in tops)
    return r


def test_long_titles_clear_of_top_labels_2x2():
    panels = [_titled(LONG, xlabel="x") for _ in range(4)]
    fig, axes = _render(panels, grid=(2, 2), sharex=True)
    _title_clear_of_top_labels(fig, axes)


def test_long_titles_clear_of_top_labels_3x1_middle_inside():
    panels = [_titled(LONG, xlabel="x") for _ in range(3)]
    panels[1]["series"][0]["x"] = np.logspace(1, 4, 50)
    fig, axes = _render(panels, grid=(3, 1))
    _title_clear_of_top_labels(fig, axes)


def test_long_titles_side_by_side_do_not_overlap_2x2():
    panels = [_titled(LONG, xlabel="x") for _ in range(4)]
    fig, axes = _render(panels, grid=(2, 2), sharex=True)
    r = _draw(fig)
    for row in axes:
        left, right = (ax.title.get_window_extent(r) for ax in row)
        assert left.x1 <= right.x0 + 0.5


def test_mixed_roles_3x2_with_hidden_cell():
    panels, _ = _col(5)
    fig, axes = _render(panels, grid=(3, 2), sharex=True)
    _draw(fig)
    (a, b), (c, d), (e, hidden) = axes
    assert _xlabels(a, top=True) and not _xlabels(a)  # col 0: top
    assert not _xlabels(c) and not _xlabels(c, top=True)  # col 0: middle, same range
    assert _xlabels(e) and not _xlabels(e, top=True)  # col 0: bottom
    assert _xlabels(b, top=True) and not _xlabels(b)  # col 1: top
    assert _xlabels(d) and not _xlabels(d, top=True)  # col 1: bottom above hidden
    assert not hidden.axison


# --- §16: empty rows/columns collapse, solo colour, kappa label --------------


def test_empty_rows_collapse_keep_axes_shape():
    panels = [_titled(f"P{k}", xlabel="x", ylabel="y") for k in range(8)]
    fig, axes = _render(panels, grid=(4, 4))
    ref, _ = _render(panels, grid=(2, 4))
    assert axes.shape == (4, 4)
    assert tuple(fig.get_size_inches()) == pytest.approx(tuple(ref.get_size_inches()))
    assert not any(ax.get_visible() for ax in axes[2:].ravel())
    r = _draw(fig)
    span = [ax.get_window_extent(r) for ax in axes[:2].ravel()]
    mid_y = (min(b.y0 for b in span) + max(b.y1 for b in span)) / 2
    sup = fig._supylabel.get_window_extent(r)
    assert (sup.y0 + sup.y1) / 2 == pytest.approx(mid_y, abs=2)


def test_empty_columns_collapse():
    panels = [_titled(f"P{k}") for k in range(3)]
    fig, axes = _render(panels, grid=(1, 4))
    ref, _ = _render(panels, grid=(1, 3))
    assert axes.shape == (1, 4) and not axes[0, 3].get_visible()
    assert tuple(fig.get_size_inches()) == pytest.approx(tuple(ref.get_size_inches()))


def test_hidden_cell_in_used_row_still_hosts_legend():
    panels = [_titled(f"P{k}", n=2) for k in range(5)]
    fig, axes = _render(panels, grid=(3, 3))
    assert not any(ax.get_visible() for ax in axes[2])
    assert axes[1, 2].get_legend() is not None


def test_solo_series_uses_solo_color():
    fig, axes = _render([{"title": "A", "series": [{"x": XS, "y": XS}]}])
    assert mcolors.to_hex(axes[0, 0].lines[0].get_color()) == "#00798c"
    assert resolve_style("talk").solo_color == "#00798c"


def test_solo_color_avoids_legend_colours():
    lone = {"title": "single", "series": [{"x": XS, "y": XS}]}
    pair = {
        "series": [{"x": XS, "y": XS, "label": "A"}, {"x": XS, "y": 2 * XS, "label": "B"}]
    }
    # solo colour equal to legend entry A's colour -> next unused palette colour
    _, axes = _render([lone, pair], grid=(1, 2), style={"solo_color": "#0072b2"})
    colour = mcolors.to_hex(axes[0, 0].lines[0].get_color())
    assert colour not in ("#0072b2", "#d55e00")
    assert colour == themes.PALETTES["jaff"].cycle[2]


def test_solo_color_not_for_labelled_or_explicit():
    _, axes = _render([{"series": [{"x": XS, "y": XS, "color": "#123456"}]}])
    assert mcolors.to_hex(axes[0, 0].lines[0].get_color()) == "#123456"


@pytest.mark.parametrize("bad", ["notacolour", 3])
def test_solo_color_validated(bad):
    with pytest.raises(ValueError, match="style.solo_color"):
        resolve_style({"solo_color": bad})


def test_plot_rates_kappa_label():
    import sympy

    from jaff.plotting import plot_rates

    tgas = sympy.Symbol("tgas")
    from types import SimpleNamespace

    r = SimpleNamespace(
        rate=SimpleNamespace(coefficient=1e-10 * tgas),
        tmin=10,
        tmax=1000,
        get_latex=lambda: "R",
    )
    _, axes = plot_rates([r], show=False)
    assert axes[0, 0].get_ylabel() == r"$\log_{10}$ $\kappa$"


def test_inside_y_labels_keep_clear_of_inside_titles():
    # GOW: long reaction titles inside 2 x 4 panels next to inside y labels
    from jaff import Network
    from jaff.plotting import plot_rates
    from jaff.plotting.adapters import rate_series

    net = Network(str(REPO_ROOT / "networks/GOW/GOW.jet"))
    ok = [r for r in net.reactions if rate_series(r, quiet=True) is not None][:8]
    fig, axes = plot_rates(ok, grid=(4, 4), show=False, close=False)
    r = _draw(fig)
    clear = 3 * fig.dpi / 72
    for ax in axes[:2].ravel():
        title = ax.title.get_window_extent(r)
        for t in _ylabels(ax):
            e = t.get_window_extent(r)
            in_band = e.y1 > title.y0 and e.y0 < title.y1
            assert not in_band or e.x1 + clear <= title.x0


# --- §17: scale with column count ---------------------------------------------


def _row(n_cols, **kw):
    panels = [_titled(f"P{k}", xlabel="x", ylabel="y") for k in range(n_cols)]
    return _render(panels, grid=(1, n_cols), **kw)


def _tick_size(ax):
    return ax.yaxis.get_major_ticks()[0].label1.get_fontsize()


def test_six_columns_scale_one_and_a_half():
    f4, a4 = _row(4)
    f6, a6 = _row(6)
    assert _tick_size(a6[0, 0]) == pytest.approx(1.5 * _tick_size(a4[0, 0]))
    lw4, lw6 = a4[0, 0].lines[0].get_linewidth(), a6[0, 0].lines[0].get_linewidth()
    assert lw6 == pytest.approx(1.5 * lw4)
    # titles / tags do not scale (user decision on §17)
    assert a6[0, 0].title.get_fontsize() == pytest.approx(a4[0, 0].title.get_fontsize())
    _draw(f4)
    _draw(f6)
    w4 = a4[0, 0].get_window_extent().width / f4.dpi
    w6 = a6[0, 0].get_window_extent().width / f6.dpi
    assert w6 == pytest.approx(w4, rel=1e-3)


def test_four_or_fewer_columns_unscaled():
    _, a2 = _row(2)
    _, a4 = _row(4)
    assert _tick_size(a2[0, 0]) == pytest.approx(_tick_size(a4[0, 0]))
    assert _tick_size(a4[0, 0]) == pytest.approx(
        resolve_style("paper").rc["ytick.labelsize"]
    )


def test_explicit_figsize_unscaled():
    _, a6 = _row(6, figsize=(12, 2))
    assert _tick_size(a6[0, 0]) == pytest.approx(
        resolve_style("paper").rc["ytick.labelsize"]
    )


def test_registered_styles_untouched():
    import copy

    before = copy.deepcopy(themes.STYLES["paper"])
    _row(8)
    assert themes.STYLES["paper"] == before


def test_ref_cols_field():
    assert resolve_style("paper").ref_cols == 4 and resolve_style("talk").ref_cols == 3
    for bad in (0, 2.5, True):
        with pytest.raises(ValueError, match="style.ref_cols"):
            resolve_style({"ref_cols": bad})


def test_scale_is_not_a_user_field():
    with pytest.raises(ValueError, match="unknown field"):
        resolve_style({"scale": 2})


def test_scaled_style_helper():
    from jaff.plotting.themes import scaled

    paper = resolve_style("paper")
    big = scaled(paper, 2.0)
    assert big.rc["font.size"] == 2 * paper.rc["font.size"]
    assert big.lw == 2 * paper.lw and big.gap == 2 * paper.gap
    assert big.ticks["major_length"] == 2 * paper.ticks["major_length"]
    assert big.marker["gap_pt"] == 2 * paper.marker["gap_pt"]
    assert big.tag == paper.tag  # titles and tags keep their size
    assert big.panel_size == paper.panel_size and big.scale == 2.0
    assert scaled(paper, 1.0) is paper


def test_tag_offset_unscaled_with_columns():
    def gap_in(n):
        panel = {"title": "T", "tag": "t", "series": [{"x": XS, "y": XS, "label": "a"}]}
        fig, axes = _render([panel] * n, grid=(1, n))
        r = _draw(fig)
        title = axes[0, 0].title.get_window_extent(r)
        (tag,) = axes[0, 0].texts
        return (title.y0 - tag.get_window_extent(r).y1) / fig.dpi

    assert gap_in(8) == pytest.approx(gap_in(4), abs=0.005)
