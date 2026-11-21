# ABOUTME: Tests for log y ranges in jaff.plotting.render: data spanning hundreds of
# ABOUTME: decades (subnormal rates) is capped to LOG_MAX_DECADES below its top, never inf

import importlib

import numpy as np
import pytest

from jaff.plotting.render import render
from jaff.plotting.spec import FigureSpec

_render_module = importlib.import_module("jaff.plotting.render")
LOG_MAX_DECADES = _render_module.LOG_MAX_DECADES

X = np.linspace(1.0, 50.0, 50)
# 1e-8 down to 1e-320 (subnormal): what an exp(-E/T) rate does at low T.
DEEP = 10.0 ** np.linspace(-8.0, -320.0, 50)


def _ylim(y, **axes):
    fig = render(
        FigureSpec.from_dict(
            {"axes": [{"series": [{"x": X, "y": y}], "yscale": "log", **axes}]}
        )
    )
    lo, hi = fig[1][0, 0].get_ylim()

    return float(lo), float(hi)


@pytest.mark.parametrize(
    ("y", "span"), [((3e-320, 3e-8), LOG_MAX_DECADES), ((1e-20, 1e-8), 12.0)]
)
def test_log_span_is_finite_and_capped(y, span):
    import matplotlib.pyplot as plt

    _, ax = plt.subplots()
    ax.set_yscale("log")
    ax.plot([1.0, 2.0], y)

    assert _render_module._span(ax, True) == pytest.approx(span)


@pytest.mark.filterwarnings("error")
@pytest.mark.parametrize("axes", [{}, {"title": "CO"}], ids=["plain", "inside_text"])
def test_deep_data_is_capped_below_its_top(axes):
    lo, hi = _ylim(DEEP, **axes)

    assert np.isfinite([lo, hi]).all()
    assert hi >= 1e-8
    assert lo == pytest.approx(1e-8 * 10.0**-LOG_MAX_DECADES, rel=1e-6)


@pytest.mark.filterwarnings("error")
def test_explicit_ylim_wins():
    assert _ylim(DEEP, ylim=(1e-300, 1e-5)) == pytest.approx((1e-300, 1e-5), rel=1e-9)


def test_data_within_the_cap_is_unchanged():
    y = 10.0 ** np.linspace(-8.0, -20.0, 50)
    lo, _ = _ylim(y)

    # The usual autoscale margin below the data, no cap.
    assert 1e-22 < lo < 1e-20
