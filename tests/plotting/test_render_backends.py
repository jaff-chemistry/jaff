# ABOUTME: jaff.plotting.render works on every backend, not only Agg: text measuring
# ABOUTME: and layout use the figure's renderer, not the Agg-only canvas.get_renderer

import matplotlib.pyplot as plt
import numpy as np
import pytest

from jaff.plotting import render
from jaff.plotting.spec import FigureSpec

X = np.logspace(0, 2, 20)


@pytest.fixture(params=["pdf", "svg", "ps"])
def backend(request):
    before = plt.get_backend()
    plt.switch_backend(request.param)
    yield request.param
    plt.close("all")
    plt.switch_backend(before)


def test_render_on_a_non_agg_backend(backend, tmp_path):
    spec = {
        "axes": [
            {
                "series": [{"x": X, "y": X, "label": "a"}, {"x": X, "y": 2 * X}],
                "title": "CO",
                "xscale": "log",
                "yscale": "log",
            }
        ]
        * 2,
        "grid": (1, 2),
    }

    fig, _ = render(FigureSpec.from_dict(spec))
    fig.savefig(tmp_path / f"out.{backend}")

    assert (tmp_path / f"out.{backend}").stat().st_size > 0
