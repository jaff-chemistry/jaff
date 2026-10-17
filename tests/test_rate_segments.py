# ABOUTME: Tests for RateSegments.evaluate_equivalent_rate multi-range assembly
# ABOUTME: Constant/temperature-dependent ranges, adjacent+separated, both clip modes

from pathlib import Path

import pytest
import numpy as np
from sympy import Float, Integer, lambdify, sqrt, symbols

from jaff import Network
from jaff.core.reaction.types import RateSegment, RateSegments

TGAS = symbols("tgas")


def _eval(segments, mode, t):
    """Assemble the piecewise rate and evaluate it at temperature *t*."""
    expr = RateSegments(list(segments), mode).sort().evaluate_equivalent_rate("r")
    return expr.subs(TGAS, t) if hasattr(expr, "subs") else expr


# --------------------------------------------------------------------------- #
# two constant ranges                                                         #
# --------------------------------------------------------------------------- #
def test_two_constant_adjacent_clip():
    segs = [RateSegment(Integer(1), 10, 100), RateSegment(Integer(2), 100, 200)]
    assert _eval(segs, "clip", 50) == 1
    assert _eval(segs, "clip", 150) == 2  # bug: returned 1
    assert _eval(segs, "clip", 5) == 1    # clipped low
    assert _eval(segs, "clip", 300) == 2  # clipped high


def test_two_constant_adjacent_extrapolate():
    segs = [RateSegment(Integer(1), 10, 100), RateSegment(Integer(2), 100, 200)]
    assert _eval(segs, "extrapolate", 50) == 1
    assert _eval(segs, "extrapolate", 150) == 2
    assert _eval(segs, "extrapolate", 5) == 1
    assert _eval(segs, "extrapolate", 300) == 2


def test_two_constant_separated_clip_interpolates_gap():
    segs = [RateSegment(Integer(1), 10, 100), RateSegment(Integer(2), 150, 250)]
    assert _eval(segs, "clip", 50) == 1
    assert _eval(segs, "clip", 200) == 2
    # gap [100, 150]: geometric interp, midpoint 125 -> sqrt(1*2)
    assert _eval(segs, "clip", 125) == sqrt(2)


def test_two_constant_separated_extrapolate():
    segs = [RateSegment(Integer(1), 10, 100), RateSegment(Integer(2), 150, 250)]
    assert _eval(segs, "extrapolate", 50) == 1
    assert _eval(segs, "extrapolate", 200) == 2
    assert _eval(segs, "extrapolate", 125) == sqrt(2)


# --------------------------------------------------------------------------- #
# constant followed by temperature-dependent range                           #
# --------------------------------------------------------------------------- #
def test_constant_then_tempdep_adjacent_clip():
    segs = [RateSegment(Integer(1), 10, 100), RateSegment(TGAS, 100, 200)]
    assert _eval(segs, "clip", 50) == 1
    assert _eval(segs, "clip", 150) == 150
    assert _eval(segs, "clip", 5) == 1     # clipped low -> constant
    assert _eval(segs, "clip", 300) == 200  # clipped high -> tmax value


def test_constant_then_tempdep_separated_extrapolate():
    segs = [RateSegment(Integer(1), 10, 100), RateSegment(TGAS, 150, 250)]
    assert _eval(segs, "extrapolate", 50) == 1
    assert _eval(segs, "extrapolate", 200) == 200
    # gap [100,150]: interp between the two BOUNDING rates -- const 1 at a=100
    # and tgas evaluated at b=150 (=150); geometric midpoint at 125 -> sqrt(150)
    assert _eval(segs, "extrapolate", 125) == sqrt(150)


# --------------------------------------------------------------------------- #
# separated temperature-dependent ranges: gap interpolates BOUNDING rates      #
# --------------------------------------------------------------------------- #
def test_two_tempdep_separated_interp_uses_boundary_rates():
    """Gap bridges the two boundary constants, not the extrapolated fits.

    seg1 rate=tgas valid [10,100] -> boundary rate at a=100 is 100.
    seg2 rate=2*tgas valid [200,300] -> boundary rate at b=200 is 400.
    Gap [100,200] is a straight line in log k from 100 to 400; midpoint 150
    -> 200.  (The old code blended both fits extrapolated into the gap.)
    """
    segs = [RateSegment(TGAS, 10, 100), RateSegment(2 * TGAS, 200, 300)]
    assert _eval(segs, "clip", 100) == 100   # left boundary rate
    assert _eval(segs, "clip", 200) == 400   # right boundary rate
    assert _eval(segs, "clip", 150) == 200   # geometric midpoint of gap
    # interior fraction: 100 * 4**((125-100)/(200-100)) = 100*sqrt(2)
    assert _eval(segs, "clip", 125) == 100 * sqrt(2)


def test_gap_nonpositive_boundary_falls_back_to_linear():
    """Log-space bridging needs positive rates; otherwise interpolate linearly."""
    segs = [RateSegment(Integer(-1), 10, 100), RateSegment(Integer(3), 150, 250)]
    assert _eval(segs, "clip", 125) == 1


def test_steep_gap_is_smooth_in_log_and_finite_everywhere():
    """KIDA S + C+ -> C + S+: a ~600x drop across the 1 K gap 499-500 K.

    A linear bridge falls to ~0 just below 500 K (a log singularity that rate
    table refinement can never resolve).  The log-space bridge must be log-
    linear and, since numpy/C evaluate every Piecewise branch, finite when
    evaluated far outside the gap.
    """
    hi, lo = Float(1.3e-9), Float(2.2e-12)
    segs = [RateSegment(hi, 10, 499), RateSegment(lo, 500, 50000)]
    expr = RateSegments(segs, "clip").sort().evaluate_equivalent_rate("r")
    gap_expr = expr.args[1].expr  # the bridge branch
    f = lambdify(TGAS, gap_expr, "numpy")

    t = np.linspace(499.0, 500.0, 11)
    np.testing.assert_allclose(
        np.log(f(t)), np.log(1.3e-9) + (t - 499.0) * np.log(2.2e-12 / 1.3e-9)
    )
    with np.errstate(over="ignore", under="ignore"):
        far = f(np.array([5.0, 5000.0]))
    assert not np.any(np.isnan(far))


# --------------------------------------------------------------------------- #
# Network-level regression (the reported repro)                               #
# --------------------------------------------------------------------------- #
def test_two_touching_tempdep_clip():
    """Two touching temperature-dependent ranges assemble each body (issue 045)."""
    segs = [RateSegment(TGAS, 10, 100), RateSegment(2 * TGAS, 100, 200)]
    assert _eval(segs, "clip", 50) == 50
    assert _eval(segs, "clip", 150) == 300   # bug returned 400 (upper clip)
    assert _eval(segs, "clip", 100) == 200   # shared boundary -> upper range
    assert _eval(segs, "clip", 5) == 10      # clipped low -> first tmin value
    assert _eval(segs, "clip", 300) == 400   # clipped high -> last tmax value


def test_two_touching_tempdep_extrapolate():
    segs = [RateSegment(TGAS, 10, 100), RateSegment(2 * TGAS, 100, 200)]
    assert _eval(segs, "extrapolate", 50) == 50
    assert _eval(segs, "extrapolate", 150) == 300
    assert _eval(segs, "extrapolate", 100) == 200
    assert _eval(segs, "extrapolate", 5) == 5      # first body extends down
    assert _eval(segs, "extrapolate", 300) == 600  # last body extends up


def test_three_touching_tempdep_clip():
    segs = [
        RateSegment(TGAS, 10, 100),
        RateSegment(2 * TGAS, 100, 200),
        RateSegment(3 * TGAS, 200, 300),
    ]
    # inside every segment
    assert _eval(segs, "clip", 50) == 50
    assert _eval(segs, "clip", 150) == 300
    assert _eval(segs, "clip", 250) == 750
    # shared boundaries -> upper range governs (strict tgas < tmax)
    assert _eval(segs, "clip", 100) == 200
    assert _eval(segs, "clip", 200) == 600
    # clipped ends
    assert _eval(segs, "clip", 5) == 10
    assert _eval(segs, "clip", 350) == 900


def test_three_touching_tempdep_extrapolate():
    segs = [
        RateSegment(TGAS, 10, 100),
        RateSegment(2 * TGAS, 100, 200),
        RateSegment(3 * TGAS, 200, 300),
    ]
    assert _eval(segs, "extrapolate", 50) == 50
    assert _eval(segs, "extrapolate", 150) == 300
    assert _eval(segs, "extrapolate", 250) == 750
    assert _eval(segs, "extrapolate", 100) == 200
    assert _eval(segs, "extrapolate", 200) == 600
    assert _eval(segs, "extrapolate", 5) == 5       # first extends down
    assert _eval(segs, "extrapolate", 350) == 1050  # last extends up


def test_network_two_touching_tempdep_ranges(tmp_path):
    dat = tmp_path / "tempdep.dat"
    dat.write_text(
        "@format:idx,R,R,P,tmin,tmax,rate\n"
        "1,H,H,H2,10,100,tgas\n"
        "2,H,H,H2,100,200,2*tgas\n"
    )
    net = Network(str(dat), funcfile=False)
    assert net.reactions[0].rate.subs(TGAS, 150) == 300


def test_network_two_constant_ranges(tmp_path):
    dat = tmp_path / "constant.dat"
    dat.write_text(
        "@format:idx,R,R,P,tmin,tmax,rate\n"
        "1,H,H,H2,10,100,1\n"
        "2,H,H,H2,100,200,2\n"
    )
    net = Network(str(dat), funcfile=False)
    rate = net.reactions[0].rate
    assert rate.subs(TGAS, 150) == 2
    assert rate.subs(TGAS, 50) == 1


# --------------------------------------------------------------------------- #
# open lower bound (tmin None) sorts first                                    #
# --------------------------------------------------------------------------- #
def test_sort_open_lower_bound_first():
    """A ``tmin=None`` segment sorts below any finite bound (was a TypeError)."""
    high = RateSegment(Integer(2), 100, None)
    low = RateSegment(Integer(1), None, 100)
    segs = RateSegments([high, low], "clip").sort()
    assert list(segs) == [low, high]


@pytest.mark.parametrize(
    "bounds",
    [
        [(None, 100), (None, 200)],  # two open lower bounds
        [(100, None), (200, None)],  # two open upper bounds
        [(10, None), (300, 500)],  # open upper bound below another range
    ],
)
def test_sort_rejects_open_bound_off_the_edge(bounds):
    """Only the lowest range may lack ``tmin``; only the highest may lack ``tmax``."""
    segs = [RateSegment(Integer(i), lo, hi) for i, (lo, hi) in enumerate(bounds)]
    with pytest.raises(ValueError, match="open"):
        RateSegments(segs, "clip").sort()


def test_sort_allows_lone_fully_open_segment():
    seg = RateSegment(Integer(1), None, None)
    assert list(RateSegments([seg], "clip").sort()) == [seg]


def test_network_open_ended_ranges(tmp_path):
    dat = tmp_path / "open.dat"
    dat.write_text(
        "@format:idx,R,R,P,tmin,tmax,rate\n"
        "1,H,H,H2,100,NONE,2\n"
        "2,H,H,H2,NONE,100,1\n"
    )
    net = Network(str(dat), funcfile=False)
    rate = net.reactions[0].rate
    assert rate.subs(TGAS, 50) == 1
    assert rate.subs(TGAS, 150) == 2
