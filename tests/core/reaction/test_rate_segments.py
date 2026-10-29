# ABOUTME: Tests for RateSegments.evaluate_equivalent_rate multi-range assembly
# ABOUTME: Constant/temperature-dependent ranges, adjacent+separated, both clip modes


from types import SimpleNamespace

import pytest
from sympy import Integer, Rational, symbols

from jaff import Network
from jaff.core.reaction.types import RateSegment, RateSegments, ReactionSource

TGAS = symbols("tgas")


def _owner(mode):
    """A stand-in reaction whose ``t_cutoff`` setting is *mode*."""
    return SimpleNamespace(
        config=SimpleNamespace(t_cutoff=mode), verbatim="H + H -> H2", type="unknown"
    )


def _eval(segments, mode, t):
    """Assemble the piecewise rate and evaluate it at temperature *t*."""
    expr = RateSegments(_owner(mode), list(segments)).sort().evaluate_equivalent_rate()
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
    # gap [100, 150]: linear interp, midpoint 125 -> 1.5
    assert _eval(segs, "clip", 125) == Rational(3, 2)


def test_two_constant_separated_extrapolate():
    segs = [RateSegment(Integer(1), 10, 100), RateSegment(Integer(2), 150, 250)]
    assert _eval(segs, "extrapolate", 50) == 1
    assert _eval(segs, "extrapolate", 200) == 2
    assert _eval(segs, "extrapolate", 125) == Rational(3, 2)


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
    # and tgas evaluated at b=150 (=150); at 125 ->
    # (1*(150-125) + 150*(125-100)) / (150-100) = (25 + 3750)/50 = 75.5
    assert _eval(segs, "extrapolate", 125) == Rational(25 + 150 * 25, 50)


# --------------------------------------------------------------------------- #
# separated temperature-dependent ranges: gap interpolates BOUNDING rates      #
# --------------------------------------------------------------------------- #
def test_two_tempdep_separated_interp_uses_boundary_rates():
    """Gap bridges the two boundary constants, not the extrapolated fits.

    seg1 rate=tgas valid [10,100] -> boundary rate at a=100 is 100.
    seg2 rate=2*tgas valid [200,300] -> boundary rate at b=200 is 400.
    Gap [100,200] is a straight line 100->400; midpoint 150 -> 250.
    (The old code blended both fits extrapolated into the gap -> 225.)
    """
    segs = [RateSegment(TGAS, 10, 100), RateSegment(2 * TGAS, 200, 300)]
    assert _eval(segs, "clip", 100) == 100   # left boundary rate
    assert _eval(segs, "clip", 200) == 400   # right boundary rate
    assert _eval(segs, "clip", 150) == 250   # midpoint of straight gap
    # interior fraction: 100 + (400-100)*(125-100)/(200-100) = 175
    assert _eval(segs, "clip", 125) == 175


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
    assert net.reactions[0].rate_coefficient.subs(TGAS, 150) == 300


def test_network_two_constant_ranges(tmp_path):
    dat = tmp_path / "constant.dat"
    dat.write_text(
        "@format:idx,R,R,P,tmin,tmax,rate\n"
        "1,H,H,H2,10,100,1\n"
        "2,H,H,H2,100,200,2\n"
    )
    net = Network(str(dat), funcfile=False)
    rate = net.reactions[0].rate_coefficient
    assert rate.subs(TGAS, 150) == 2
    assert rate.subs(TGAS, 50) == 1


# --------------------------------------------------------------------------- #
# open lower bound (tmin None) sorts first                                    #
# --------------------------------------------------------------------------- #
def test_sort_open_lower_bound_first():
    """A ``tmin=None`` segment sorts below any finite bound (was a TypeError)."""
    high = RateSegment(Integer(2), 100, None)
    low = RateSegment(Integer(1), None, 100)
    segs = RateSegments(_owner("clip"), [high, low]).sort()
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
        RateSegments(_owner("clip"), segs).sort()


def test_sort_allows_lone_fully_open_segment():
    seg = RateSegment(Integer(1), None, None)
    assert list(RateSegments(_owner("clip"), [seg]).sort()) == [seg]


def test_network_open_ended_ranges(tmp_path):
    dat = tmp_path / "open.dat"
    dat.write_text(
        "@format:idx,R,R,P,tmin,tmax,rate\n"
        "1,H,H,H2,100,NONE,2\n"
        "2,H,H,H2,NONE,100,1\n"
    )
    net = Network(str(dat), funcfile=False)
    rate = net.reactions[0].rate_coefficient
    assert rate.subs(TGAS, 50) == 1
    assert rate.subs(TGAS, 150) == 2


def test_segment_is_immutable():
    from dataclasses import FrozenInstanceError

    with pytest.raises(FrozenInstanceError):
        RateSegment(Integer(1), 10, 100).rate = Integer(2)


def test_map_rates_replaces_every_rate():
    segs = RateSegments(
        _owner("clip"), [RateSegment(Integer(1), 10, 100), RateSegment(Integer(2), 100, 200)]
    )
    segs.map_rates(lambda rate: 3 * rate)

    assert [s.rate for s in segs] == [3, 6]
    assert segs._by_prop[(100, 200)].rate == 6


def test_equivalent_rate_tracks_changes():
    segs = RateSegments(_owner("clip"), [RateSegment(Integer(1), 10, 100)])
    assert segs.equivalent_rate.subs(TGAS, 50) == 1

    segs.add(RateSegment(Integer(2), 100, 200))
    assert segs.equivalent_rate.subs(TGAS, 150) == 2

    segs.map_rates(lambda rate: 10 * rate)
    assert segs.equivalent_rate.subs(TGAS, 150) == 20


def test_reaction_temperature_range_is_derived(tmp_path):
    dat = tmp_path / "r.dat"
    dat.write_text(
        "@format:idx,R,R,P,tmin,tmax,rate\n1,H,H,H2,10,100,1.0\n2,H,H,H2,100,200,2.0\n"
    )
    rxn = Network(str(dat), funcfile=False).reactions[0]

    assert (rxn.tmin, rxn.tmax) == (10, 200)
    with pytest.raises(AttributeError):
        rxn.tmin = 5


# --------------------------------------------------------------------------- #
# range errors name the reaction and each segment's source line               #
# --------------------------------------------------------------------------- #
def _sourced(bounds):
    """Segments over *bounds*, the i-th read from file line ``10 + i``."""
    return [
        RateSegment(Integer(i), lo, hi, ReactionSource(f"rate {i}", 10 + i, i, "krome"))
        for i, (lo, hi) in enumerate(bounds)
    ]


def _error(bounds, mode="clip"):
    segs = RateSegments(_owner(mode), _sourced(bounds))
    with pytest.raises(ValueError) as info:
        segs.sort().evaluate_equivalent_rate()

    return str(info.value)


@pytest.mark.parametrize(
    "bounds",
    [
        [(None, 400), (None, None), (400, None)],  # open below, not lowest
        [(10, None), (300, 500)],  # open above, not highest
        [(10, 200), (100, 300)],  # overlapping ranges
    ],
)
def test_range_error_names_reaction(bounds):
    assert "H + H -> H2 [type=unknown]" in _error(bounds)


@pytest.mark.parametrize(
    "bounds",
    [
        [(None, 400), (None, None), (400, None)],
        [(10, None), (300, 500)],
        [(10, 200), (100, 300)],
    ],
)
def test_range_error_lists_every_source_line(bounds):
    msg = _error(bounds)

    assert all(f"line {10 + i}: rate {i}" in msg for i in range(len(bounds)))


def test_range_error_lists_ranges_in_sorted_order():
    msg = _error([(400, None), (None, 400), (None, None)])
    rows = [ln.strip() for ln in msg.splitlines()[1:]]

    assert [r.split(")")[0] + ")" for r in rows] == [
        "(None, 400)",
        "(None, None)",
        "(400, None)",
    ]


def test_range_error_without_source():
    segs = RateSegments(
        _owner("clip"),
        [RateSegment(Integer(1), 10, 200), RateSegment(Integer(2), 100, 300)],
    )
    with pytest.raises(ValueError, match="no source line"):
        segs.sort().evaluate_equivalent_rate()


def test_network_range_error_names_file_lines(tmp_path):
    dat = tmp_path / "open.dat"
    dat.write_text(
        "@format:idx,R,R,P,tmin,tmax,rate\n"
        "1,H,H,H2,NONE,400,1\n"
        "2,H,H,H2,400,NONE,2\n"
        "@format:idx,R,R,P,rate\n"
        "3,H,H,H2,3\n"
    )
    with pytest.raises(ValueError) as info:
        Network(str(dat), funcfile=False)

    msg = str(info.value)
    assert "H + H -> H2" in msg
    assert "line 2: 1,H,H,H2,NONE,400,1" in msg
    assert "line 3: 2,H,H,H2,400,NONE,2" in msg
    assert "line 5: 3,H,H,H2,3" in msg
