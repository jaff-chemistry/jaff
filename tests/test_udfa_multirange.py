# ABOUTME: Regression tests for UDFA/RATE22 multi-temperature-range parsing.
# ABOUTME: One reaction line with N ranges must yield N segment records.

from pathlib import Path

import pytest

from jaff.core.parsers.network._formats.udfa.reaction import UdfaReaction
from jaff.errors import ParserError

REPO = Path(__file__).resolve().parent.parent
BUNDLED = REPO / "networks" / "rate22_final" / "rate22_final.rates.jet"


def _bundled_line(index: str) -> str:
    """Return the raw bundled RATE22 line whose leading field is *index*."""
    with open(BUNDLED, encoding="utf-8") as f:
        for line in f:
            if line.strip().split(":", 1)[0] == index:
                return line
    raise AssertionError(f"reaction {index} not found in {BUNDLED}")


def _segments(line: str) -> list[dict]:
    """Parse a UDFA line into its list of per-range segment dicts."""
    result = UdfaReaction().parse(line, 1, {}, BUNDLED)
    assert isinstance(result, list), "parse() must return a list of segments"
    return result


def test_rxn74_yields_both_fits() -> None:
    """Bundled reaction 74 has two fits; both must survive parsing (pre-merge)."""
    segs = _segments(_bundled_line("74"))
    assert len(segs) == 2

    bounds = {(s["tmin"], s["tmax"]) for s in segs}
    assert bounds == {(10.0, 100.0), (101.0, 3000.0)}

    by_bounds = {(s["tmin"], s["tmax"]): s for s in segs}
    # First fit: alpha=4.82e-9, beta=0.02, gamma=4.3
    assert "4.82e-09" in by_bounds[(10.0, 100.0)]["rate"]
    # Second fit: alpha=4.32e-9, beta=-0.39, gamma=39.4 — previously dropped
    assert "4.32e-09" in by_bounds[(101.0, 3000.0)]["rate"]


def test_rxn74_shared_participants() -> None:
    """Every segment of a multi-range line shares reactants/products/type."""
    segs = _segments(_bundled_line("74"))
    first = segs[0]
    for s in segs[1:]:
        assert s["r"] == first["r"]
        assert s["p"] == first["p"]
        assert s["type"] == first["type"]
    assert first["r"] == ["H-", "H"]
    assert first["p"] == ["H2", "e-"]


@pytest.mark.parametrize("index,count", [("1", 1), ("74", 2), ("3177", 3), ("7298", 4)])
def test_segment_count_matches_range_flag(index: str, count: int) -> None:
    """The number of emitted segments equals the declared range-count flag."""
    assert len(_segments(_bundled_line(index))) == count


def test_single_range_unchanged() -> None:
    """A single-range line still yields one segment with its full bounds."""
    (seg,) = _segments(_bundled_line("1"))
    assert (seg["tmin"], seg["tmax"]) == (10.0, None)  # tmax 41000 clamps to None
    assert "5.00e-10" in seg["rate"]


def test_quoted_note_with_comma_preserved() -> None:
    """Metadata quotes containing commas must not corrupt range tokenisation."""
    # Reaction 7298 note: "West N et al, PCCP, 25, 7719 (2023)"
    segs = _segments(_bundled_line("7298"))
    assert len(segs) == 4
    bounds = [(s["tmin"], s["tmax"]) for s in segs]
    assert (4.0, 20.0) in bounds
    assert (20.1, 100.0) in bounds


def test_range_count_too_high_rejected() -> None:
    """A flag that promises more blocks than present is rejected, not truncated."""
    # flag says 2 ranges but only one parameter block follows.
    line = '1:AD:H-:H:H2:e-:::2:4.82e-09:0.02:4.3:10:100:M:A:"doi":"note":\n'
    with pytest.raises(ParserError):
        UdfaReaction().parse(line, 1, {}, BUNDLED)
