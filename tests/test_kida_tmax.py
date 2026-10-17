# ABOUTME: KIDA Tmax handling: 9999 is "no upper limit" unless another range of the
# ABOUTME: same reaction continues above it; every other Tmax is a real bound

from pathlib import Path

from jaff import Network
from jaff.core.parsers.network._formats.kida.parser import _resolve_open_tmax

REPO = Path(__file__).resolve().parent.parent
KIDA_NET = REPO / "networks" / "kida_uva_2024" / "gas_reactions_kida.uva.2024.jet"


def _row(tmin, tmax, r=("H", "Cl"), p=("HCl",)):
    return {"r": list(r), "p": list(p), "tmin": tmin, "tmax": tmax}


def test_lone_9999_is_open():
    rows = [_row(None, 9999.0)]
    _resolve_open_tmax(rows)
    assert rows[0]["tmax"] is None


def test_9999_followed_by_higher_range_is_real():
    rows = [_row(3000.0, 9999.0), _row(10000.0, 20000.0)]
    _resolve_open_tmax(rows)
    assert [r["tmax"] for r in rows] == [9999.0, 20000.0]


def test_9999_duplicate_ranges_stay_open():
    """Exact duplicates (resolved later by duplicate_policy) do not continue."""
    rows = [_row(None, 9999.0), _row(None, 9999.0)]
    _resolve_open_tmax(rows)
    assert [r["tmax"] for r in rows] == [None, None]


def test_other_reaction_does_not_continue_range():
    rows = [_row(3000.0, 9999.0), _row(10000.0, 20000.0, r=("H", "O"), p=("OH",))]
    _resolve_open_tmax(rows)
    assert rows[0]["tmax"] is None


def test_bounds_above_9999_are_real(tmp_path):
    """Rows bounded at 15000/30000 used to be read as open-ended and overlap."""
    lines = [
        line
        for line in KIDA_NET.read_text().splitlines()
        if line.startswith("H          O                      OH         Photon")
    ]
    assert len(lines) == 7
    dat = tmp_path / "h_o.dat"
    dat.write_text("\n".join(lines) + "\n")

    rxn = Network(str(dat), funcfile=False).reactions[0]
    bounds = [(s.tmin, s.tmax) for s in rxn.rate_segments]
    assert bounds[-2:] == [(5001.0, 15000.0), (15001.0, 30000.0)]
