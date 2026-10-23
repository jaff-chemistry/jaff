# ABOUTME: Characterization snapshot of NetworkParser.get_parsed() over all networks
# ABOUTME: Regression net guarding the format-ownership refactor (behavior-preserving)

from pathlib import Path

import pytest

from jaff.core.parsers.network._engine import NetworkParser

# Every bundled network file the parser can read. `*.rates.jet` is a subset of
NETWORK_FILES = sorted(
    {
        p.as_posix()
        for pat in ("*.jet", "*.dat")
        for p in Path("networks").glob(f"**/{pat}")
    }
)


def _fingerprint(path):
    with NetworkParser(path) as p:
        parsed, _ = p.get_parsed()
    return [
        (i, tuple(r["r"]), tuple(r["p"]), r["tmin"], r["tmax"], r["rate"], r["type"])
        for i, r in enumerate(parsed)
    ]


@pytest.mark.parametrize("path", NETWORK_FILES)
def test_parser_output_is_stable(path, snapshot_store):
    fp = _fingerprint(path)
    snapshot_store.check(path, fp)
