# ABOUTME: Every bundled network under networks/ loads into a full Network object
# ABOUTME: Covers parsing plus Network post-processing (rate segments, species checks)

from pathlib import Path
from typing import Any, Dict, List

import pytest

from jaff import Network
from jaff.physics import DustProps, RadiationProps

NETWORKS_DIR = Path(__file__).resolve().parent.parent / "networks"

NETWORK_FILES: List[str] = sorted(
    p.relative_to(NETWORKS_DIR).as_posix() for p in NETWORKS_DIR.glob("**/*.jet")
)

# Networks that need constructor arguments to load (see the docs for each).
LOAD_KWARGS: Dict[str, Dict[str, Any]] = {
    "GOW++/GOW++.jet": {
        "radiation_props": RadiationProps(bands=[6.0, 11.2, 13.6]),
        "dust_props": DustProps(rv=3.1),
    },
}

# Networks whose source data is currently broken. strict=True turns a fix into
# a failure, so the marker is removed once the network loads again.
KNOWN_FAILURES: Dict[str, pytest.MarkDecorator] = {
    "rate22_final/rate22_final.rates.jet": pytest.mark.xfail(
        strict=True,
        raises=ValueError,
        reason="H2CCO+ + e- -> CO + CH2 duplicated with overlapping temperature "
        "ranges (lines 1452, 2061)",
    ),
}


def _params() -> List[Any]:
    return [
        pytest.param(name, marks=KNOWN_FAILURES.get(name, ()), id=name)
        for name in NETWORK_FILES
    ]


@pytest.mark.parametrize("name", _params())
def test_network_loads(name: str):
    net = Network(str(NETWORKS_DIR / name), **LOAD_KWARGS.get(name, {}))

    assert len(net.reactions) > 0


def test_special_cases_name_existing_networks():
    """A renamed network must not silently lose its kwargs or xfail marker."""
    stale = (set(LOAD_KWARGS) | set(KNOWN_FAILURES)) - set(NETWORK_FILES)

    assert not stale
