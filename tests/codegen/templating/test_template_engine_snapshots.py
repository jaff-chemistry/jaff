# ABOUTME: Snapshots of the $$JAFF template engine: render fixture templates against
# ABOUTME: small networks and compare the whole output, echo included, with committed renders

# Regression net for the template engine.  Each template in
# tests/fixtures/templates/engine/ is rendered for the networks listed in CASES and
# compared line for line with the committed render in
# tests/snapshots/template_engine/<network>/.  The snapshots were first produced by
# the old engine; the rewrite reproduced them byte for byte.
#
# They pin current behaviour, including the remaining quirks listed in
# parser_bugs.md.
# After an intentional change refresh them with
#   JAFF_UPDATE_SNAPSHOTS=1 uv run pytest tests/test_template_engine_snapshots.py
# and review the diff.

import os
import re
from functools import cache
from pathlib import Path
from typing import Callable, Dict, List, Tuple

import pytest

from jaff import Network
from jaff.codegen import TemplateParser
from jaff.physics.photo_reactions.config import RadiationConfig
from tests._paths import FIXTURES_DIR, REPO_ROOT, SNAPSHOTS_DIR

REPO = REPO_ROOT
FIXTURES = FIXTURES_DIR
TEMPLATES = FIXTURES / "templates" / "engine"
SNAPSHOTS = SNAPSHOTS_DIR / "template_engine"
UPDATE = bool(os.environ.get("JAFF_UPDATE_SNAPSHOTS"))

# Numeric literals (incl. Fortran's 1.0d0). Floats from numerical integration (the
# photoionisation rates) differ in the last digit between platforms, so numbers are
# compared with pytest.approx and the text around them exactly.
NUMBER = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?")
REL = 1e-12


def _numbers(text: str) -> List[float]:
    return [float(n.replace("d", "e").replace("D", "e")) for n in NUMBER.findall(text)]


def assert_same_render(output: str, expected: str) -> None:
    """Text identical outside numeric literals; numbers equal within ``REL``."""
    if output == expected:
        return

    assert NUMBER.split(output) == NUMBER.split(expected)
    assert _numbers(output) == pytest.approx(_numbers(expected), rel=REL)


NETWORKS: Dict[str, Callable[[], Network]] = {
    # neutral H/D/He with an energy equation (dE/dt from the .jfunc)
    "dedt": lambda: Network(str(FIXTURES / "test_jac_dedt.dat")),
    # charged species and electrons (He++), no radiation
    "cie_hepp": lambda: Network(str(FIXTURES / "react_cie_hepp.jet")),
    # photoionisation with one radiation band
    "h_photo": lambda: Network(
        str(REPO / "networks" / "h_photoionization" / "h_photo.jet"),
        radiation=RadiationConfig(
            bands=[13.6, float("inf")], profile_index=0, mode="nph"
        ),
    ),
}

ALL = list(NETWORKS)
CASES: Dict[str, List[str]] = {
    "directives.py": ALL,
    "repeat_lists.py": ALL,
    "repeat_expressions.py": ALL,
    "reduce.py": ALL,
    "language.cpp": ALL,
    "language.f90": ALL,
    "modifiers_thermal.py": ["dedt"],
    "modifiers_radiation.py": ["h_photo"],
}


@cache
def network(name: str) -> Network:
    return NETWORKS[name]()


def render(template: str, net: str, echo: bool = True) -> str:
    path = TEMPLATES / template
    return TemplateParser(network(net), path, echo_directives=echo).parse_file()


def _params() -> List[Tuple[str, str]]:
    return [(t, n) for t, nets in CASES.items() for n in nets]


@pytest.mark.parametrize("template, net", _params())
def test_render_matches_snapshot(template: str, net: str) -> None:
    snapshot = SNAPSHOTS / net / template
    output = render(template, net)
    if UPDATE:
        snapshot.parent.mkdir(parents=True, exist_ok=True)
        snapshot.write_text(output)

    assert snapshot.exists(), f"no snapshot {snapshot}; run JAFF_UPDATE_SNAPSHOTS=1"
    assert_same_render(output, snapshot.read_text())


def test_echo_off_leaves_no_directive_lines() -> None:
    for template, net in _params():
        output = render(template, net, echo=False)
        leaked = [line for line in output.splitlines() if "$JAFF" in line]
        assert not leaked, f"{template}/{net}: directive lines in output: {leaked}"


# Echo comment of each fixture's language, picked from its extension
ECHO_COMMENTS: Dict[str, str] = {".py": "#", ".cpp": "//", ".f90": "!"}


@pytest.mark.parametrize("template, net", _params())
def test_echo_only_adds_commented_directive_lines(template: str, net: str) -> None:
    echoed = render(template, net).splitlines()
    plain = render(template, net, echo=False).splitlines()

    # The plain render must be a subsequence of the echoed one ...
    extra = []
    remaining = iter(plain)
    expected = next(remaining, None)
    for line in echoed:
        if line == expected:
            expected = next(remaining, None)
        else:
            extra.append(line)
    assert expected is None, f"echoed render lost {expected!r}"

    # ... and everything added is a directive line behind the language's comment.
    comment = ECHO_COMMENTS[Path(template).suffix]
    assert extra
    assert all(line.lstrip().startswith(f"{comment} ") for line in extra), extra


def test_snapshot_set_matches_cases() -> None:
    expected = {(n, t) for t, nets in CASES.items() for n in nets}
    found = {(p.parent.name, p.name) for p in SNAPSHOTS.glob("*/*") if p.is_file()}
    assert found == expected, "stale or missing snapshots; run JAFF_UPDATE_SNAPSHOTS=1"


def test_every_template_has_cases() -> None:
    templates = {p.name for p in TEMPLATES.iterdir()}
    assert templates == set(CASES)


def test_same_render_tolerates_last_digit_noise_only() -> None:
    assert_same_render("k = 5.0572844545056314e-8*x;", "k = 5.0572844545056307e-8*x;")
    assert_same_render("y = 1.0d0", "y = 1.0D0")
    with pytest.raises(AssertionError):
        assert_same_render("k = 5.06e-8*x;", "k = 5.05e-8*x;")
    with pytest.raises(AssertionError):
        assert_same_render("k = 5.0e-8*x;", "k = 5.0e-8*y;")
