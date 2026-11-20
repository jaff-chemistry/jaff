# ABOUTME: Codegen reports each optimization stage on the progress display, in order
# ABOUTME: (a stand-in jaff_progress records the descriptions)

from contextlib import contextmanager

import pytest

import jaff.codegen.codegen as codegen_module
from jaff.codegen import Codegen

NETWORK = [
    "H + H -> H2 [10,1e4] 1e-10*tgas**0.7",
    "H2 + H -> H + H + H [10,1e4] 2e-10*tgas**(-1.3)*exp(-5e4/tgas)",
]


class Recorder:
    """Stands in for jaff_progress: records each description, shows nothing."""

    def __init__(self) -> None:
        self.descriptions: list[str] = []

    @contextmanager
    def indeterminate(self, description: str):
        self.descriptions.append(description)
        yield None

    def track(self, sequence, *args, description: str = "Working...", **kwargs):
        self.descriptions.append(description)
        yield from sequence


@pytest.fixture
def recorder(monkeypatch):
    recorder = Recorder()
    monkeypatch.setattr(codegen_module, "jaff_progress", recorder)

    return recorder


def test_cse_stages_are_reported_in_order(make_network, recorder):
    Codegen(make_network(NETWORK), "c++").cse_group([("rhs", {}), ("jacobian", {})])
    stages = [
        "Rewriting expressions",
        "Sharing logarithms",
        "Finding common subexpressions",
    ]
    positions = [recorder.descriptions.index(stage) for stage in stages]

    assert positions == sorted(positions)
    assert recorder.descriptions.index("Generating RHS code") > positions[-1]


def test_without_cse_only_the_rewrite_is_reported(make_network, recorder):
    Codegen(make_network(NETWORK), "c++").get_indexed_rates(use_cse=False)

    assert "Rewriting expressions" in recorder.descriptions
    assert "Sharing logarithms" not in recorder.descriptions
    assert "Finding common subexpressions" not in recorder.descriptions
