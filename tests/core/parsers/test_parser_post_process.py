# ABOUTME: Tests for the Parser.post_process hook: abstract contract, which parsers
# ABOUTME: are called, and frozen reaction/species catalogues while it runs

from unittest.mock import patch

import pytest

from jaff.core.parsers.network._formats._parser import Parser
from jaff.core.parsers.network._formats.kida.parser import KidaParser
from jaff.core.parsers.network._formats.prizmo.parser import PrizmoParser
from jaff.core.species import Specie

NETWORK = "H + H -> H2 [] 1e-17\n"


def test_post_process_is_abstract():
    class Incomplete(Parser):
        name = "incomplete"
        priority = 0
        handlers = []

        def process(self, records):
            return None

    with pytest.raises(TypeError, match="post_process"):
        Incomplete()


def test_hook_runs_for_contributing_parser(make_network):
    seen = []

    def spy(self, network):
        seen.append(network.reactions.count)

    with patch.object(PrizmoParser, "post_process", spy):
        make_network(NETWORK, name="r.jet")

    assert seen == [1]


def test_hook_not_called_for_parsers_without_records(make_network):
    calls = []

    with patch.object(KidaParser, "post_process", lambda self, network: calls.append(1)):
        make_network(NETWORK, name="r.jet")

    assert calls == []


def test_catalogues_frozen_during_hook(make_network):
    errors = []

    def spy(self, network):
        for add in (
            lambda: network.reactions.add(network.reactions[0]),
            lambda: network.species.add(Specie("He", 99)),
        ):
            try:
                add()
            except RuntimeError as e:
                errors.append(str(e))

    with patch.object(PrizmoParser, "post_process", spy):
        make_network(NETWORK, name="r.jet")

    assert len(errors) == 2
    assert all("frozen" in e for e in errors)


def test_hook_runs_before_rates_are_standardised(make_network):
    seen = []

    def spy(self, network):
        seen.append(str(network.reactions[0].rate_segments[0].rate))

    with patch.object(PrizmoParser, "post_process", spy):
        make_network("H + H -> H2 [] 1e-10*n_H\n", name="r.jet")

    # Still the raw symbol: normalisation (n_H -> nden[i]) has not run yet.
    assert "n_h" in seen[0]


def test_elements_available_during_hook(make_network):
    seen = []

    def spy(self, network):
        seen.append(network.elements.count)

    with patch.object(PrizmoParser, "post_process", spy):
        make_network(NETWORK, name="r.jet")

    assert seen == [1]
