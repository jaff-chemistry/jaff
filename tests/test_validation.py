# ABOUTME: Behavioural tests for Network validation diagnostics
# ABOUTME: Assert on Sink/Source/Recombination/Isomer reports, not on log text

import pytest


class TestSinkSources:
    """Species never produced (sink) or never consumed (source)."""

    def test_balanced_network_is_ok(self, make_network):
        net = make_network(
            [
                "H + H -> H2 [10,1000] 1e-10",
                "H2 -> H + H [10,1000] 1e-15",
            ]
        )
        report = net.check_sink_sources(errors=False)
        assert report.ok
        assert report.sinks == frozenset()
        assert report.sources == frozenset()

    def test_sink_detected(self, make_network):
        # He appears only as a reactant -> sink, never a source.
        net = make_network(
            [
                "H + He -> H2 [10,1000] 1e-10",
                "H2 -> H + H [10,1000] 1e-15",
            ]
        )
        report = net.check_sink_sources(errors=False)
        assert "He" in report.sinks
        assert "He" not in report.sources
        assert not report.ok

    def test_source_detected(self, make_network):
        # He appears only as a product -> source, never a sink.
        net = make_network(
            [
                "H + H -> H2 + He [10,1000] 1e-10",
                "H2 -> H + H [10,1000] 1e-15",
            ]
        )
        report = net.check_sink_sources(errors=False)
        assert "He" in report.sources
        assert "He" not in report.sinks

    def test_dummy_species_ignored(self, make_network):
        net = make_network(
            [
                "H + H -> H2 + dummy [10,1000] 1e-10",
                "dummy + O -> O [10,1000] 1e-15",
            ]
        )
        report = net.check_sink_sources(errors=False)
        flagged = report.sinks | report.sources
        assert not any("dummy" in name.lower() for name in flagged)

    def test_errors_true_exits(self, make_network):
        # He=sink, Ne=source: with errors=True the loader aborts.
        with pytest.raises(SystemExit):
            make_network(["He -> Ne [10,1000] 1e-10"], strict=True)


class TestRecombinations:
    """Cations lacking an electron recombination reaction."""

    def test_proper_recombination_is_ok(self, make_network):
        net = make_network(
            [
                "H -> H+ + e- [10,1000] 1e-10",
                "H+ + e- -> H [10,1000] 1e-12",
            ]
        )
        assert net.check_recombinations(errors=False).ok

    def test_missing_recombination_flagged(self, make_network):
        net = make_network(
            [
                "H -> H+ + e- [10,1000] 1e-10",
                "C+ + H2 -> CH+ + H [10,1000] 1e-11",  # no recombination for C+
            ]
        )
        report = net.check_recombinations(errors=False)
        assert "C+" in report.missing

    def test_errors_true_exits(self, make_network):
        net = make_network(
            [
                "H -> H+ + e- [10,1000] 1e-10",
                "C+ + H2 -> CH+ + H [10,1000] 1e-11",
            ]
        )
        with pytest.raises(SystemExit):
            net.check_recombinations(errors=True)


class TestIsomers:
    """Species sharing an atomic composition (e.g. HCO / HOC)."""

    def test_no_isomers_is_ok(self, make_network):
        net = make_network(
            [
                "H + H -> H2 [10,1000] 1e-10",
                "C + O -> CO [10,1000] 1e-11",
            ]
        )
        assert net.check_isomers(errors=False).ok

    def test_isomers_detected(self, make_network):
        net = make_network(
            [
                "C + O -> CO [10,1000] 1e-11",
                "H + CO -> HCO [10,1000] 1e-11",
                "H + CO -> HOC [10,1000] 1e-11",  # same composition as HCO
            ]
        )
        report = net.check_isomers(errors=False)
        assert report.groups == (("HCO", "HOC"),)
        assert not report.ok

    def test_errors_true_exits(self, make_network):
        net = make_network(
            [
                "C + O -> CO [10,1000] 1e-11",
                "H + CO -> HCO [10,1000] 1e-11",
                "H + CO -> HOC [10,1000] 1e-11",
            ]
        )
        with pytest.raises(SystemExit):
            net.check_isomers(errors=True)


class TestValidationRunsDuringLoad:
    """A clean fixture network loads without triggering the error exit."""

    def test_valid_network_loads_with_errors_true(self, fixtures_dir):
        from jaff import Network

        # sample_kida_valid.dat is balanced; errors=True must not exit.
        Network(str(fixtures_dir / "sample_kida_valid.dat"), strict=True)
