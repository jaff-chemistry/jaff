# ABOUTME: Tests for network file-format detection and format-specific parsing
# ABOUTME: Covers KIDA, UDFA, PRIZMO, KROME, UCLCHEM plus comment/malformed handling

import pytest

from jaff import Network
from jaff.errors import ParserError

# Each bundled sample exercises one auto-detected input format.
SAMPLE_FORMATS = ["kida", "udfa", "prizmo", "krome", "uclchem"]


@pytest.fixture
def load_sample(fixtures_dir):
    """Return a loader for a bundled ``sample_<fmt>.dat`` fixture."""

    def _load(fmt):
        return Network(str(fixtures_dir / f"sample_{fmt}.dat"))

    return _load


@pytest.mark.parametrize("fmt", SAMPLE_FORMATS)
def test_format_auto_detected_and_parsed(load_sample, fmt):
    """Every supported format auto-detects and yields at least one reaction."""
    net = load_sample(fmt)
    assert net.reactions.count > 0


class TestFormatSpecifics:
    """Format-specific parsing details worth pinning down."""

    def test_kida_multibody_reaction_and_limits(self, load_sample):
        # H + H + H -> H2 + H, bounded [10, 1000].
        net = load_sample("kida")
        match = [
            r
            for r in net.reactions
            if r.reactants.count == 3
            and all(x.name == "H" for x in r.reactants)
            and r.products.count == 2
            and "H2" in r.products
        ]
        assert match, "expected H + H + H -> H2 + H"
        assert (match[0].tmin, match[0].tmax) == (10, 1000)

    def test_udfa_photodissociation_carries_av(self, load_sample):
        net = load_sample("udfa")
        match = [
            r
            for r in net.reactions
            if r.reactants.core.count == 1
            and r.reactants.core[0].name == "H2"
            and r.products.core.count == 2
            and "H" in r.products
        ]
        assert match, "expected H2 -> H + H photodissociation"
        assert (match[0].tmin, match[0].tmax) == (10, 3000)
        assert "av" in str(match[0].rate_coefficient)

    def test_prizmo_variable_substitution(self, load_sample):
        # y = tgas/300 is substituted into the O + H -> OH coefficient.
        net = load_sample("prizmo")
        match = [
            r
            for r in net.reactions
            if r.reactants.count == 2 and "O" in r.reactants and "H" in r.reactants
        ]
        assert match, "expected O + H -> OH"
        rate = str(match[0].rate_coefficient)
        assert "tgas" in rate.lower()
        assert "8.648" in rate or "8.65" in rate  # 9.9e-11 * 300**0.38

    def test_prizmo_cosmic_ray_reaction_present(self, load_sample):
        net = load_sample("prizmo")
        assert any(
            r.reactants.count == 2 and "H2" in r.reactants and "_CR" in r.reactants
            for r in net.reactions
        ), "expected H2 + CR reaction"

    def test_prizmo_photo_reaction_typed(self, load_sample):
        net = load_sample("prizmo")
        photo = [r for r in net.reactions if "photorates" in str(r.rate_coefficient)]
        assert photo, "expected a photochemistry reaction"
        assert all(r.type == "photo" for r in photo)

    def test_krome_variable_substitution(self, load_sample):
        # inv_tgas @var is expanded in C+ + H2 -> CH+ + H.
        net = load_sample("krome")
        match = [
            r
            for r in net.reactions
            if r.reactants.count == 2 and "C+" in r.reactants and "H2" in r.reactants
        ]
        assert match, "expected C+ + H2 reaction"
        rate = str(match[0].rate_coefficient).lower()
        assert "tgas" in rate and "exp" in rate

    def test_uclchem_species_parsed_despite_nan_fields(self, load_sample):
        net = load_sample("uclchem")
        names = {s.name for s in net.species}
        assert {"H", "H2", "OH"} <= names


class TestParserRobustness:
    """Comments, empty files, malformed lines, mixed-format inputs."""

    def test_empty_fixture_has_no_reactions(self, fixtures_dir):
        net = Network(str(fixtures_dir / "empty_network.dat"))
        assert net.reactions.count == 0

    def test_only_comments(self, make_network):
        net = make_network(
            [
                "# comment only",
                "! KIDA style comment",
                "# nothing to parse",
            ]
        )
        assert net.reactions.count == 0

    def test_malformed_fixture_raises_parsererror(self, fixtures_dir):
        # A line with no rate expression is a hard parse error, not a skip.
        with pytest.raises(ParserError):
            Network(str(fixtures_dir / "malformed_network.dat"))

    def test_mixed_format_lines_both_parsed(self, make_network):
        net = make_network(
            [
                "H + H -> H2 [10,1000] 1e-10",
                "1:RR:H:e-:H-::::1:1e-16:0:0:10:10000",
            ]
        )
        assert net.reactions.count >= 2

    def test_krome_shortcuts_expanded(self, make_network):
        net = make_network(
            [
                "@format:idx,R,R,P,P,tmin,tmax,rate",
                "1,H+,e-,H,,10,10000,2.59e-13*invte",
                "2,O,H,OH,,10,41000,9.9e-11*t32**(-0.38)",
            ]
        )
        for r in net.reactions:
            rate = str(r.rate_coefficient)
            assert "invte" not in rate and "t32" not in rate
            assert "tgas" in rate.lower()


def test_species_are_consistent_with_reactions(sample_network):
    """Every reactant/product core species is registered in the catalogue."""
    names = {s.name for s in sample_network.species}
    for r in sample_network.reactions:
        for sp in list(r.reactants.core) + list(r.products.core):
            assert sp.name in names
