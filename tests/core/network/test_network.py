# ABOUTME: Tests for Network construction, metadata, matrices, lookups, and robustness
# ABOUTME: Public-behaviour only - no patching of private methods

from pathlib import Path

import pytest

from jaff import Network, Reactions, Species


class TestConstruction:
    """Basic construction, labelling, and initialised attributes."""

    def test_attributes_initialised(self, fixtures_dir):
        net = Network(str(fixtures_dir / "sample_kida.dat"))
        assert net.filepath == (fixtures_dir / "sample_kida.dat").resolve()
        assert net.label == "sample_kida"
        assert isinstance(net.species, Species)
        assert isinstance(net.reactions, Reactions)
        assert net.reactions.count > 0 and net.species.count > 0

    def test_custom_label(self, fixtures_dir):
        net = Network(str(fixtures_dir / "sample_kida.dat"), label="custom")
        assert net.label == "custom"

    def test_label_defaults_to_file_stem(self, make_network):
        net = make_network(["H + H -> H2 [10,1000] 1e-10"], name="my_net.dat")
        assert net.label == "my_net"

    def test_nonexistent_file_raises(self):
        with pytest.raises(FileNotFoundError):
            Network(str(Path("no") / "such" / "file.dat"))


class TestEmptyNetwork:
    """A file with no reactions still yields a well-formed Network."""

    def test_no_reactions(self, fixtures_dir):
        net = Network(str(fixtures_dir / "empty_network.dat"))
        assert net.reactions.count == 0
        assert isinstance(net.species, Species)

    def test_matrices_have_zero_rows(self, make_network):
        net = make_network(["# comments only", "# no reactions"])
        assert net.reactant_matrix is not None
        assert net.product_matrix is not None
        assert net.reactant_matrix.shape[0] == 0


class TestLookups:
    """Missing-key lookups raise KeyError."""

    def test_missing_species_name(self, sample_network):
        with pytest.raises(KeyError):
            sample_network.species["NONEXISTENT"]

    def test_missing_reaction_serialized(self, sample_network):
        with pytest.raises(KeyError):
            sample_network.reactions["NONEXISTENT_SERIALIZED"]

    def test_missing_species_from_serialized(self, sample_network):
        with pytest.raises(KeyError):
            sample_network.species.from_serialized("NONEXISTENT_SERIALIZED")


class TestRobustness:
    """The loader tolerates awkward but valid inputs."""

    def test_long_species_name(self, make_network):
        long_name = "C10H20O5N3S2P1"
        net = make_network([f"H + {long_name} -> H2 + {long_name} [10,1000] 1e-10"])
        assert long_name in {s.name for s in net.species}

    def test_charged_and_organic_species(self, make_network):
        net = make_network(
            [
                "H+ + e- -> H [10,1000] 1e-12",
                "C2H5OH + OH -> C2H4OH + H2O [10,1000] 1e-11",
                "H3O+ + NH3 -> NH4+ + H2O [10,1000] 1e-9",
            ]
        )
        names = {s.name for s in net.species}
        assert {"H+", "e-", "C2H5OH", "H3O+", "NH4+"} <= names

    def test_unicode_in_comments(self, make_network):
        net = make_network(
            [
                "# Greek letters: α + β → γ",
                "H + H -> H2 [10,1000] 1e-10",
            ]
        )
        assert net.reactions.count == 1

    def test_circular_dependencies(self, make_network):
        net = make_network(
            [
                "H -> H+ + e- [10,1000] 1e-10",
                "H+ + e- -> H [10,1000] 1e-12",
                "H + H -> H2 [10,1000] 1e-15",
            ]
        )
        assert net.reactions.count == 3

    def test_moderately_large_distinct_network(self, make_network):
        import itertools

        elements = [
            "H",
            "He",
            "C",
            "N",
            "O",
            "F",
            "Ne",
            "Na",
            "Mg",
            "Al",
            "Si",
            "P",
            "S",
            "Cl",
            "Ar",
            "K",
        ]
        pairs = list(itertools.combinations(elements, 2))[:50]
        lines = [f"{a} + {b} -> {a}{b} [10,1000] 1e-10" for a, b in pairs]
        net = make_network(lines)
        assert net.reactions.count == 50
        assert net.reactant_matrix.shape[0] == 50
        assert net.product_matrix.shape[0] == 50


class TestStoichiometryMatrices:
    """``reactant_matrix`` / ``product_matrix``: built on first use, cached, read-only."""

    LINES = [
        "H + H -> H2 [10,1000] 1e-10",
        "H2 + H+ -> H2+ + H [10,1000] 1e-9",
    ]

    def test_counts_per_reaction_and_species(self, make_network) -> None:
        net = make_network(self.LINES)
        h, h2 = net.species["H"].index, net.species["H2"].index

        assert net.reactant_matrix.shape == (2, net.species.count)
        assert net.reactant_matrix[0, h] == 2
        assert net.product_matrix[0, h2] == 1
        assert net.product_matrix[1, h] == 1
        assert net.reactant_matrix.dtype == int

    def test_not_built_until_used(self, make_network) -> None:
        net = make_network(self.LINES)

        assert "reactant_matrix" not in vars(net)
        assert "product_matrix" not in vars(net)

    def test_cached(self, make_network) -> None:
        net = make_network(self.LINES)

        assert net.reactant_matrix is net.reactant_matrix
        assert net.product_matrix is net.product_matrix

    @pytest.mark.parametrize("name", ["reactant_matrix", "product_matrix"])
    def test_read_only(self, make_network, name) -> None:
        matrix = getattr(make_network(self.LINES), name)

        with pytest.raises(ValueError, match="read-only"):
            matrix[0, 0] += 1
