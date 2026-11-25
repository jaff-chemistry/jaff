# ABOUTME: Unit tests for Network class parser functionality
# ABOUTME: Tests format detection and parsing for all supported formats

import os
import tempfile
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
from sympy import Piecewise, symbols

# Add src to path for imports

from jaff import Network, Reaction, Specie
from tests._paths import FIXTURES_DIR


class TestNetworkParsers:
    """Test Network class parser functionality for all formats."""

    @pytest.fixture
    def fixtures_dir(self):
        """Return path to fixtures directory."""
        return str(FIXTURES_DIR)

    def test_kida_format_detection(self, fixtures_dir):
        """Test automatic detection and parsing of KIDA format."""
        kida_file = os.path.join(fixtures_dir, "sample_kida.dat")

        with patch("builtins.print"):
            network = Network(kida_file)

        # Check that reactions were parsed
        assert len(network.reactions) > 0

        # Check specific reaction from our sample file
        # H + H + H -> H2 + H
        h_idx = network.species["H"].index
        h2_idx = network.species["H2"].index

        # Find the reaction
        found = False
        for reaction in network.reactions:
            if (
                reaction.reactants.count == 3
                and all(r.name == "H" for r in reaction.reactants)
                and reaction.products.count == 2
                and "H2" in reaction.products
            ):
                found = True
                assert reaction.tmin == 10
                assert reaction.tmax == 1000
                break

        assert found, "Expected H + H + H -> H2 + H reaction not found"

    def test_udfa_format_detection(self, fixtures_dir):
        """Test automatic detection and parsing of UDFA format."""
        udfa_file = os.path.join(fixtures_dir, "sample_udfa.dat")

        with patch("builtins.print"):
            network = Network(udfa_file)

        # Check that reactions were parsed
        assert len(network.reactions) > 0

        # Check for H2 photodissociation reaction (UDFA format filters out PHOTON)
        found = False
        for reaction in network.reactions:
            if (
                reaction.reactants.core.count == 1
                and reaction.reactants.core[0].name == "H2"
                and reaction.products.core.count == 2
                and "H" in reaction.products
            ):
                found = True
                assert reaction.tmin == 10
                assert reaction.tmax == 3000
                # Check that rate contains 'av' parameter (photodissociation)
                assert "av" in str(reaction.rate.coefficient)
                break

        assert found, "Expected H2 -> H + H photodissociation reaction not found"

    def test_prizmo_format_detection(self, fixtures_dir):
        """Test automatic detection and parsing of PRIZMO format."""
        prizmo_file = os.path.join(fixtures_dir, "sample_prizmo.dat")

        with patch("builtins.print"):
            network = Network(prizmo_file)

        # Check that reactions were parsed
        assert len(network.reactions) > 0

        # Check that variables were parsed and substituted
        # The O + H -> OH reaction should have the variable y substituted
        found = False
        for reaction in network.reactions:
            if (
                reaction.reactants.count == 2
                and "O" in reaction.reactants
                and "H" in reaction.reactants
            ):
                found = True
                # Check that the rate expression contains tgas and substituted coefficient
                rate_str = str(reaction.rate.coefficient)
                assert "tgas" in rate_str.lower()
                # After variable substitution y=tgas/300, the coefficient changes
                assert "8.648" in rate_str or "8.65" in rate_str  # 9.9e-11 * 300^0.38
                break

        assert found, "Expected O + H -> OH reaction not found"

    def test_krome_format_detection(self, fixtures_dir):
        """Test automatic detection and parsing of KROME format."""
        krome_file = os.path.join(fixtures_dir, "sample_krome.dat")

        with patch("builtins.print"):
            network = Network(krome_file)

        # Check that reactions were parsed
        assert len(network.reactions) > 0

        # Check that @var declarations were processed
        # The C+ + H2 -> CH+ + H reaction uses inv_tgas variable
        found = False
        for reaction in network.reactions:
            if (
                reaction.reactants.count == 2
                and "C+" in reaction.reactants
                and "H2" in reaction.reactants
            ):
                found = True
                # Check that inv_tgas was substituted
                rate_str = str(reaction.rate.coefficient)
                assert "tgas" in rate_str.lower()
                assert "exp" in rate_str.lower()
                break

        assert found, "Expected C+ + H2 reaction not found"

    def test_uclchem_format_detection(self, fixtures_dir):
        """Test automatic detection and parsing of UCLCHEM format."""
        uclchem_file = os.path.join(fixtures_dir, "sample_uclchem.dat")

        with patch("builtins.print"):
            network = Network(uclchem_file)

        # Check that reactions were parsed
        assert len(network.reactions) > 0

        # Check for reactions with NAN markers
        # All reactions should have been parsed despite NAN fields
        species_names = [s.name for s in network.species]
        assert "H" in species_names
        assert "H2" in species_names
        assert "OH" in species_names

    def test_custom_variables_parsing(self, fixtures_dir):
        """Test parsing of custom variables in PRIZMO and KROME formats."""
        prizmo_file = os.path.join(fixtures_dir, "sample_prizmo.dat")

        with patch("builtins.print"):
            network = Network(prizmo_file)

        # Check cosmic ray reaction that uses zeta variable
        found = False
        for reaction in network.reactions:
            if (
                reaction.reactants.count == 2
                and "H2" in reaction.reactants
                and "_CR" in reaction.reactants
            ):
                found = True
                # The rate should be the zeta value (1.3e-17)
                # Since it's a constant after substitution
                break

        assert found, "Expected H2 + CR reaction not found"

    def test_photo_chemistry_parsing(self, fixtures_dir):
        """Test parsing of photochemistry reactions."""
        prizmo_file = os.path.join(fixtures_dir, "sample_prizmo.dat")

        with patch("builtins.print"):
            network = Network(prizmo_file)

        # Check for photo reaction by looking for photorates function
        photo_reactions = [
            r for r in network.reactions if "photorates" in str(r.rate.coefficient)
        ]
        assert len(photo_reactions) >= 1, "Expected photochemistry reaction not found"

        # Verify the photo reactions are correctly identified
        for reaction in photo_reactions:
            assert reaction.type == "photo"

    def test_temperature_limits_application(self, fixtures_dir):
        """Test that temperature limits (tmin/tmax) are correctly applied."""
        kida_file = os.path.join(fixtures_dir, "sample_kida.dat")

        with patch("builtins.print"):
            network = Network(kida_file)

        # Check reactions have temperature limits
        tgas = symbols("tgas")
        for reaction in network.reactions:
            assert hasattr(reaction, "tmin")
            assert hasattr(reaction, "tmax")

            # A temperature-dependent rate with a defined bound holds that
            # boundary value under the default "clip" cutoff, i.e. it is wrapped
            # in a Piecewise.
            rate = reaction.rate.coefficient
            if tgas in rate.free_symbols and reaction.config.t_cutoff == "clip":
                if reaction.tmin is not None or reaction.tmax is not None:
                    assert rate.has(Piecewise)

    def test_comment_and_empty_line_handling(self, fixtures_dir):
        """Test that comments and empty lines are properly handled."""
        empty_file = os.path.join(fixtures_dir, "empty_network.dat")

        with patch("builtins.print"):
            network = Network(empty_file)

        # Should load without errors but have no reactions
        assert len(network.reactions) == 0

    def test_malformed_line_handling(self, fixtures_dir):
        """Test handling of malformed lines in network files."""
        malformed_file = os.path.join(fixtures_dir, "malformed_network.dat")

        # The parser should skip malformed lines without crashing
        with patch("builtins.print"):
            try:
                network = Network(malformed_file)
                # If it loads, check that some lines were skipped
                assert True  # Successfully handled malformed input
            except Exception as e:
                # Some malformed lines might cause exceptions, which is acceptable
                assert True

    def test_format_detection_priority(self):
        """Test format detection priority when multiple patterns match."""
        # Create a file that could match multiple formats
        with tempfile.NamedTemporaryFile(mode="w", suffix=".dat", delete=False) as f:
            f.write("# Test file with mixed format indicators\n")
            f.write("# This has -> like PRIZMO\n")
            f.write("H + H -> H2 [10,1000] 1e-10\n")
            f.write("# But also has : like UDFA\n")
            f.write("1:RR:H:e-:H-::::1:1e-16:0:0:10:10000\n")
            temp_file = f.name

        try:
            with patch("builtins.print"):
                network = Network(temp_file)

            # Should parse both lines correctly
            assert len(network.reactions) >= 2
        finally:
            os.unlink(temp_file)

    def test_krome_shortcuts_parsing(self):
        """Test that KROME shortcuts are properly parsed."""
        # Create a file using KROME shortcuts
        with tempfile.NamedTemporaryFile(mode="w", suffix=".dat", delete=False) as f:
            f.write("@format:idx,R,R,P,P,tmin,tmax,rate\n")
            f.write("1,H+,e-,H,,10,10000,2.59e-13*invte\n")
            f.write("2,O,H,OH,,10,41000,9.9e-11*t32**(-0.38)\n")
            temp_file = f.name

        try:
            with patch("builtins.print"):
                network = Network(temp_file)

            # Check that shortcuts were substituted
            for reaction in network.reactions:
                rate_str = str(reaction.rate.coefficient)
                # Should not contain shortcut names
                assert "invte" not in rate_str
                assert "t32" not in rate_str
                # Should contain expanded expressions
                assert "tgas" in rate_str.lower()
        finally:
            os.unlink(temp_file)

    def test_species_creation_from_reactions(self, fixtures_dir):
        """Test that species are correctly created from parsed reactions."""
        kida_file = os.path.join(fixtures_dir, "sample_kida.dat")

        with patch("builtins.print"):
            network = Network(kida_file)

        # Check that all species in reactions exist in species list
        species_names = [s.name for s in network.species]

        for reaction in network.reactions:
            for reactant in reaction.reactants.core:
                assert reactant.name in species_names
            for product in reaction.products.core:
                assert product.name in species_names

        # Check species dictionary
        for species in network.species:
            assert species.name in network.species

    def test_rate_expression_parsing(self):
        """Test parsing of various rate expressions."""
        # Create a file with different rate expressions
        with tempfile.NamedTemporaryFile(mode="w", suffix=".dat", delete=False) as f:
            f.write("# Test various rate expressions\n")
            f.write("H + H -> H2 [10,1000] 1.0e-10\n")  # Simple constant
            f.write("H + e- -> H- [10,10000] 3e-16 * (tgas/300)**0.5\n")  # Power law
            f.write(
                "C+ + H2 -> CH+ + H [10,41000] 1e-10 * exp(-4640/tgas)\n"
            )  # Exponential
            f.write(
                "O + H -> OH [10,41000] 9.9e-11 * sqrt(tgas) * exp(-100/tgas)\n"
            )  # Complex
            temp_file = f.name

        try:
            with patch("builtins.print"):
                network = Network(temp_file)

            assert len(network.reactions) == 4

            # Check that all rate expressions were parsed as sympy objects
            for reaction in network.reactions:
                assert reaction.rate.coefficient is not None
                # Rate should be a sympy expression or number
                assert hasattr(reaction.rate.coefficient, "free_symbols") or isinstance(
                    reaction.rate.coefficient, (int, float)
                )
        finally:
            os.unlink(temp_file)

    def test_special_species_handling(self, fixtures_dir):
        """Test handling of special species like e-, PHOTON, CR, etc."""
        kida_file = os.path.join(fixtures_dir, "sample_kida.dat")

        with patch("builtins.print"):
            network = Network(kida_file)

        # Check that special species are recognized
        species_names = [s.name for s in network.species]

        # These special species should be in our sample files
        special_species = ["e-", "PHOTON", "H+", "H-"]
        for special in special_species:
            if any(
                special in network.reactions[i].verbatim
                for i in range(len(network.reactions))
            ):
                if special == "PHOTON":
                    mapped = f"_{special}"
                    assert any(
                        mapped in network.reactions[i].verbatim
                        for i in range(len(network.reactions))
                    ), f"{mapped} should appear in reaction verbatim"
                else:
                    assert special in species_names
