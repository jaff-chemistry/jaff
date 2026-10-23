# ABOUTME: Unit tests for Network class initialization
# ABOUTME: Tests constructor, mass dict loading, and label handling

import os
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from jaff import Reactions, Species

# Add src to path for imports

from jaff import Network
from jaff.core.network.config import NetworkConfig
from tests._paths import FIXTURES_DIR


class TestNetworkInitialization:
    """Test Network class initialization functionality."""

    @pytest.fixture
    def fixtures_dir(self):
        """Return path to fixtures directory."""
        return str(FIXTURES_DIR)

    @pytest.fixture
    def sample_kida_file(self, fixtures_dir):
        """Return path to sample KIDA file."""
        return os.path.join(fixtures_dir, "sample_kida.dat")

    def test_successful_initialization(self, sample_kida_file, capsys):
        """Test successful network initialization with valid file."""
        # Suppress print output during initialization
        with patch("builtins.print"):
            network = Network(sample_kida_file)

        # Check basic attributes are initialized
        assert network.filename == Path(sample_kida_file).resolve()
        assert network.label == "sample_kida"
        assert isinstance(network.species, Species)
        assert isinstance(network.reactions, Reactions)
        assert isinstance(network.mass_dict, dict)

        # Check that some reactions were loaded
        assert len(network.reactions) > 0
        assert len(network.species) > 0

    def test_initialization_with_custom_label(self, sample_kida_file):
        """Test initialization with custom label parameter."""
        custom_label = "my_custom_network"

        with patch("builtins.print"):
            network = Network(sample_kida_file, label=custom_label)

        assert network.label == custom_label
        assert network.filename == Path(sample_kida_file).resolve()

    def test_initialization_with_nonexistent_file(self):
        """Test initialization fails gracefully with non-existent file."""
        with pytest.raises(FileNotFoundError):
            Network(os.path.join("path", "to", "nonexistent", "file.dat"))

    def test_default_label_extraction(self, sample_kida_file):
        """Test default label extraction from filename."""
        with patch("builtins.print"):
            network = Network(sample_kida_file)

        # Should extract 'sample_kida' from 'sample_kida.dat'
        assert network.label == "sample_kida"

        # Test with path containing directories
        test_path = os.path.join("some", "long", "path", "to", "network_file.txt")
        with (
            patch("builtins.print"),
            patch("builtins.open", MagicMock()),
            patch("pathlib.Path.exists", return_value=True),
        ):
            with (
                patch(
                    "jaff.core.network.network.load_network_config",
                    return_value=NetworkConfig(),
                ),
                patch("jaff.core.network._spec.NetworkSpec._load_aux_funcs", return_value={}),
                patch.object(Network, "_Network__load_network", MagicMock()),
            ):
                with patch.object(Network, "check_sink_sources", MagicMock()):
                    with patch.object(Network, "check_recombinations", MagicMock()):
                        with patch.object(Network, "check_isomers", MagicMock()):
                            with patch.object(
                                Network,
                                "_Network__generate_reaction_matrices",
                                MagicMock(),
                            ):
                                network = Network(test_path)

        assert network.label == "network_file"

    def test_mass_dict_loading(self, fixtures_dir):
        """Test mass dictionary loading"""
        fpath = os.path.join(fixtures_dir, "empty_network.dat")
        net = Network(fpath)

        assert net.mass_dict["H"]["mass"] == pytest.approx(1.673773e-24)
        assert net.mass_dict["He"]["mass"] == pytest.approx(6.646473e-24)
        assert net.mass_dict["C"]["mass"] == pytest.approx(1.994473e-23)
        assert net.mass_dict["O"]["mass"] == pytest.approx(2.656763e-23)

    def test_errors_parameter_true(self):
        """Test initialization with errors=True parameter."""
        # Use a valid file that shouldn't trigger errors
        valid_file = FIXTURES_DIR / "sample_kida_valid.dat"

        with patch("builtins.print"):
            # Mock sys.exit to prevent test from exiting
            with patch("sys.exit") as mock_exit:
                network = Network(str(valid_file), strict=True)

                # If there are any validation errors, sys.exit should be called
                # In our valid sample file, we don't expect errors, so it shouldn't exit
                mock_exit.assert_not_called()

    def test_errors_parameter_false(self, sample_kida_file):
        """Test initialization with errors=False parameter (default)."""
        with patch("builtins.print"):
            network = Network(sample_kida_file, strict=False)

        # Should complete without raising exceptions
        assert network is not None

    def test_initialization_workflow(self, sample_kida_file):
        """Test that all initialization steps are called in correct order."""
        with patch("builtins.print"):
            with patch.object(Network, "_Network__load_network") as mock_load:
                with patch.object(Network, "check_sink_sources") as mock_sink:
                    with patch.object(Network, "check_recombinations") as mock_recomb:
                        with patch.object(Network, "check_isomers") as mock_isomers:
                            with patch.object(
                                Network, "_Network__generate_reaction_matrices"
                            ) as mock_matrices:
                                network = Network(sample_kida_file, strict=True)

        # Verify all methods were called
        mock_load.assert_called_once_with()
        mock_sink.assert_called_once_with(True)
        mock_recomb.assert_called_once_with(True)
        mock_isomers.assert_called_once_with(True)
        mock_matrices.assert_called_once()

    def test_empty_network_file(self, fixtures_dir):
        """Test initialization with empty network file."""
        empty_file = os.path.join(fixtures_dir, "empty_network.dat")

        with patch("builtins.print"):
            network = Network(empty_file)

        # Should initialize but with no reactions
        assert len(network.reactions) == 0
        # May have default species or none
        assert isinstance(network.species, Species)

    def test_initial_data_structures(self, sample_kida_file):
        """Test that data structures are properly initialized."""
        with patch("builtins.print"):
            network = Network(sample_kida_file)

        # Check data structure types and initial states
        assert isinstance(network.species, Species)
        assert isinstance(network.reactions, Reactions)
        assert isinstance(network.mass_dict, dict)

        # Check that reactant_matrix and product_matrix are generated
        assert network.reactant_matrix is not None
        assert network.product_matrix is not None
        assert hasattr(network.reactant_matrix, "shape")  # Should be numpy array
        assert hasattr(network.product_matrix, "shape")  # Should be numpy array
