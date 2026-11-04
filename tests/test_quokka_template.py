# ABOUTME: The quokka template writes one C++ header, network.hpp, for Quokka's photochemistry solver:
# ABOUTME: species, bands and masses as constants, and the RHS and Jacobian over Quokka's 0-based solver vector

from pathlib import Path
from types import SimpleNamespace

import pytest

from jaff.cli import JaffGen

REPO = Path(__file__).resolve().parent.parent
H_PHOTO = REPO / "networks" / "h_photoionization" / "h_photo.jet"
CONFIG = """
[network.radiation]
bands = [13.6, 26.0]
profile_index = 0
mode = "nph"
rsl = "c_hat"
"""


@pytest.fixture(scope="module")
def generated(tmp_path_factory) -> str:
    """Render the quokka template for the hydrogen photoionization network once."""
    tmp = tmp_path_factory.mktemp("quokka")
    config = tmp / "jaffgen.toml"
    config.write_text(CONFIG)
    args = SimpleNamespace(
        network=str(H_PHOTO),
        config=str(config),
        label=None,
        funcfile=None,
        duplicate_policy=None,
        expand_nuclei=None,
        errors=None,
        network_config=None,
        outdir=str(tmp / "out"),
        indir=None,
        files=None,
        template="quokka",
    )
    JaffGen(args)
    return (tmp / "out" / "network.hpp").read_text()


def test_species_bands_and_masses(generated):
    assert "constexpr int NumSpec = 3;" in generated
    # identifiers keep the case of the species names: "+" -> "_p", the electron -> "e"
    assert "H = 0," in generated
    assert "H_p = 1," in generated
    assert "e = 2," in generated
    assert "constexpr int NumChemBands = 1;" in generated
    assert "return {{13.6, 26.0}};" in generated
    assert "constexpr double ChemBandsPowerLawIndex_ = 0;" in generated
    assert "spmasses[NumSpec] = {" in generated
    assert "gammas[NumSpec] = {" in generated


def test_rhs_and_jacobian_use_the_quokka_solver_vector(generated):
    assert "void network_rhs(State const &state, NetworkParams const &p, rosenbrock::Vector<NumEqs> &ydot)" in generated
    assert "void network_jac(State const &state, NetworkParams const &p, rosenbrock::Matrix<NumEqs> &J)" in generated
    assert "ydot[0] = " in generated  # 0-based rows
    assert "J(0, 0) = " in generated  # 0-based entries
    assert "state.xn[0]" in generated
    assert "state.rn[NumChemRadVarsPerGroup*0]" in generated
    assert "state.rn[NumChemRadVarsPerGroup*0+1]" in generated  # the flux factor
    assert "state.T" in generated and "state.c_hat" in generated


def test_no_microphysics(generated):
    for token in ("burn_t", "Array1D", "network_rp", "extern_parameters", "actual_rhs", "jacobian_utilities"):
        assert token not in generated, token
