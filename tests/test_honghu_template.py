# ABOUTME: The honghu template writes a Julia network for Honghu.jl: one module per network,
# ABOUTME: one-based indices, no SymPy infinity, and Julia that parses (when julia is installed)

import re
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

from jaff.cli import JaffGen

REPO = Path(__file__).resolve().parent.parent
H_PHOTO = REPO / "networks" / "h_photoionization" / "h_photo.jet"
CONFIG = """
[network]
label = "HPhoto"

[network.radiation]
bands = [13.6, "inf"]
profile_index = 0
mode = "nph"
rsl = "c_hat"
"""


@pytest.fixture(scope="module")
def generated(tmp_path_factory) -> str:
    """Render the honghu template for the hydrogen photoionization network once."""
    tmp = tmp_path_factory.mktemp("honghu")
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
        template="honghu",
    )
    JaffGen(args)
    return (tmp / "out" / "network.jl").read_text()


def test_network_is_a_module_named_after_the_label(generated):
    assert "module JAFF_HPhoto" in generated
    assert "struct HPhoto{T<:AbstractFloat} <: Honghu.AbstractNetwork" in generated
    assert "const HPhoto = JAFF_HPhoto.HPhoto" in generated


def test_species_and_bands(generated):
    assert "const NSPEC = 3" in generated
    assert "const NBANDS = 1" in generated
    assert 'tuple("H", "H+", "e-")' in generated
    assert "Honghu.band_edges(::Type{<:Network}) = tuple(13.6, Inf)" in generated


def test_indices_are_one_based(generated):
    """JAFF indexes from zero; every density, photon and flux index is shifted by one."""
    code = [ln for ln in generated.splitlines() if not ln.lstrip().startswith("#")]
    for line in code:
        assert not re.search(r"\b(nden|photden|rflux)\[\d+\]", line), line
    # the state is [n_1..n_3, E, photden_1, rflux_1]: six rows of the right-hand side
    assert re.findall(r"du\[(\d+)\] =", generated) == ["1", "2", "3", "4", "5", "6"]


def test_no_sympy_infinity(generated):
    code = "\n".join(
        ln for ln in generated.splitlines() if not ln.lstrip().startswith("#")
    )
    assert not re.search(r"\boo\b", code)


@pytest.mark.skipif(shutil.which("julia") is None, reason="julia is not installed")
def test_output_parses_as_julia(generated, tmp_path):
    path = tmp_path / "network.jl"
    path.write_text(generated)
    check = (
        "ex = Meta.parseall(read(ARGS[1], String)); "
        "bad = [e for e in ex.args if e isa Expr && e.head in (:error, :incomplete)]; "
        "isempty(bad) || (foreach(println, bad); exit(1))"
    )
    proc = subprocess.run(
        ["julia", "--startup-file=no", "-e", check, str(path)],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
