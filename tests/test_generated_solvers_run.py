# ABOUTME: Generated solvers actually run: jaffgen output for each template is built
# ABOUTME: (where needed) and executed, and must exit 0 without solver errors

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
DEMO2 = REPO / "networks" / "demos" / "demo2.jet"


def _jaffgen(network: Path, template: str, outdir: Path) -> None:
    """Run the jaffgen CLI the way a user would and require success."""
    # The console script installed next to this interpreter (jaffgen.exe on Windows).
    jaffgen = shutil.which("jaffgen", path=str(Path(sys.executable).parent))
    assert jaffgen is not None, "jaffgen console script not found"

    proc = subprocess.run(
        [
            jaffgen,
            "--network",
            str(network),
            "--template",
            template,
            "--outdir",
            str(outdir),
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr


def test_python_solve_ivp_driver_runs(tmp_path):
    _jaffgen(DEMO2, "python_solve_ivp", tmp_path)

    proc = subprocess.run(
        [sys.executable, "main.py"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        env={**os.environ, "MPLBACKEND": "Agg"},
    )

    assert proc.returncode == 0, proc.stdout + proc.stderr


@pytest.mark.skipif(
    not sys.platform.startswith("linux")
    or shutil.which("gfortran") is None
    or shutil.which("make") is None,
    reason="needs Linux with gfortran and make",
)
def test_fortran_dlsodes_builds_and_runs(tmp_path):
    _jaffgen(DEMO2, "fortran_dlsodes", tmp_path)

    build = subprocess.run(["make"], cwd=tmp_path, capture_output=True, text=True)
    assert build.returncode == 0, build.stdout + build.stderr

    run = subprocess.run(["./test"], cwd=tmp_path, capture_output=True, text=True)
    output = run.stdout + run.stderr

    assert run.returncode == 0, output
    assert "DLSODES-" not in output, output

    # The driver holds the gas temperature fixed at its initial 100 K.
    tgas = re.search(r"tgas\s*=\s*(\S+)", output)
    assert tgas, output
    assert float(tgas.group(1).replace("D", "E")) == pytest.approx(100.0)


@pytest.mark.skipif(
    not sys.platform.startswith("linux")
    or shutil.which("cmake") is None
    or shutil.which("clang++") is None,
    reason="needs Linux with cmake and clang++ (and network access for FetchContent)",
)
def test_kokkos_ode_builds_and_runs(tmp_path):
    _jaffgen(DEMO2, "kokkos_ode", tmp_path)

    for cmd in (
        ["cmake", "-S", ".", "-B", "build", "-DCMAKE_CXX_COMPILER=clang++"],
        ["cmake", "--build", "build", "-j"],
        ["./build/chemistry_ode"],
    ):
        proc = subprocess.run(cmd, cwd=tmp_path, capture_output=True, text=True)
        assert proc.returncode == 0, f"{cmd}:\n{proc.stdout}{proc.stderr}"

    assert "Number of reactions: 2" in proc.stdout, proc.stdout
    # The driver reports element and charge conservation from injected metadata.
    assert "net_charge:" in proc.stdout, proc.stdout
