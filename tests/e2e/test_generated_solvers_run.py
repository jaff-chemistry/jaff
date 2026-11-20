# ABOUTME: Generated solvers actually run: jaffgen output for each template is built
# ABOUTME: (where needed) and executed, and must exit 0 without solver errors

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from tests._paths import REPO_ROOT
from tests.jaffgen_run import run_folder

REPO = REPO_ROOT
DEMO2 = REPO / "networks" / "demos" / "demo2.jet"


def _jaffgen(network: Path, template: str, outdir: Path) -> Path:
    """Run the jaffgen CLI the way a user would, require success; return the run folder."""
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

    return run_folder(outdir)


def test_python_solve_ivp_driver_runs(tmp_path):
    out = _jaffgen(DEMO2, "python_solve_ivp", tmp_path)

    proc = subprocess.run(
        [sys.executable, "main.py"],
        cwd=out,
        capture_output=True,
        text=True,
        env={**os.environ, "MPLBACKEND": "Agg"},
    )

    assert proc.returncode == 0, proc.stdout + proc.stderr

    # The driver saves time, temperature and every species' density per output time.
    lines = (out / "solution.dat").read_text().splitlines()
    assert lines[0].split() == ["#", "t[s]", "T[K]", "e-", "CO+", "CO"]
    assert len(lines) == 1 + 101
    assert "evolution saved to solution.dat" in proc.stdout


@pytest.mark.skipif(
    not sys.platform.startswith("linux")
    or shutil.which("gfortran") is None
    or shutil.which("make") is None,
    reason="needs Linux with gfortran and make",
)
def test_fortran_dlsodes_builds_and_runs(tmp_path):
    out = _jaffgen(DEMO2, "fortran_dlsodes", tmp_path)

    build = subprocess.run(["make"], cwd=out, capture_output=True, text=True)
    assert build.returncode == 0, build.stdout + build.stderr

    run = subprocess.run(["./test"], cwd=out, capture_output=True, text=True)
    output = run.stdout + run.stderr

    assert run.returncode == 0, output
    assert "DLSODES-" not in output, output

    # The driver holds the gas temperature fixed at its initial 100 K.
    tgas = re.search(r"T\s*=\s*(\S+)\s*K", output)
    assert tgas, output
    assert float(tgas.group(1)) == pytest.approx(100.0)
    assert "evolution saved to solution.dat" in output


@pytest.mark.skipif(
    not sys.platform.startswith("linux")
    or shutil.which("cmake") is None
    or shutil.which("clang++") is None,
    reason="needs Linux with cmake and clang++ (and network access for FetchContent)",
)
def test_kokkos_ode_builds_and_runs(tmp_path):
    out = _jaffgen(DEMO2, "kokkos_ode", tmp_path)

    for cmd in (
        ["cmake", "-S", ".", "-B", "build", "-DCMAKE_CXX_COMPILER=clang++"],
        ["cmake", "--build", "build", "-j"],
        ["./build/chemistry_ode"],
    ):
        proc = subprocess.run(cmd, cwd=out, capture_output=True, text=True)
        assert proc.returncode == 0, f"{cmd}:\n{proc.stdout}{proc.stderr}"

    # The driver reports element and charge conservation from injected metadata.
    assert "net_charge:" in proc.stdout, proc.stdout
    assert "evolution saved to solution.dat" in proc.stdout, proc.stdout
