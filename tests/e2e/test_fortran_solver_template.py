# ABOUTME: The fortran_dlsodes template as a standalone solver: generated for real
# ABOUTME: networks it builds, its Jacobian matches the rhs and DLSODES integrates
"""End-to-end tests of the generated Fortran (``fortran_dlsodes``).

Each network is generated and built with its ``Makefile``. Small test drivers,
compiled against the generated modules, evaluate the rhs and Jacobian or call
``solve`` with chosen inputs. Unknown functions (GOW's tables) are made to return
1, and their partial derivatives 0, by rewriting the generated stubs.
"""

import re
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from jaff import Network
from tests._paths import FIXTURES_DIR, REPO_ROOT
from tests.python_solver import generate

pytestmark = pytest.mark.skipif(
    not sys.platform.startswith("linux")
    or shutil.which("gfortran") is None
    or shutil.which("make") is None,
    reason="needs Linux with gfortran and make",
)

NETWORKS = {
    "demo2": str(REPO_ROOT / "networks" / "demos" / "demo2.jet"),
    "COthin": "COthin",
    "GOW": str(FIXTURES_DIR / "GOW_benchmark" / "GOW.jet"),
}

# Parameter values by name (1 for any other); as in the Python template tests
VALUES = {
    "av": 1.0,
    "chi": 1.0,
    "crate": 1e-16,
    "d2g": 0.01,
    "gradv": 1e-14,
    "tdust": 20.0,
    "vdisp": 1e5,
}
NCOL = 1e14
TEMPERATURE = 1000.0
FFLAGS = ["-ffree-line-length-none", "-O2", "-fallow-argument-mismatch"]


def _stub_unknown_functions(folder: Path) -> None:
    """Make the generated stubs return 1 (functions) and 0 (partial derivatives)."""
    path = folder / "user_functions.f90"
    lines = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.lstrip().startswith("function "):
            value = "0d0" if "_partial_" in line.split("(")[0] else "1d0"
            line = re.sub(r"r = 0d0; error stop \"[^\"]*\"", f"r = {value}", line)
        lines.append(line)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _value(name: str) -> float:
    return NCOL if name.startswith("ncol_") else VALUES.get(name, 1.0)


def _set_parameters(net: Network) -> str:
    """Fortran statements setting every parameter of the chemistry module."""
    names = net.ode_system(thermal="dtdt").parameters
    return "\n".join(f"    {p} = {_value(str(p)):.17e}".replace("e", "d") for p in names)


def _compile_driver(folder: Path, name: str, source: str) -> Path:
    """Compile *source* against the generated modules; return the executable."""
    (folder / f"{name}.f90").write_text(source, encoding="utf-8")
    objects = [str(p) for p in folder.glob("*.o") if p.name != "main.o"]
    exe = folder / name
    subprocess.run(
        ["gfortran", *FFLAGS, f"{name}.f90", *objects, "-o", str(exe)],
        cwd=folder,
        check=True,
        capture_output=True,
        text=True,
    )

    return exe


@pytest.fixture(scope="module", params=list(NETWORKS))
def built(request, tmp_path_factory) -> SimpleNamespace:
    network = NETWORKS[request.param]
    net = Network(network)
    folder = generate(network, tmp_path_factory.mktemp(request.param), "fortran_dlsodes")
    _stub_unknown_functions(folder)
    build = subprocess.run(["make"], cwd=folder, capture_output=True, text=True)

    return SimpleNamespace(net=net, folder=folder, build=build)


def test_writes_the_solver_files(built):
    for name in ("chemistry.f90", "user_functions.f90", "solver.f90", "main.f90"):
        assert (built.folder / name).is_file()
    for name in ("commons.f90", "reactions.f90", "fluxes.f90", "ode.f90"):
        assert not (built.folder / name).exists()


def test_builds(built):
    assert built.build.returncode == 0, built.build.stdout + built.build.stderr


def test_jacobian_matches_finite_differences(built):
    nspec = built.net.species.count
    driver = f"""
program check_jacobian
    use chemistry
    implicit none
    real*8 :: y(nspecs+1), up(nspecs+1), down(nspecs+1)
    real*8 :: jac(nspecs+1, nspecs+1), fd(nspecs+1, nspecs+1), fu(nspecs+1), fdn(nspecs+1)
    real*8 :: h
    integer :: i, j
{_set_parameters(built.net)}
    do i = 1, nspecs
        y(i) = 1d2*(1d0 + 0.37d0*(i - 1))
    end do
    y(nspecs+1) = {TEMPERATURE}d0
    call jacobian(y(1:nspecs), y(nspecs+1), jac)
    do j = 1, nspecs+1
        h = 1d-6*abs(y(j))
        up = y
        down = y
        up(j) = up(j) + h
        down(j) = down(j) - h
        call rhs(up(1:nspecs), up(nspecs+1), fu)
        call rhs(down(1:nspecs), down(nspecs+1), fdn)
        fd(:, j) = (fu - fdn)/(2*h)
    end do
    do i = 1, nspecs+1
        print '(*(es25.16e3))', jac(i, :)
    end do
    do i = 1, nspecs+1
        print '(*(es25.16e3))', fd(i, :)
    end do
end program check_jacobian
"""
    exe = _compile_driver(built.folder, "check_jacobian", driver)
    out = subprocess.run([str(exe)], capture_output=True, text=True, check=True).stdout
    values = np.array([float(v) for v in out.split()]).reshape(2, nspec + 1, nspec + 1)
    jac, fd = values

    for row in range(nspec + 1):
        scale = max(np.abs(fd[row]).max(), np.abs(jac[row]).max(), 1e-300)
        assert jac[row] == pytest.approx(fd[row], rel=1e-5, abs=1e-7 * scale)


def _neutral_start(net: Network) -> str:
    """Every species at 1 cm^-3, electrons equal to the total ion charge."""
    lines = []
    electron = net.species.e_idx()
    charge = sum(s.charge for s in net.species if s.charge > 0)
    for i in range(net.species.count):
        value = float(charge) if i == electron else 1.0
        lines.append(f"    n({i + 1}) = {value:.17e}".replace("e", "d"))

    return "\n".join(lines)


def _has_unknown_functions(net: Network) -> bool:
    return bool(net.symbols.applied_functions)


@pytest.mark.parametrize("evolve", [False, True])
def test_solve_integrates_with_dlsodes(built, evolve):
    if evolve and _has_unknown_functions(built.net):
        # Stubbed tables (= 1) make GOW cool at ~1e14 K/s: T -> 0 within 1e-11 s.
        pytest.skip("temperature evolution needs the real tables (benchmarks)")
    flag = ".true." if evolve else ".false."
    driver = f"""
program run_solve
    use chemistry
    use solver
    implicit none
    real*8 :: n(nspecs), temperature
{_set_parameters(built.net)}
{_neutral_start(built.net)}
    temperature = {TEMPERATURE}d0
    evolve_temperature = {flag}
    call solve(n, temperature, 1d8, "run.dat")
    print '(*(es25.16e3))', temperature, n
end program run_solve
"""
    exe = _compile_driver(built.folder, "run_solve", driver)
    proc = subprocess.run([str(exe)], cwd=built.folder, capture_output=True, text=True)

    assert proc.returncode == 0, proc.stdout + proc.stderr
    final = np.array([float(v) for v in proc.stdout.split()])
    assert np.all(np.isfinite(final))
    if not evolve:
        assert final[0] == pytest.approx(TEMPERATURE)
    lines = (built.folder / "run.dat").read_text().splitlines()
    names = built.net.species.names()
    assert lines[0].split() == ["#", "t[s]", "T[K]", *names]
    assert len(lines) == 1 + 101


def test_unset_parameters_are_named(built):
    names = [str(p) for p in built.net.ode_system(thermal="dtdt").parameters]
    driver = """
program unset_check
    use chemistry
    use solver
    implicit none
    real*8 :: n(nspecs), temperature
    n = 1d0
    temperature = 1d2
    call solve(n, temperature, 1d0, "unset.dat")
end program unset_check
"""
    exe = _compile_driver(built.folder, "unset", driver)
    proc = subprocess.run([str(exe)], cwd=built.folder, capture_output=True, text=True)

    assert proc.returncode != 0
    for name in names:
        assert f"parameter not set: {name}" in proc.stdout


def test_generated_driver_runs(built):
    if _has_unknown_functions(built.net):
        pytest.skip("the generated driver needs the real table functions")
    proc = subprocess.run(["./test"], cwd=built.folder, capture_output=True, text=True)

    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "evolution saved to solution.dat" in proc.stdout
    assert len((built.folder / "solution.dat").read_text().splitlines()) == 1 + 101
