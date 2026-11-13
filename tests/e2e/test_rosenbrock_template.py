# ABOUTME: The rosenbrock_cpp template: rhs_jac() (rhs and Jacobian sharing one CSE group)
# ABOUTME: matches rhs() and finite differences of it, and the integrator builds and runs

import re
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

from tests._paths import FIXTURES_DIR
from tests.jaffgen_run import run_jaffgen

NETWORK = FIXTURES_DIR / "react_cie_hepp.jet"
GXX = shutil.which("g++")
needs_gxx = pytest.mark.skipif(GXX is None, reason="g++ is not installed")

HARNESS = r"""
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "integrator_data.hpp"
#include "parameters.hpp"
#include "types.hpp"
#include "rosenbrock.hpp"
#include "jac.hpp"
#include "rhs.hpp"
#include "state.hpp"
#include "eos.hpp"

static void rhs_at(const std::vector<double> &y, DArray1D &out) {
  state<double> s;
  for (int i = 0; i < params::NumSpec; ++i) s.rho[i] = y[i];
  s.e = y[params::NumSpec];
  rhs(s, out);
}

int main() {
  constexpr int n = INT_NEQS;
  state<double> s;
  for (int i = 0; i < params::NumSpec; ++i) s.rho[i] = 1.0e2 * (1.0 + 0.37 * i);
  s.e = eos::energy_from_temperature(800.0, s.n_tot());

  DArray1D f(n), g(n);
  DArray2D J(n, DArray1D(n, 0.0));
  rhs(s, f);
  rhs_jac(s, g, J);

  double scale = 0.0;
  for (int i = 0; i < n; ++i) scale = std::max(scale, std::abs(f[i]));
  double rhs_err = 0.0;
  for (int i = 0; i < n; ++i) rhs_err = std::max(rhs_err, std::abs(g[i] - f[i]) / scale);

  std::vector<double> y(n);
  for (int i = 0; i < params::NumSpec; ++i) y[i] = s.rho[i];
  y[params::NumSpec] = s.e;
  std::vector<std::vector<double>> fd(n, std::vector<double>(n));
  DArray1D up(n), down(n);
  for (int j = 0; j < n; ++j) {
    double h = 1.0e-6 * std::abs(y[j]);
    std::vector<double> yu = y, yd = y;
    yu[j] += h;
    yd[j] -= h;
    rhs_at(yu, up);
    rhs_at(yd, down);
    for (int i = 0; i < n; ++i) fd[i][j] = (up[i] - down[i]) / (2.0 * h);
  }

  // Each entry's contribution J_ij * y_j against the row's largest one.
  double jac_err = 0.0;
  for (int i = 0; i < n; ++i) {
    double row = 0.0;
    for (int j = 0; j < n; ++j)
      row = std::max(row, std::max(std::abs(fd[i][j]), std::abs(J[i][j])) * std::abs(y[j]));
    if (row == 0.0) continue;
    for (int j = 0; j < n; ++j)
      jac_err = std::max(jac_err, std::abs(J[i][j] - fd[i][j]) * std::abs(y[j]) / row);
  }

  std::printf("rhs %.6e\njac %.6e\n", rhs_err, jac_err);
  return 0;
}
"""


@pytest.fixture(scope="module")
def rendered(tmp_path_factory) -> Path:
    out = tmp_path_factory.mktemp("rosenbrock")
    args = SimpleNamespace(
        network=str(NETWORK),
        config=None,
        label=None,
        funcfile=None,
        duplicate_policy=None,
        expand_nuclei=None,
        strict=None,
        network_config=None,
        outdir=str(out),
        indir=None,
        files=None,
        template="rosenbrock_cpp",
    )

    return run_jaffgen(args).output_dir


def _compile(source: Path, include: Path, binary: Path) -> None:
    subprocess.run(
        [GXX, "-std=c++17", "-O1", f"-I{include}", "-o", str(binary), str(source)],
        check=True,
        capture_output=True,
        text=True,
    )


def test_rhs_jac_replaces_jac(rendered):
    jac = (rendered / "jac.hpp").read_text()

    assert "inline void rhs_jac(" in jac
    assert not re.search(r"inline void jac\(", jac)


def test_evaluate_jacobian_calls_rhs_jac(rendered):
    rosenbrock = (rendered / "rosenbrock.hpp").read_text()

    assert "rhs_jac(state, rstate.ak1, rstate.jac);" in rosenbrock


def test_rhs_jac_shares_one_set_of_temporaries(rendered):
    jac = (rendered / "jac.hpp").read_text()
    body = jac[jac.index("inline void rhs_jac(") :]

    declared = set(re.findall(r"const double (cse\d+) =", body))
    read = set(re.findall(r"\bcse\d+\b", body))
    assert declared and read <= declared


@needs_gxx
def test_rhs_jac_matches_rhs_and_finite_differences(rendered, tmp_path):
    source = tmp_path / "harness.cpp"
    source.write_text(HARNESS)
    binary = tmp_path / "harness"
    _compile(source, rendered, binary)
    run = subprocess.run([str(binary)], check=True, capture_output=True, text=True)
    errors = dict(line.split() for line in run.stdout.splitlines())

    assert float(errors["rhs"]) == pytest.approx(0.0, abs=1e-12)
    assert float(errors["jac"]) == pytest.approx(0.0, abs=1e-5)


@needs_gxx
def test_integrator_builds_and_runs(rendered, tmp_path):
    binary = tmp_path / "rosenbrock"
    _compile(rendered / "testing" / "main.cpp", rendered, binary)
    run = subprocess.run(
        [str(binary)], check=True, capture_output=True, text=True, cwd=tmp_path
    )

    assert "Integration successful" in run.stdout
