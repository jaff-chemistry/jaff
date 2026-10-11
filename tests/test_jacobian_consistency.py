# ABOUTME: Checks generated analytic Jacobians against central finite differences of the
# ABOUTME: generated RHS, for every internal-energy normalisation, on random inputs

# The state is (n_0 .. n_{N-1}, e[, radiation]) where e is the evolved internal
# energy of each variant.  Species columns are taken at constant e, so they
# exercise the EOS chain rule; the e column exercises the de/dT normalisation.
# Radiation columns are not checked (their state ordering is not exposed to the
# rendered module); radiation rows are.

from pathlib import Path
from typing import Callable, Dict

import numpy as np
import pytest

from jaff.physics.constants import k_B
from jaff.physics import EosProps
from tests.codegen_render import NETWORKS, draw_inputs, evaluate, free_names, load

KB = k_B.cgs.value
GAMMA = EosProps("ideal").gamma
FD_STEP = 1e-6
TOL = 1e-5


def _volumetric(n: np.ndarray, T: float, m: np.ndarray) -> float:
    return n.sum() * KB * T / (GAMMA - 1.0)


ENERGY: Dict[str, Callable[[np.ndarray, float, np.ndarray], float]] = {
    "volumetric": _volumetric,
    "specific_mass": lambda n, T, m: _volumetric(n, T, m) / (m @ n),
    "specific_number": lambda n, T, m: _volumetric(n, T, m) / n.sum(),
}

VARIANTS = list(ENERGY)


@pytest.mark.parametrize("variant", VARIANTS)
@pytest.mark.parametrize("name", list(NETWORKS))
def test_jacobian_matches_finite_differences(
    name: str, variant: str, rendered: Callable[[str], Path], render_seed: int
) -> None:
    path = rendered(name) / "codegen_outputs.py"
    module = load(path, f"fd_{name}_{variant}")
    inputs = draw_inputs(free_names(path.read_text()), render_seed)
    energy = ENERGY[variant]

    masses = np.array([module.SPECIE_MASSES[i] for i in range(len(inputs["nden"]))])
    n0 = np.array(inputs["nden"])
    T0 = inputs["tgas"]
    e0 = energy(n0, T0, masses)

    def rhs(n: np.ndarray, e: float) -> np.ndarray:
        T = T0 * e / energy(n, T0, masses)  # every variant is linear in T
        state = {**inputs, "nden": tuple(float(v) for v in n), "tgas": T}
        out = evaluate(module, f"rhs_{variant}", state, render_seed)
        return np.array([out[i] for i in range(len(out))])

    jac = evaluate(module, f"jacobian_{variant}", inputs, render_seed)
    n_rows = len(rhs(n0, e0))
    n_cols = len(n0) + 1
    y0 = np.append(n0, e0)

    fd = np.zeros((n_rows, n_cols))
    for col in range(n_cols):
        h = FD_STEP * y0[col]
        up, down = y0.copy(), y0.copy()
        up[col] += h
        down[col] -= h
        fd[:, col] = (rhs(up[:-1], up[-1]) - rhs(down[:-1], down[-1])) / (2.0 * h)

    analytic = np.array(
        [[jac.get((row, col), 0.0) for col in range(n_cols)] for row in range(n_rows)]
    )

    # Compare each entry's contribution J_ij * y_j against the row's largest one.
    scale = np.maximum(np.abs(fd), np.abs(analytic)) * np.abs(y0)
    row_scale = scale.max(axis=1, keepdims=True)
    row_scale[row_scale == 0.0] = 1.0
    err = np.abs(analytic - fd) * np.abs(y0) / row_scale

    bad = np.argwhere(err > TOL)
    rows = "\n".join(
        f"  J[{i}, {j}]: analytic={analytic[i, j]:.6e} fd={fd[i, j]:.6e} "
        f"err={err[i, j]:.1e}"
        for i, j in bad[:20]
    )
    assert not bad.size, (
        f"{name}/{variant} (seed={render_seed}): {len(bad)} Jacobian entries disagree "
        f"with finite differences (tol={TOL})\n{rows}"
    )
