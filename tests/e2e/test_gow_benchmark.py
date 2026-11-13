# ABOUTME: Benchmark: the generated Python solver evolves the GOW network to equilibrium
# ABOUTME: at n_H = 100 and chi = 0.1, 1, 10 and matches the reference within 1e-6
"""GOW one-zone equilibrium against reference runs.

The network is ``tests/fixtures/GOW_benchmark``: ``networks/GOW`` without the
He + CR heating term (``deltaE2``), as in the reference runs. The references
(``reference/state_over_time_n100_00_chi*_met1_00.txt``, Microphysics runs of the
same network) start from the Gong+17 abundances at n_H = 100 cm^-3 and 100 K and
reach equilibrium by ~2e13 s; their last row (7e20 s) is the equilibrium state.

The RHS and Jacobian are generated from the network by the current code; the
test supplies only the table functions (:mod:`tests.gow_tables`), the start state
(the reference's first row) and the parameters below. The gas sits at the cloud
face: no extinction, all column densities 0, so every shielding factor is 1.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from tests import gow_tables
from tests._paths import FIXTURES_DIR
from tests.python_solver import generate, load

BENCHMARK = FIXTURES_DIR / "GOW_benchmark"

#: Reference file name tag -> FUV field chi (Draine units)
CHI = {"0_10": 0.1, "1_00": 1.0, "10_00": 10.0}

#: Reference column -> species name in the network
SPECIES = {
    "H": "H", "Hp": "H+", "E": "e-", "H2": "H2", "H2p": "H2+", "HE": "He",
    "HEp": "He+", "Carbon": "C", "Cp": "C+", "CO": "CO", "HCOp": "HCO+", "O": "O",
    "SI": "Si", "SIp": "Si+", "CH": "CHx", "OH": "OHx", "H3p": "H3+", "Op": "O+",
}  # fmt: skip

#: End of the integration [s]: 5x the slowest approach to equilibrium (chi = 10)
T_END = 1.0e14

RTOL = 1e-12
ATOL = 1e-20
#: Required agreement with the reference equilibrium (relative)
AGREEMENT = 1e-6


def _parameters(chi: float) -> dict[str, float]:
    return {
        "av": 0.0,
        "chi": chi,
        "crate": 2e-16,
        "d2g": 0.01,
        "gradv": 3.0 * 3.0e-14,
        "vdisp": 3.0 * 1.0e5 / 2**0.5,
        "ncol_C": 0.0,
        "ncol_CO": 0.0,
        "ncol_H": 0.0,
        "ncol_H2": 0.0,
    }


def _reference(tag: str) -> dict[str, np.ndarray]:
    """The reference run for *tag* as ``{column: values over time}``."""
    path = BENCHMARK / "reference" / f"state_over_time_n100_00_chi{tag}_met1_00.txt"
    header = path.read_text(encoding="utf-8").splitlines()[0].replace("#", "").split()
    data = np.loadtxt(path)

    return {name: data[:, i] for i, name in enumerate(header)}


@pytest.fixture(scope="module")
def solver(tmp_path_factory) -> SimpleNamespace:
    folder = generate(str(BENCHMARK / "GOW.jet"), tmp_path_factory.mktemp("gow"))
    modules = load(folder)
    gow_tables.install(modules.chemistry, folder)

    return modules


@pytest.mark.parametrize("tag", list(CHI))
def test_equilibrium_matches_the_reference(solver, tag):
    ref = _reference(tag)
    n0 = {name: ref[column][0] for column, name in SPECIES.items()}

    result = solver.solver.solve(
        n0,
        ref["Tgas"][0],
        T_END,
        _parameters(CHI[tag]),
        evolve_temperature=True,
        t_eval=[T_END],
        rtol=RTOL,
        atol=ATOL,
    )

    assert result.temperature[-1] == pytest.approx(ref["Tgas"][-1], rel=AGREEMENT)
    for column, name in SPECIES.items():
        assert result.n[name][-1] == pytest.approx(ref[column][-1], rel=AGREEMENT), name
