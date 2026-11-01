# ABOUTME: get_indexed_rhs_jacobian agrees numerically with get_indexed_rhs and get_indexed_jacobian
# ABOUTME: Both are rendered as Python and evaluated at random states of real networks

import math
import random
from pathlib import Path

import pytest

from jaff import Network
from jaff.codegen import Codegen

REPO = Path(__file__).resolve().parents[1]
NETWORKS = [
    REPO / "networks" / "cie_h" / "react_cie_h.jet",
    REPO / "networks" / "demos" / "demo2.jet",
    Path(__file__).parent / "fixtures" / "test_jac_dedt.dat",
]


def _evaluate(ir_parts, nden, tgas):
    """Run CSE temporaries, then return {index tuple: value} for the expressions."""
    # av and crate are the extinction and cosmic-ray rate that some networks read.
    scope = {"math": math, "nden": nden, "tgas": tgas, "av": 2.0, "crate": 1.3e-17}
    out = {}
    for ir in ir_parts:
        for [idx], expr in ir["extras"]["cse"]:
            scope[f"cse{idx}"] = eval(expr, scope)
        for idx, expr in ir["expressions"]:
            out[tuple(idx)] = eval(expr, scope)
    return out


def _close(a, b, scale):
    return abs(a - b) <= 1e-9 * scale + 1e-300


@pytest.mark.parametrize("path", NETWORKS, ids=lambda p: p.name)
def test_shared_rhs_jacobian_matches_separate_getters(path):
    net = Network(str(path))
    cg = Codegen(net, "python")
    old_rhs = cg.get_indexed_rhs(thermal="dedt")
    old_jac = cg.get_indexed_jacobian(thermal="dedt")
    new_rhs = cg.get_indexed_rhs_jacobian(part="rhs", thermal="dedt")
    new_jac = cg.get_indexed_rhs_jacobian(part="jacobian", thermal="dedt")
    assert not new_jac["extras"]["cse"]

    rng = random.Random(3)
    n_species = net.species.count
    for _ in range(25):
        nden = [10 ** rng.uniform(-6, 8) for _ in range(n_species)]
        tgas = 10 ** rng.uniform(1.3, 4.0)

        f_old = _evaluate([old_rhs], nden, tgas)
        f_new = _evaluate([new_rhs], nden, tgas)
        assert f_old.keys() == f_new.keys()
        for key in f_old:
            assert _close(f_new[key], f_old[key], abs(f_old[key])), key

        j_old = _evaluate([old_jac], nden, tgas)
        j_new = _evaluate([new_rhs, new_jac], nden, tgas)
        j_new = {k: v for k, v in j_new.items() if len(k) == 2}
        rows = {i for i, _ in j_old} | {i for i, _ in j_new}
        for i in rows:
            scale = max(abs(v) for (r, _), v in j_old.items() if r == i)
            for j in range(n_species + 1):
                assert _close(j_new.get((i, j), 0.0), j_old.get((i, j), 0.0), scale), (
                    i,
                    j,
                )
