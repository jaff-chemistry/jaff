# ABOUTME: CSE temporaries read only inside a rate's temperature range are guarded: the
# ABOUTME: printed rhs and Jacobian run outside that range and match the symbolic values

import math

import pytest
import sympy as sp

from jaff.codegen import Codegen

# exp(2e6/T) is valid above 5000 K; below ~2700 K it overflows a double.
NETWORK = [
    "H + H -> H2 [5e3,1e8] 1e-10*exp(2e6/tgas)",
    "H2 + H -> H + H + H [5e3,1e8] 2e-10*exp(2e6/tgas)",
]
DENSITIES = [1.0e2, 3.0e1]


@pytest.fixture
def network(make_network):
    return make_network(NETWORK)


def _run(temporaries, expressions, tgas):
    """The printed Python of *expressions* evaluated at *tgas*."""
    namespace = {"math": math, "tgas": tgas, "nden": DENSITIES}
    for idx, code in temporaries:
        exec(f"cse{idx[0]} = {code}", namespace)

    return [eval(code, namespace) for _, code in expressions]


def _symbolic(network, tgas):
    """The rhs and the non-zero Jacobian entries of the ODE system at *tgas*."""
    system = network.ode_system()
    entries = [e for _, e in sorted(system.jacobian.todok().items())]
    values = []
    for expr in [*system.rhs.flat(), *entries]:
        subs = {i: DENSITIES[int(i.indices[0])] for i in expr.atoms(sp.Indexed)}
        values.append(float(expr.xreplace({**subs, sp.Symbol("tgas"): tgas})))

    return values


def test_branch_only_temporary_is_guarded(network):
    group = Codegen(network, "python").cse_group([("rhs", {}), ("jacobian", {})])
    exp_lines = [code for _, code in group.temporaries if "math.exp" in code]

    assert exp_lines
    assert all(code.endswith("else (0.0))") for code in exp_lines)


@pytest.mark.parametrize("tgas", [1000.0, 6000.0, 2.0e4])
def test_rhs_and_jacobian_run_outside_the_rate_range(network, tgas):
    group = Codegen(network, "python").cse_group([("rhs", {}), ("jacobian", {})])
    printed = [
        v for part in group.expressions for v in _run(group.temporaries, part, tgas)
    ]

    assert printed == pytest.approx(_symbolic(network, tgas), rel=1e-12)
