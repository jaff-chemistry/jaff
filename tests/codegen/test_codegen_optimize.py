# ABOUTME: Codegen's optimizations: no integer, numeric-base or shared-base std::pow left
# ABOUTME: in C++ (one log(tgas) instead), and the printed Python keeps its values

import math
import re

import pytest
import sympy as sp

from jaff.codegen import Codegen

NETWORK = [
    "H + H -> H2 [10,1e4] 1e-10*(tgas/300)**3",
    "H2 + H -> H + H + H [10,1e4] 2e-10*(tgas/300)**(-2)*exp(-5e4/tgas)",
    "H2 + H2 -> H2 + H + H [10,1e4] 10**(-12 - 1.5*log10(tgas))",
]
INTEGER_POW = re.compile(r"pow\([^()]*, -?\d+\)")
NUMERIC_BASE_POW = re.compile(r"pow\(\d")


@pytest.fixture
def network(make_network):
    return make_network(NETWORK)


def test_cxx_has_no_integer_pow(network):
    group = Codegen(network, "c++").cse_group([("rhs", {}), ("jacobian", {})])
    code = [c for _, c in group.temporaries] + [
        c for part in group.expressions for _, c in part
    ]

    assert not [c for c in code if INTEGER_POW.search(c)]
    assert not [c for c in code if NUMERIC_BASE_POW.search(c)]


@pytest.mark.parametrize("tgas", [50.0, 900.0])
def test_python_rates_keep_their_values(network, tgas):
    rates = Codegen(network, "python").get_indexed_rates(use_cse=False)
    namespace = {"math": math, "tgas": tgas}
    printed = [eval(code, namespace) for _, code in rates["expressions"]]
    want = [
        1e-10 * (tgas / 300) ** 3,
        2e-10 * (tgas / 300) ** -2 * math.exp(-5e4 / tgas),
        10 ** (-12 - 1.5 * math.log10(tgas)),
    ]

    assert printed == pytest.approx(want, rel=1e-12)


POWER_NETWORK = [
    "H + H -> H2 [10,1e4] 1e-10*tgas**0.7",
    "H2 + H -> H + H + H [10,1e4] 2e-10*tgas**(-1.3)*exp(-5e4/tgas)",
    "H2 + H2 -> H2 + H + H [10,1e4] 3e-11*tgas**2.7",
]
DENSITIES = [1.0e2, 3.0e1]


@pytest.fixture
def power_network(make_network):
    return make_network(POWER_NETWORK)


def test_powers_of_tgas_share_one_log(power_network):
    group = Codegen(power_network, "c++").cse_group([("rhs", {}), ("jacobian", {})])
    code = " ".join(
        [c for _, c in group.temporaries]
        + [c for part in group.expressions for _, c in part]
    )

    assert "std::pow(tgas" not in code
    assert code.count("std::log(tgas)") == 1


@pytest.mark.parametrize("tgas", [50.0, 900.0])
def test_python_rhs_and_jacobian_keep_their_values(power_network, tgas):
    group = Codegen(power_network, "python").cse_group([("rhs", {}), ("jacobian", {})])
    namespace = {"math": math, "tgas": tgas, "nden": DENSITIES}
    for idx, code in group.temporaries:
        exec(f"cse{idx[0]} = {code}", namespace)
    printed = [eval(code, namespace) for part in group.expressions for _, code in part]

    system = power_network.ode_system()
    entries = [e for _, e in sorted(system.jacobian.todok().items())]
    want = []
    for expr in [*system.rhs.flat(), *entries]:
        subs = {i: DENSITIES[int(i.indices[0])] for i in expr.atoms(sp.Indexed)}
        want.append(float(expr.xreplace({**subs, sp.Symbol("tgas"): tgas})))

    assert printed == pytest.approx(want, rel=1e-13)
