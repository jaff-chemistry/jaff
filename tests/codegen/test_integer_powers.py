# ABOUTME: Codegen prints small integer powers as products (no std::pow(x, 3) left) and
# ABOUTME: the printed Python keeps its values

import math
import re

import pytest

from jaff.codegen import Codegen

NETWORK = [
    "H + H -> H2 [10,1e4] 1e-10*(tgas/300)**3",
    "H2 + H -> H + H + H [10,1e4] 2e-10*(tgas/300)**(-2)*exp(-5e4/tgas)",
]
INTEGER_POW = re.compile(r"pow\([^()]*, -?\d+\)")


@pytest.fixture
def network(make_network):
    return make_network(NETWORK)


def test_cxx_has_no_integer_pow(network):
    group = Codegen(network, "c++").cse_group([("rhs", {}), ("jacobian", {})])
    code = [c for _, c in group.temporaries] + [
        c for part in group.expressions for _, c in part
    ]

    assert not [c for c in code if INTEGER_POW.search(c)]


@pytest.mark.parametrize("tgas", [50.0, 900.0])
def test_python_rates_keep_their_values(network, tgas):
    rates = Codegen(network, "python").get_indexed_rates(use_cse=False)
    namespace = {"math": math, "tgas": tgas}
    printed = [eval(code, namespace) for _, code in rates["expressions"]]
    want = [
        1e-10 * (tgas / 300) ** 3,
        2e-10 * (tgas / 300) ** -2 * math.exp(-5e4 / tgas),
    ]

    assert printed == pytest.approx(want, rel=1e-12)
