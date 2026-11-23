# ABOUTME: Every Codegen *_str method passes cse_suffix on to its get_indexed_* getter,
# ABOUTME: so its CSE temporaries are named <cse_var><n><suffix> as in templates

import re

import pytest

from jaff import Network
from jaff.codegen import Codegen


@pytest.fixture(scope="module")
def cg():
    net = Network("h_photoionization", radiation={"bands": [1, 13.6, 100, "inf"]})

    return Codegen(net, lang="cxx")


@pytest.fixture(scope="module")
def gow():
    # h_photoionization's rates share no subexpression, so no rate temporaries.
    return Codegen(Network("GOW"), lang="cxx")


@pytest.mark.parametrize(
    ("method", "cse_var", "options"),
    [
        ("get_ode_str", "cse", {}),
        ("get_rhs_str", "cse", {"thermal": "dedt"}),
        ("get_radode_str", "rcse", {}),
        ("get_jacobian_str", "cse", {"thermal": "dtdt"}),
    ],
)
def test_temporaries_carry_the_suffix(cg, method, cse_var, options):
    code = getattr(cg, method)(use_cse=True, cse_suffix="_s", **options)
    names = set(re.findall(rf"\b{cse_var}\d+\w*", code))

    assert names, f"{method} produced no CSE temporaries"
    assert all(name.endswith("_s") for name in names), sorted(names)


def test_rate_temporaries_carry_the_suffix(gow):
    code = gow.get_rates_str(use_cse=True, cse_suffix="_s")
    names = set(re.findall(r"\bx\d+\w*", code))

    assert names
    assert all(name.endswith("_s") for name in names), sorted(names)
