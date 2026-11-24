# ABOUTME: Tests for Codegen.cse_group: one CSE over several expression kinds, run as
# ABOUTME: Python and compared with the separate getters; one member prints as its getter

import math
import random
import re

import pytest

from jaff import Network
from jaff.codegen import Codegen

BANDS = {"bands": [1, 13.6, 100, "inf"]}
TGAS_NETWORK = [
    "H+ + e- -> H [10,1e4] 2e-10*tgas**(-0.7)",
    "H + H -> H2 [10,1e4] 1e-17*(tgas/300)**0.5",
    "H2 + e- -> H + H + e- [10,1e4] 3e-9*exp(-5e4/tgas)",
    "H2 + H -> H + H + H [10,1e4] 3e-9*exp(-5e4/tgas)",
]
GETTERS = {
    "rates": "get_indexed_rates",
    "odes": "get_indexed_odes",
    "radodes": "get_indexed_radodes",
    "rhs": "get_indexed_rhs",
    "jacobian": "get_indexed_jacobian",
}

PAIRS = [
    pytest.param("tgas", {}, id="none"),
    pytest.param("tgas", {"thermal": "dtdt"}, id="dtdt"),
    *(
        pytest.param("tgas", {"thermal": "dedt", "energy": e}, id=f"dedt-{e}")
        for e in ("volumetric", "specific", "per_particle", "molar")
    ),
    pytest.param("photo", {"thermal": "dtdt"}, id="photorates"),
    pytest.param("rad", {"thermal": "dtdt", "radiation": True, "rad_order": 2}, id="rad"),
    # dE/dt with radiation, as a solver evolving specific energy and photon bands uses.
    *(
        pytest.param(
            "rad",
            {"thermal": "dedt", "energy": e, "radiation": True, "rad_order": 2},
            id=f"rad-dedt-{e}",
        )
        for e in ("volumetric", "specific")
    ),
]


@pytest.fixture(scope="module")
def networks(tmp_path_factory):
    path = tmp_path_factory.mktemp("cse_group") / "tgas.dat"
    path.write_text("\n".join(TGAS_NETWORK) + "\n", encoding="utf-8")

    return {
        "tgas": Network(str(path), funcfile=False),
        "photo": Network("h_photoionization"),
        "rad": Network("h_photoionization", radiation=BANDS),
    }


def _state(seed: int) -> dict:
    rng = random.Random(seed)
    vector = [rng.uniform(1.0, 10.0) for _ in range(16)]

    return {
        "math": math,
        "tgas": rng.uniform(100.0, 3000.0),
        "nden": vector,
        "photden": vector[::-1],
        "radeden": vector[::-1],
        "rflux": [0.5 * v for v in vector],
        "photorates": lambda *args: 1e-10 * args[0],
    }


def _define(temporaries, namespace: dict, cse_var: str = "cse", suffix: str = "") -> None:
    for (idx,), code in temporaries:
        namespace[f"{cse_var}{idx}{suffix}"] = eval(code, namespace)


def _values(expressions, namespace: dict) -> dict:
    return {tuple(idx): eval(code, namespace) for idx, code in expressions}


def _getter(cg, kind, options, seed) -> dict:
    printed = getattr(cg, GETTERS[kind])(**options, cse_var="cse")
    namespace = _state(seed)
    _define(printed["extras"]["cse"], namespace)

    return _values(printed["expressions"], namespace)


def _close(got: dict, want: dict) -> None:
    assert set(got) == set(want)
    scale = max(abs(v) for v in want.values())
    for key, value in want.items():
        assert got[key] == pytest.approx(value, rel=1e-9, abs=1e-12 * scale), key


@pytest.mark.parametrize(("which", "options"), PAIRS)
@pytest.mark.parametrize("seed", [0, 1])
def test_rhs_and_jacobian_match_the_separate_getters(networks, which, options, seed):
    cg = Codegen(networks[which], lang="python")
    group = cg.cse_group([("rhs", options), ("jacobian", options)])
    namespace = _state(seed)
    _define(group.temporaries, namespace)
    rhs, jac = (_values(e, namespace) for e in group.expressions)

    _close(rhs, _getter(cg, "rhs", options, seed))
    _close(jac, _getter(cg, "jacobian", options, seed))


def test_odes_and_radodes_match_the_separate_getters(networks):
    cg = Codegen(networks["rad"], lang="python")
    group = cg.cse_group([("odes", {}), ("radodes", {"order": 1})])
    namespace = _state(3)
    _define(group.temporaries, namespace)
    odes, radodes = (_values(e, namespace) for e in group.expressions)

    _close(odes, _getter(cg, "odes", {}, 3))
    _close(radodes, _getter(cg, "radodes", {"order": 1}, 3))


@pytest.mark.parametrize(
    ("which", "kind", "options"),
    [
        ("tgas", "rates", {}),
        ("tgas", "odes", {}),
        ("rad", "radodes", {"order": 2}),
        ("tgas", "rhs", {"thermal": "dedt"}),
        ("rad", "jacobian", {"thermal": "dtdt", "radiation": True}),
    ],
)
def test_one_member_prints_as_its_getter(networks, which, kind, options):
    cg = Codegen(networks[which], lang="cxx")
    group = cg.cse_group([(kind, options)], cse_var="cse")
    printed = getattr(cg, GETTERS[kind])(**options, cse_var="cse")

    assert group.temporaries == printed["extras"]["cse"]
    assert group.expressions == (printed["expressions"],)


def test_rates_are_evaluated_once(networks):
    # Each exp(-5e4/tgas) rate is evaluated once for the rhs and the Jacobian.
    cg = Codegen(networks["tgas"], lang="python")
    group = cg.cse_group(
        [("rhs", {"thermal": "dtdt"}), ("jacobian", {"thermal": "dtdt"})]
    )
    code = [c for _, c in group.temporaries]
    code += [c for expressions in group.expressions for _, c in expressions]

    assert sum(c.count("math.exp") for c in code) == 1


def test_expressions_read_only_the_group_temporaries(networks):
    cg = Codegen(networks["tgas"], lang="python")
    group = cg.cse_group(
        [("rhs", {"thermal": "dtdt"}), ("jacobian", {"thermal": "dtdt"})]
    )
    names = {f"cse{idx}" for (idx,), _ in group.temporaries}
    read = {
        n
        for expressions in group.expressions
        for _, code in expressions
        for n in re.findall(r"\bcse\d+", code)
    }

    assert read and read <= names


def test_uses_cse_var_and_suffix(networks):
    cg = Codegen(networks["tgas"], lang="python")
    group = cg.cse_group(
        [("rhs", {"thermal": "dtdt"}), ("jacobian", {"thermal": "dtdt"})],
        cse_var="tmp",
        cse_suffix="_v",
    )
    namespace = _state(0)
    _define(group.temporaries, namespace, "tmp", "_v")

    assert all(re.fullmatch(r"\d+", str(idx[0])) for idx, _ in group.temporaries)
    assert _values(group.expressions[1], namespace)


def test_is_cached(networks):
    cg = Codegen(networks["tgas"], lang="python")
    members = [("rhs", {"thermal": "dtdt"}), ("jacobian", {"thermal": "dtdt"})]

    assert cg.cse_group(members) is cg.cse_group(list(members))


def test_unknown_kind_raises(networks):
    with pytest.raises(ValueError, match="kind"):
        Codegen(networks["tgas"], lang="python").cse_group([("fluxes", {})])


def test_empty_group_raises(networks):
    with pytest.raises(ValueError, match="member"):
        Codegen(networks["tgas"], lang="python").cse_group([])


def test_radodes_without_radiation_raise(networks):
    with pytest.raises(RuntimeError, match="radiation"):
        Codegen(networks["tgas"], lang="python").cse_group([("radodes", {})])
