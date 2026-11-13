# ABOUTME: Tests for Codegen.get_indexed_rhs_jacobian: the rhs and Jacobian sharing one
# ABOUTME: set of temporaries, run as Python and compared with the separate getters

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

CASES = [
    pytest.param("tgas", {}, id="none"),
    pytest.param("tgas", {"thermal": "dtdt"}, id="dtdt"),
    *(
        pytest.param("tgas", {"thermal": "dedt", "energy": e}, id=f"dedt-{e}")
        for e in ("volumetric", "specific", "per_particle", "molar")
    ),
    pytest.param("photo", {"thermal": "dtdt"}, id="photorates"),
    pytest.param("rad", {"thermal": "dtdt", "radiation": True, "rad_order": 2}, id="rad"),
]


@pytest.fixture(scope="module")
def networks(tmp_path_factory):
    path = tmp_path_factory.mktemp("rhs_jac") / "tgas.dat"
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


def _run(printed: dict, namespace: dict, cse_var: str = "cse") -> dict:
    """Run *printed*'s temporaries into *namespace*; its expressions by index."""
    for (idx,), code in printed["extras"]["cse"]:
        namespace[f"{cse_var}{idx}"] = eval(code, namespace)

    return {tuple(idx): eval(code, namespace) for idx, code in printed["expressions"]}


def _close(got: dict, want: dict) -> None:
    assert set(got) == set(want)
    scale = max(abs(v) for v in want.values())
    for key, value in want.items():
        assert got[key] == pytest.approx(value, rel=1e-9, abs=1e-12 * scale), key


@pytest.mark.parametrize(("which", "options"), CASES)
@pytest.mark.parametrize("seed", [0, 1])
def test_matches_the_separate_rhs_and_jacobian(networks, which, options, seed) -> None:
    cg = Codegen(networks[which], lang="python")
    namespace = _state(seed)
    rhs = _run(cg.get_indexed_rhs_jacobian("rhs", **options), namespace)
    jac = _run(cg.get_indexed_rhs_jacobian("jacobian", **options), namespace)

    _close(rhs, _run(cg.get_indexed_rhs(**options), _state(seed)))
    _close(jac, _run(cg.get_indexed_jacobian(**options), _state(seed)))


def test_jacobian_part_reads_the_rhs_temporaries(networks) -> None:
    cg = Codegen(networks["tgas"], lang="python")
    rhs = cg.get_indexed_rhs_jacobian("rhs", thermal="dtdt")
    jac = cg.get_indexed_rhs_jacobian("jacobian", thermal="dtdt")
    names = {f"cse{idx}" for (idx,), _ in rhs["extras"]["cse"]}
    read = {n for _, code in jac["expressions"] for n in re.findall(r"cse\d+", code)}

    assert rhs["extras"]["cse"]
    assert not jac["extras"]["cse"]
    assert read and read <= names


def test_rates_are_not_written_out_twice(networks) -> None:
    # Each exp(-5e4/tgas) rate is evaluated once for the rhs and the Jacobian.
    cg = Codegen(networks["tgas"], lang="python")
    printed = cg.get_indexed_rhs_jacobian("rhs", thermal="dtdt")
    code = [c for _, c in printed["extras"]["cse"]]
    code += [
        c
        for _, c in cg.get_indexed_rhs_jacobian("jacobian", thermal="dtdt")["expressions"]
    ]

    assert sum(c.count("math.exp") for c in code) == 1


def test_uses_cse_var_and_suffix(networks) -> None:
    cg = Codegen(networks["tgas"], lang="python")
    printed = cg.get_indexed_rhs_jacobian(
        "jacobian", thermal="dtdt", cse_var="tmp", cse_suffix="_v"
    )

    assert all(
        name.startswith("tmp") and name.endswith("_v")
        for _, code in printed["expressions"]
        for name in re.findall(r"\btmp\w*", code)
    )


def test_is_cached(networks) -> None:
    cg = Codegen(networks["tgas"], lang="python")

    assert cg.get_indexed_rhs_jacobian("rhs", thermal="dtdt") is (
        cg.get_indexed_rhs_jacobian("rhs", thermal="dtdt")
    )


def test_unknown_part_raises(networks) -> None:
    with pytest.raises(ValueError, match="part"):
        Codegen(networks["tgas"], lang="python").get_indexed_rhs_jacobian("both")


def test_needs_cse(networks) -> None:
    with pytest.raises(ValueError, match="temporaries"):
        Codegen(networks["tgas"], lang="python").get_indexed_rhs_jacobian(
            "rhs", use_cse=False
        )
