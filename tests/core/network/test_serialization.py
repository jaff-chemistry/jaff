# ABOUTME: Tests for .jaff network serialization and the SymPy JSON codec
# ABOUTME: Writing .jaff payloads (loading is unsupported) and SymPy codec round trips

import gzip
import json

import pytest
import sympy
from sympy.matrices.expressions.matexpr import MatrixElement

from jaff import Network
from jaff.common._sympy_json import dumps, from_jsonable, loads, to_jsonable
from tests._paths import FIXTURES_DIR, REPO_ROOT

FIXTURES = FIXTURES_DIR
REPO_ROOT = REPO_ROOT


# --------------------------------------------------------------------------- #
# Network .jaff serialization                                                  #
# --------------------------------------------------------------------------- #


def _make_valid_jaff(tmp_path, name="corrupt.jaff"):
    """Serialize the minimal valid ``H + H -> H2`` network to a ``.jaff`` file."""
    dat = tmp_path / "src.dat"
    dat.write_text("@format:idx,R,R,P,rate\n1,H,H,H2,1\n")
    net = Network(str(dat), funcfile=False)
    path = tmp_path / name
    net.to_jaff(str(path))
    return path


def _read_payload(path):
    with gzip.open(path, "rt", encoding="utf-8") as f:
        return json.load(f)


def _write_payload(path, payload):
    with gzip.open(path, "wt", encoding="utf-8") as f:
        json.dump(payload, f)


def test_to_jaff_writes_gzip_json(make_network, tmp_path):
    net = make_network("H + H -> H2 [] 1e-17\n", name="r.jet")
    net.to_jaff(tmp_path / "r.jaff")

    with open(tmp_path / "r.jaff", "rb") as fb:
        assert fb.read(2) == b"\x1f\x8b"

    payload = _read_payload(tmp_path / "r.jaff")

    assert payload["format"] == "jaff.network_json"
    assert payload["schema_version"] == 1
    assert [r["format"] for r in payload["reactions"]] == ["prizmo"]


def test_loading_a_jaff_file_is_not_supported(make_network, tmp_path):
    from jaff.errors import ParserError

    make_network("H + H -> H2 [] 1e-17\n", name="r.jet").to_jaff(tmp_path / "r.jaff")

    with pytest.raises(ParserError, match="loading .jaff files is not supported"):
        Network(str(tmp_path / "r.jaff"))


def test_jaff_does_not_store_drad(make_network, tmp_path):
    """Reactions carry no radiation source term, so ``.jaff`` stores none."""
    net = make_network("H + H -> H2 [] 1e-17\n", name="r.jet")
    net.to_jaff(tmp_path / "r.jaff")

    payload = _read_payload(tmp_path / "r.jaff")

    assert all("dRad" not in r for r in payload["reactions"])


# --------------------------------------------------------------------------- #
# SymPy JSON codec                                                             #
# --------------------------------------------------------------------------- #


def _rt(expr: sympy.Basic) -> sympy.Basic:
    return loads(dumps(expr))


def test_roundtrip_basic_arithmetic_and_numbers():
    x = sympy.Symbol("x", real=True)
    y = sympy.Symbol("y")
    expr = sympy.Add(
        sympy.Integer(2),
        sympy.Rational(1, 3),
        sympy.Float("1.0e-10"),
        x * y,
        evaluate=False,
    )
    expr2 = _rt(expr)
    assert abs(float(sympy.simplify(expr2 - expr).evalf())) < 1e-15
    x2 = next(s for s in expr2.free_symbols if s.name == "x")
    assert x2.assumptions0.get("real") is True


def test_roundtrip_commutative_order_is_deterministic():
    x, y = sympy.Symbol("x"), sympy.Symbol("y")
    expr1 = sympy.Add(x, y, sympy.Integer(1), evaluate=False)
    expr2 = sympy.Add(y, sympy.Integer(1), x, evaluate=False)
    j1, j2 = dumps(expr1), dumps(expr2)
    assert j1 != j2
    assert loads(j1) == expr1
    assert loads(j2) == expr2


def test_roundtrip_piecewise_strict_lessthan():
    x = sympy.Symbol("x")
    pw = sympy.Piecewise(
        (x, sympy.StrictLessThan(x, sympy.Integer(0))),
        (sympy.Integer(0), sympy.true),
        evaluate=False,
    )
    assert _rt(pw) == pw


def test_roundtrip_min_max_log_exp_pow_mul():
    x = sympy.Symbol("x")
    expr = sympy.Mul(
        sympy.Max(x, sympy.Integer(1), evaluate=False),
        sympy.Min(x, sympy.Integer(2), evaluate=False),
        sympy.exp(x),
        sympy.log(sympy.Integer(10)),
        sympy.Pow(x, sympy.Integer(2), evaluate=False),
        evaluate=False,
    )
    assert _rt(expr) == expr


def test_roundtrip_matrix_symbol_and_element():
    nden = sympy.MatrixSymbol("nden", sympy.Integer(3), sympy.Integer(1))
    expr = nden[0, 0] + nden[2, 0]
    assert _rt(expr) == expr


def test_direct_jsonable_api_roundtrip():
    x = sympy.Symbol("x", positive=True)
    expr = sympy.Add(sympy.Integer(1), sympy.exp(x), evaluate=False)
    assert from_jsonable(to_jsonable(expr)) == expr


def test_compact_float_is_number():
    expr = sympy.Add(sympy.Float("2.0e-12"), sympy.Float("1.5"), evaluate=False)
    node = to_jsonable(expr)
    assert isinstance(node, list) and node[0] == "Add"
    assert all(isinstance(arg, (int, float)) for arg in node[1])
    assert abs(float(sympy.simplify(from_jsonable(node) - expr).evalf())) < 1e-15
