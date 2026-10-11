# ABOUTME: Tests for .jaff network serialization and the SymPy JSON codec
# ABOUTME: Round-trip fidelity for reactions, rates, ODEs, and SymPy node types

import gzip
import json
from pathlib import Path

import pytest
import sympy
from sympy.tensor.indexed import Indexed
from sympy.matrices.expressions.matexpr import MatrixElement

from jaff import Network
from jaff.common._sympy_json import dumps, from_jsonable, loads, to_jsonable

FIXTURES = Path(__file__).parent / "fixtures"
REPO_ROOT = Path(__file__).resolve().parent.parent


# --------------------------------------------------------------------------- #
# Network .jaff serialization                                                  #
# --------------------------------------------------------------------------- #


def _assert_no_symbol_assumptions(node):
    """Rate expression nodes must not embed per-symbol assumptions."""
    if isinstance(node, list):
        if node and node[0] == "S" and len(node) > 2:
            raise AssertionError("Symbol node should not include assumptions")
        for item in node:
            _assert_no_symbol_assumptions(item)
    elif isinstance(node, dict):
        if node.get("type") == "Symbol" and "assumptions" in node:
            raise AssertionError("Symbol node should not include assumptions")
        for value in node.values():
            _assert_no_symbol_assumptions(value)


def test_network_json_roundtrip_sample_kida_valid(tmp_path):
    net = Network(str(FIXTURES / "sample_kida_valid.dat"))
    json_path = str(tmp_path / "net.jaff")
    net.to_jaff(json_path)

    # `.jaff` files are gzip-compressed by default.
    with open(json_path, "rb") as fb:
        assert fb.read(2) == b"\x1f\x8b"

    with gzip.open(json_path, "rt", encoding="utf-8") as f:
        payload = json.load(f)

    # rate_symbols carry the assumptions of every free symbol used in a rate.
    rate_symbols = payload.get("rate_symbols")
    assert isinstance(rate_symbols, list)
    rate_symbols_by_name = {
        item.get("name"): item.get("assumptions")
        for item in rate_symbols
        if isinstance(item, dict)
    }
    expected_symbols = {
        s
        for r in net.reactions
        if isinstance(r.rate, sympy.Basic)
        for s in r.rate.free_symbols
    }
    assert set(rate_symbols_by_name) == {s.name for s in expected_symbols}
    for sym in expected_symbols:
        expected = {
            k: v
            for k, v in (sym.assumptions0 or {}).items()
            if isinstance(k, str) and isinstance(v, bool)
        }
        assert rate_symbols_by_name.get(sym.name) == expected

    # Rate expression nodes themselves must stay assumption-free.
    for rj in payload.get("reactions") or []:
        rate_node = rj.get("rate")
        if isinstance(rate_node, dict) and rate_node.get("kind") == "string":
            continue
        if rate_node is not None:
            _assert_no_symbol_assumptions(rate_node)

    # Reload and compare reaction-by-reaction.
    net2 = Network(json_path)
    assert net2.label == net.label
    assert len(net2.species) == len(net.species)
    assert len(net2.reactions) == len(net.reactions)

    for r1, r2 in zip(net.reactions, net2.reactions):
        assert r2.get_verbatim() == r1.get_verbatim()
        assert (r2.tmin, r2.tmax) == (r1.tmin, r1.tmax)

        if isinstance(r1.rate, str):
            assert r2.rate == r1.rate
        else:
            assert isinstance(r2.rate, sympy.Basic)
            diff = sympy.simplify(r2.rate - r1.rate)
            free = sorted(diff.free_symbols, key=lambda s: s.name)
            if not free:
                r1_symbols = getattr(r1.rate, "free_symbols", set())
                if r1_symbols:
                    ref = abs(float(sympy.N(r1.rate.subs({s: 1.0 for s in r1_symbols}))))
                else:
                    ref = abs(float(sympy.N(r1.rate)))
                assert abs(float(diff.evalf())) <= 1e-15 * max(1.0, ref)
            else:
                for offset in (1.1, 10.1):
                    subs = {s: float(offset + i) for i, s in enumerate(free)}
                    val1 = float(sympy.N(r1.rate.subs(subs)))
                    val2 = float(sympy.N(r2.rate.subs(subs)))
                    assert abs(val2 - val1) <= 1e-12 * max(1.0, abs(val1))

        if isinstance(r1.dE, sympy.Basic) or isinstance(r2.dE, sympy.Basic):
            assert sympy.simplify(r2.dE - r1.dE) == 0
        else:
            assert r2.dE == r1.dE

    # Backward compatibility: a legacy *uncompressed* `.jaff` still loads.
    with gzip.open(json_path, "rt", encoding="utf-8") as f:
        raw = f.read()
    legacy_path = tmp_path / "legacy.jaff"
    legacy_path.write_text(raw)
    net3 = Network(str(legacy_path))
    assert len(net3.species) == len(net.species)
    assert len(net3.reactions) == len(net.reactions)


def test_network_json_roundtrip_preserves_nden_rates(tmp_path):
    """Rates standardized to ``nden[i, 0]`` MatrixElements survive a round-trip.

    Regression: the reload side once rebuilt the ``nden`` MatrixSymbol as a plain
    Symbol (breaking MatrixElement) and read ``dRad`` under the wrong JSON key.

    Uses a small self-contained network: the first rate references an electron
    density (``ne``) that standardizes to an ``nden[i, 0]`` MatrixElement, and a
    ``deltarad{i}`` auxiliary gives each reaction a non-trivial ``dRad``.  It
    deliberately avoids interpolation/undefined functions, which the ``.jaff``
    serializer does not support.
    """
    from sympy.matrices.expressions.matexpr import MatrixElement

    net_dir = tmp_path / "synnet"
    net_dir.mkdir()
    jet = net_dir / "syn.jet"
    jet.write_text(
        "H + H+ -> H+ + H               []         2.0e-9*ne\n"
        "H+ + E -> H                    []         2.63e-13*(Tgas/1e4)**(-0.7)*n_H\n"
    )
    (net_dir / "syn.jet.jfunc").write_text(
        "@function deltarad0()\n"
        "    @return 1.5e-11\n"
        "\n"
        "@function deltarad1()\n"
        "    @return 3.0e-11\n"
    )

    net = Network(str(jet))
    assert any("nden" in str(r.rate) for r in net.reactions)
    assert any(r.dRad != sympy.Float(0.0) for r in net.reactions)

    json_path = str(tmp_path / "gow.jaff")
    net.to_jaff(json_path)
    net2 = Network(json_path)  # must not raise

    assert [s.name for s in net2.species] == [s.name for s in net.species]
    assert len(net2.reactions) == len(net.reactions)

    # ODEs must match numerically; simplifying the symbolic difference is too
    # slow for large RHSs, so sample at a fixed point instead.
    odes1, odes2 = net.sodes(), net2.sodes()
    assert len(odes1) == len(odes2)

    names = set()
    for e in odes1 + odes2:
        names |= {
            str(t)
            for t in (e.atoms(sympy.Symbol) | e.atoms(MatrixElement) | e.atoms(Indexed))
        }
    sample = {n: 2.0 + i for i, n in enumerate(sorted(names))}

    def evaluate(expr):
        targets = (
            expr.atoms(sympy.Symbol) | expr.atoms(MatrixElement) | expr.atoms(Indexed)
        )
        return float(
            expr.xreplace({t: sympy.Float(sample[str(t)]) for t in targets}).evalf()
        )

    for e1, e2 in zip(odes1, odes2):
        v1, v2 = evaluate(e1), evaluate(e2)
        assert abs(v2 - v1) <= 1e-9 * max(1.0, abs(v1))


def test_network_json_roundtrip_preserves_dEdt_extra(tmp_path):
    """The ``heatingcoolingrate`` term (``thermodynamics.dEdt_extra``) survives a
    round-trip.

    Regression: ``to_jaff`` serialized only per-reaction energy terms, and the
    ``.jaff`` load path left the extra term at its ``Float(0.0)`` default, so the
    total thermal RHS silently lost the extra heating/cooling contribution.
    """
    net_dir = tmp_path / "heatnet"
    net_dir.mkdir()
    jet = net_dir / "heat.jet"
    jet.write_text("H + H -> H2                    []         1.0e-17\n")
    (net_dir / "heat.jet.jfunc").write_text(
        "@function heatingCoolingRate(tgas)\n    @return 2*tgas\n"
    )

    net = Network(str(jet))
    extra = net.thermodynamics.dEdt_extra.volumetric
    assert extra != sympy.Float(0.0)

    json_path = str(tmp_path / "heat.jaff")
    net.to_jaff(json_path)
    net2 = Network(json_path)

    # Total extra thermal RHS must match, not just per-reaction dE terms.
    extra2 = net2.thermodynamics.dEdt_extra.volumetric
    diff = sympy.simplify(extra2 - extra)
    assert diff == 0, f"dEdt_extra lost on round-trip: {extra} -> {extra2}"


# --------------------------------------------------------------------------- #
# Corrupted / malformed .jaff payloads must be rejected, not silently aliased  #
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


def test_jaff_rejects_missing_real_participant(tmp_path):
    """A reaction referencing a real species absent from the table is rejected.

    Regression: the missing name was silently turned into an out-of-catalogue
    pseudo-species at index -1, aliasing its stoichiometry onto the last real
    species via Python negative indexing and corrupting the generated ODE.
    """
    path = _make_valid_jaff(tmp_path)
    payload = _read_payload(path)
    payload["species"] = [s for s in payload["species"] if s["name"] != "H2"]
    _write_payload(path, payload)

    with pytest.raises(ValueError, match="H2"):
        Network(str(path), funcfile=False)


def test_jaff_retains_valid_pseudo_species(tmp_path):
    """An underscore-prefixed pseudo-species absent from the table still loads."""
    path = _make_valid_jaff(tmp_path)
    payload = _read_payload(path)
    # Add a real _PHOTON-style pseudo participant to the reaction's reactants.
    payload["reactions"][0]["reactants"].append("_PHOTON")
    _write_payload(path, payload)

    net = Network(str(path), funcfile=False)
    # Pseudo-species stay out of the integrated catalogue.
    assert "_PHOTON" not in {s.name for s in net.species}


def test_jaff_rejects_duplicate_species_names(tmp_path):
    path = _make_valid_jaff(tmp_path)
    payload = _read_payload(path)
    payload["species"].append({"name": "H", "index": 2, "mass": 1.0, "charge": 0})
    _write_payload(path, payload)

    with pytest.raises(ValueError, match="[Dd]uplicate"):
        Network(str(path), funcfile=False)


def test_jaff_rejects_negative_species_index(tmp_path):
    path = _make_valid_jaff(tmp_path)
    payload = _read_payload(path)
    payload["species"][0]["index"] = -1
    _write_payload(path, payload)

    with pytest.raises(ValueError):
        Network(str(path), funcfile=False)


def test_jaff_rejects_species_index_gap(tmp_path):
    path = _make_valid_jaff(tmp_path)
    payload = _read_payload(path)
    # Two species indexed 0 and 2 (gap at 1) — not contiguous 0..N-1.
    payload["species"][-1]["index"] = payload["species"][-1]["index"] + 1
    _write_payload(path, payload)

    with pytest.raises(ValueError):
        Network(str(path), funcfile=False)


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
