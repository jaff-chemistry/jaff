# ABOUTME: CSE temporaries in REPEAT blocks: the $cse$ line's identifier pattern
# ABOUTME: (prefix$idx$suffix) must name both declarations and references, or be rejected

import math
import re
from pathlib import Path
from typing import Any, Dict, List, Set

import pytest

from jaff import Network
from jaff.codegen import Codegen
from jaff.codegen import TemplateParser
from jaff.errors import ParserError
from tests._paths import FIXTURES_DIR

FIXTURES = FIXTURES_DIR
TGAS = 321.0

# prop -> (REPEAT variable, main line using it)
PROPS = {
    "rates": ("rate", "f[$idx$] = $rate$"),
    "odes": ("ode", "f[$idx$] = $ode$"),
    "rhses": ("rhs", "f[$idx$] = $rhs$"),
    "jacobian": ("expr", "f[$idx$, $idx$] = $expr$"),
}


@pytest.fixture(scope="module")
def cse_net() -> Network:
    """Rates share tgas**0.5, exp(-100/tgas), ... across reactions."""
    return Network(str(FIXTURES / "test_cse.dat"))


def _render(net: Network, tmp_path: Path, header: str, body: str) -> str:
    template = tmp_path / "t.py"
    template.write_text(f"$$JAFF REPEAT {header}\n{body}\n$$JAFF END\n")
    return TemplateParser(net, template).parse_file()


def _exec_rates(source: str) -> Dict[int, float]:
    scope: Dict[str, Any] = {"math": math, "tgas": TGAS, "k": {}}
    exec(source, scope)
    return scope["k"]


def _assert_names_consistent(code: str, pattern: str) -> None:
    """Every temporary matching *pattern* that is referenced is also declared."""
    ident = re.compile(pattern)
    declared: Set[str] = set()
    referenced: Set[str] = set()
    for line in code.splitlines():
        lhs, sep, rhs = line.partition("=")
        if not sep:
            continue
        m = ident.fullmatch(lhs.strip())
        if m:
            declared.add(m.group(0))
        referenced |= set(ident.findall(rhs))
    assert declared, f"no temporaries matching {pattern!r} declared:\n{code}"
    assert referenced, f"no temporaries matching {pattern!r} referenced:\n{code}"
    assert referenced <= declared, f"undeclared: {referenced - declared}\n{code}"


@pytest.mark.parametrize(
    "cse_line, rate_line, pattern",
    [
        ("cse$idx$ = $cse$", "k[$idx$] = $rate$", r"\bcse\d+\b"),
        ("tmp$idx$_value = $cse$", "k[$idx$] = $rate$", r"\btmp\d+_value\b"),
        ("cse$idx$=$cse$", "k[$idx$]=$rate$", r"\bcse\d+\b"),
        ("tmp$idx$_value=$cse$", "k[$idx$]=$rate$", r"\btmp\d+_value\b"),
        ("t_$idx$_v2 = $cse$", "k[$idx$] = $rate$", r"\bt_\d+_v2\b"),
        ("s.x$idx$ = $cse$", "k[$idx$] = $rate$", r"\bx\d+\b"),
    ],
)
def test_cse_names_match_template_pattern(
    cse_net: Network, tmp_path: Path, cse_line: str, rate_line: str, pattern: str
) -> None:
    out = _render(
        cse_net, tmp_path, "idx, rate, cse IN rates", f"{cse_line}\n{rate_line}"
    )
    _assert_names_consistent(out.replace("s.x", "x"), pattern)


@pytest.mark.parametrize(
    "cse_line",
    ["cse$idx$ = $cse$", "tmp$idx$_value = $cse$", "cse$idx$=$cse$", "a$idx$b=$cse$"],
)
def test_cse_rates_evaluate_like_plain_rates(
    cse_net: Network, tmp_path: Path, cse_line: str
) -> None:
    plain = _exec_rates(
        _render(cse_net, tmp_path, "idx, rate IN rates", "k[$idx$] = $rate$")
    )
    rate_line = "k[$idx$]=$rate$" if "=$" in cse_line else "k[$idx$] = $rate$"
    source = _render(
        cse_net, tmp_path, "idx, rate, cse IN rates", f"{cse_line}\n{rate_line}"
    )
    with_cse = _exec_rates(source)
    assert with_cse.keys() == plain.keys()
    for i, value in plain.items():
        assert with_cse[i] == pytest.approx(value, rel=1e-12)


@pytest.mark.parametrize("prop", list(PROPS))
def test_suffix_pattern_consistent_for_every_prop(
    cse_net: Network, tmp_path: Path, prop: str
) -> None:
    var, main_line = PROPS[prop]
    body = f"tmp$idx$_value = $cse$\n{main_line}"
    out = _render(cse_net, tmp_path, f"idx, {var}, cse IN {prop}", body)
    _assert_names_consistent(out, r"\btmp\d+_value\b")


@pytest.mark.parametrize(
    "body, match",
    [
        # identifier would start with a digit: 0x, 1x, ...
        ("$idx$x = $cse$\nk[$idx$] = $rate$", "identifier"),
        # array-style temporaries are not identifiers
        ("cse[$idx$] = $cse$\nk[$idx$] = $rate$", "identifier"),
        # arithmetic on the index would declare cse1 for temporary cse0
        ("cse$idx+1$ = $cse$\nk[$idx$] = $rate$", "arithmetic"),
        # temporaries must be emitted (and named) before they are used
        ("k[$idx$] = $rate$\ncse$idx$ = $cse$", r"\$cse\$"),
    ],
)
def test_unsupported_cse_patterns_rejected(
    cse_net: Network, tmp_path: Path, body: str, match: str
) -> None:
    with pytest.raises(ParserError, match=match):
        _render(cse_net, tmp_path, "idx, rate, cse IN rates", body)


def test_horizontal_rates_without_cse_do_not_need_idx(
    cse_net: Network, tmp_path: Path
) -> None:
    out = _render(cse_net, tmp_path, "rate IN rates", "k = [$rate$]")
    line = next(line for line in out.splitlines() if line.startswith("k = "))
    assert line.count("exp(") == cse_net.reactions.count


GETTERS: List[str] = [
    "get_indexed_rates",
    "get_indexed_odes",
    "get_indexed_rhs",
    "get_indexed_jacobian",
]


@pytest.mark.parametrize("getter", GETTERS)
def test_codegen_cse_suffix_names_temporaries(cse_net: Network, getter: str) -> None:
    out = getattr(Codegen(cse_net, lang="python"), getter)(
        use_cse=True, cse_var="t", cse_suffix="_v"
    )
    declared = {f"t{idx[0]}_v" for idx, _ in out["extras"]["cse"]}
    code = " ".join(str(e) for _, e in out["extras"]["cse"])
    code += " " + " ".join(str(e) for _, e in out["expressions"])
    referenced = set(re.findall(r"\bt\d+_v\b", code))
    assert referenced and referenced <= declared
    assert not re.search(r"\bt\d+\b", code)
