# ABOUTME: REPEAT option groups end to end: booleans are case-insensitive True/False,
# ABOUTME: ints stay ints, strings stay strings, and anything else is a ParserError

from pathlib import Path
from typing import List

import pytest

from jaff import Network
from jaff.codegen import TemplateParser
from jaff.errors import ParserError
from tests._paths import FIXTURES_DIR

FIXTURES = FIXTURES_DIR
TRUE_SPELLINGS = ["True", "TRUE", "true"]
FALSE_SPELLINGS = ["False", "FALSE", "false"]

# Boolean option "GROUP KEY" -> REPEAT line and body exercising it.  All default to False.
BOOL_CASES = {
    "RADIATION ENABLED": ("idx, rhs IN rhses", "f[$idx$] = $rhs$"),
}


@pytest.fixture(scope="module")
def dedt_net() -> Network:
    """Energy-carrying network with no radiation configured."""
    return Network(str(FIXTURES / "test_jac_dedt.dat"))


@pytest.fixture(scope="module")
def charged_net() -> Network:
    return Network(str(FIXTURES / "react_cie_hepp.jet"))


def _option(setting: str, value: str) -> str:
    """Option group text for ``"GROUP KEY"`` set to *value*."""
    group, key = setting.split()
    return f"{group} {{ {key} -> {value} }}"


def _render(net: Network, tmp_path: Path, repeat: str, body: str, opts: str) -> List[str]:
    block = f" {{ {opts} }}" if opts else ""
    template = tmp_path / "t.py"
    template.write_text(f"$$JAFF REPEAT {repeat}{block}\n{body}\n$$JAFF END\n")
    lines = TemplateParser(net, template).parse_file().splitlines()
    return [line for line in lines if line.startswith("f")]


@pytest.mark.parametrize("spelling", FALSE_SPELLINGS)
@pytest.mark.parametrize("modifier", list(BOOL_CASES))
def test_false_spellings_match_default(
    dedt_net: Network, tmp_path: Path, modifier: str, spelling: str
) -> None:
    repeat, body = BOOL_CASES[modifier]
    default = _render(dedt_net, tmp_path, repeat, body, "")
    opts = _option(modifier, spelling)
    assert _render(dedt_net, tmp_path, repeat, body, opts) == default


@pytest.mark.parametrize("spelling", TRUE_SPELLINGS)
def test_radiation_true_still_requires_radiation(
    dedt_net: Network, tmp_path: Path, spelling: str
) -> None:
    repeat, body = BOOL_CASES["RADIATION ENABLED"]
    with pytest.raises(ParserError, match="No radiation bands") as err:
        _render(dedt_net, tmp_path, repeat, body, _option("RADIATION ENABLED", spelling))
    assert isinstance(err.value.__cause__, RuntimeError)


@pytest.mark.parametrize("value", ["1", "0", "yes", "None", "'True'"])
@pytest.mark.parametrize("modifier", list(BOOL_CASES))
def test_non_boolean_values_rejected(
    dedt_net: Network, tmp_path: Path, modifier: str, value: str
) -> None:
    repeat, body = BOOL_CASES[modifier]
    with pytest.raises(ParserError, match=f"{modifier} expects True or False"):
        _render(dedt_net, tmp_path, repeat, body, _option(modifier, value))


def test_dedt_type_selects_energy_form(dedt_net: Network, tmp_path: Path) -> None:
    repeat, body = "idx, rhs IN rhses", "f[$idx$] = $rhs$"
    default = _render(dedt_net, tmp_path, repeat, body, "")
    specific = _render(dedt_net, tmp_path, repeat, body, "THERMAL { ENERGY -> specific }")
    per_particle = _render(
        dedt_net, tmp_path, repeat, body, "THERMAL { ENERGY -> per_particle }"
    )
    volumetric = "THERMAL { ENERGY -> volumetric }"
    assert _render(dedt_net, tmp_path, repeat, body, volumetric) == default
    assert specific != default
    assert per_particle != specific


def test_unknown_dedt_type_rejected(dedt_net: Network, tmp_path: Path) -> None:
    repeat, body = "idx, rhs IN rhses", "f[$idx$] = $rhs$"
    with pytest.raises(ParserError, match="THERMAL ENERGY expects one of .*'bogus'"):
        _render(dedt_net, tmp_path, repeat, body, "THERMAL { ENERGY -> bogus }")


@pytest.mark.parametrize("pos, expected", [("p", "idx_hep"), ("1", "idx_he1")])
def test_pos_neg_stay_strings(
    charged_net: Network, tmp_path: Path, pos: str, expected: str
) -> None:
    lines = _render(
        charged_net,
        tmp_path,
        "idx, specie_with_normalized_sign IN species_with_normalized_sign",
        "f_idx_$specie_with_normalized_sign$ = $idx$",
        f"SIGN {{ POS -> {pos}  NEG -> m }}",
    )
    assert any(line.startswith(f"f_{expected} =") for line in lines), lines


JAC = ("idx, expr IN jacobian", "f[$idx$, $idx$] = $expr$")
RHS = ("idx, rhs IN rhses", "f[$idx$] = $rhs$")


@pytest.mark.parametrize("repeat, body", [JAC, RHS])
def test_thermal_modes_change_output(
    dedt_net: Network, tmp_path: Path, repeat: str, body: str
) -> None:
    out = {
        mode: _render(dedt_net, tmp_path, repeat, body, f"THERMAL {{ MODE -> {mode} }}")
        for mode in ("none", "dedt", "dtdt")
    }
    assert out["none"] != out["dedt"] != out["dtdt"]


def test_thermal_rejects_unknown_mode(dedt_net: Network, tmp_path: Path) -> None:
    with pytest.raises(ParserError, match="THERMAL"):
        _render(dedt_net, tmp_path, *JAC, "THERMAL { MODE -> bogus }")


def test_use_dedt_is_gone(dedt_net: Network, tmp_path: Path) -> None:
    # Unknown option groups are template errors (bug 7: the old engine raised KeyError)
    with pytest.raises(ParserError, match="takes no option group USE_DEDT"):
        _render(dedt_net, tmp_path, *JAC, "USE_DEDT { ENABLED -> True }")


def test_dtdt_token_renders(dedt_net: Network, tmp_path: Path) -> None:
    template = tmp_path / "t.py"
    template.write_text("$$JAFF SUB dtdt\nx = $dtdt$\n$$JAFF END\n")
    out = TemplateParser(dedt_net, template).parse_file()
    assert "tgas" in out and "$dtdt$" not in out
