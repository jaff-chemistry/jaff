# ABOUTME: Tests for option resolution: USE layering, group checks and REPLACE handling
# ABOUTME: Pure inputs built with the header parser; no network or codegen is involved

from pathlib import Path

import pytest

from jaff.codegen.templating.errors import TemplateError
from jaff.codegen.templating.resolve import (
    apply_replacements,
    check_groups,
    compile_replacements,
    layer,
)
from jaff.codegen.templating.syntax import parse_header

FILE = Path("template.cpp")


def options(text: str, line: int = 1):
    """Options of ``$$JAFF SUB x { <text> }`` starting at template *line*."""
    return parse_header(f"$$JAFF SUB x {{{text}}}", line, FILE).options


def test_no_block_gives_empty_options():
    assert layer(None, {}).groups == {}


def test_use_sets_apply_in_order_then_local():
    defines = {
        "a": options(' REPLACE { "1" -> "a" } SORT { ENABLED -> True } '),
        "b": options(' REPLACE { "2" -> "b" } SORT { ENABLED -> False } '),
    }
    result = layer(options(' USE a USE b REPLACE { "3" -> "c" } '), defines)
    assert result.replace == [("1", "a"), ("2", "b"), ("3", "c")]
    assert result.groups == {"SORT": {"ENABLED": "False"}}
    assert result.uses == []


def test_local_keys_win_and_other_keys_merge():
    defines = {"m": options(" THERMAL { MODE -> dedt  ENERGY -> specific } ")}
    result = layer(options(" USE m THERMAL { MODE -> dtdt } "), defines)
    assert result.groups == {"THERMAL": {"MODE": "dtdt", "ENERGY": "specific"}}


def test_undefined_use_reports_the_use_line():
    with pytest.raises(TemplateError, match="undefined option set m") as err:
        layer(options("\n USE m\n", line=4), {})
    assert err.value.nline == 5


def test_groups_from_use_report_the_use_line():
    defines = {"m": options(" RADIATION { ENABLED -> True } ", line=1)}
    result = layer(options("\n\n USE m\n", line=10), defines)
    assert result.group_lines == {"RADIATION": 12}


def test_check_groups_accepts_known_groups_and_keys():
    check_groups("REPEAT", options(" SORT { ENABLED -> True } "), {"SORT": {"ENABLED"}})


def test_check_groups_rejects_unknown_group_with_its_line():
    with pytest.raises(TemplateError, match="SUB takes no option group SORT") as err:
        check_groups("SUB", options("\n SORT { ENABLED -> True }\n", line=7), {})
    assert err.value.nline == 8


def test_check_groups_rejects_unknown_key():
    with pytest.raises(TemplateError, match="SORT has no setting ON"):
        check_groups("REPEAT", options(" SORT { ON -> True } "), {"SORT": {"ENABLED"}})


def test_replacements_apply_in_order():
    compiled = compile_replacements(options(' REPLACE { "a" -> "b"  "b" -> "c" } '))
    assert apply_replacements("a", compiled) == "c"


def test_replacements_keep_regex_groups():
    compiled = compile_replacements(options(r' REPLACE { "n\[(\d+)\]" -> "y[\1]" } '))
    assert apply_replacements("n[3] + n[4]", compiled) == "y[3] + y[4]"


def test_invalid_pattern_reports_its_line():
    with pytest.raises(TemplateError, match="Invalid REPLACE pattern") as err:
        compile_replacements(options('\n REPLACE {\n "(" -> "x" } ', line=3))
    assert err.value.nline == 5


def test_replacement_with_missing_group_reports_its_line():
    with pytest.raises(TemplateError, match="Invalid REPLACE replacement") as err:
        compile_replacements(options(r' REPLACE { "tgas" -> "T\1" } ', line=2))
    assert err.value.nline == 2


def test_check_groups_rejects_group_from_use_at_the_use_line():
    defines = {"m": options(" RADIATION { ENABLED -> True } ", line=1)}
    result = layer(options("\n\n USE m\n", line=10), defines)

    with pytest.raises(TemplateError, match="takes no option group RADIATION") as err:
        check_groups("SUB", result, {})

    assert err.value.nline == 12


def test_replace_pairs_from_use_report_the_use_line():
    defines = {"m": options(' REPLACE { "a" -> "b"  "c" -> "d" } ', line=1)}
    result = layer(options('\n USE m\n REPLACE { "e" -> "f" }\n', line=20), defines)
    assert result.replace_lines == [21, 21, 22]
