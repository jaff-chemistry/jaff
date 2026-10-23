# ABOUTME: Tests for REPEAT's vertical and horizontal line expansion functions
# ABOUTME: Plain lists and IndexedLists in; expanded line contents out

import pytest

from jaff.codegen.templating.commands.repeat.expand import (
    as_flattened,
    expand_horizontal,
    expand_vertical,
    find_idx,
)
from jaff.codegen.templating.errors import TemplateError
from jaff.types import IndexedList


def vertical(items, line, token="$x$", indent=""):
    lines = expand_vertical(as_flattened(items), line, token, 0)
    return "".join(f"{indent}{text}\n" for text in lines)


def test_vertical_writes_one_line_per_item_with_indent():
    assert vertical(["a", "b"], "y[$idx$] = $x$", indent="  ") == (
        "  y[0] = a\n  y[1] = b\n"
    )


def test_vertical_applies_offsets():
    assert vertical(["a", "b"], "y[$idx+1$] = $x$") == "y[1] = a\ny[2] = b\n"
    assert vertical(["a"], "z[$idx-1$] = $x$") == "z[-1] = a\n"


def test_vertical_upper_idx_takes_first_index_and_offset():
    assert vertical(["a"], "$IDX$ $idx+2$ $x$") == "2 2 a\n"


def test_vertical_matrix_needs_one_idx_per_dimension():
    assert vertical([[1, 2], [3]], "m[$idx$][$idx$] = $x$") == (
        "m[0][0] = 1\nm[0][1] = 2\nm[1][0] = 3\n"
    )


def test_vertical_idx_count_mismatch_is_an_error_at_the_line():
    with pytest.raises(TemplateError, match="expected 2 \\$idx\\$ placeholders") as e:
        expand_vertical(as_flattened([[1]]), "m[$idx$] = $x$", "$x$", 3)
    assert e.value.body_index == 3


def test_vertical_accepts_indexed_lists():
    items = IndexedList(["p", "q"])
    assert vertical(items, "$idx$:$x$") == "0:p\n1:q\n"


@pytest.mark.parametrize(
    "line, expected",
    [
        ('n = ["$x$", ]', 'n = ["a", "b"]'),
        ("n = ('$x$', )", "n = ('a', 'b')"),
        ("n = {$x$; }", "n = {a; b}"),
        ("n = [$x$]", "n = [a, b]"),
        ("n = (/ $x$, /)", "n = (/a, b/)"),
    ],
)
def test_horizontal_takes_brackets_quotes_and_separator_from_the_line(line, expected):
    assert expand_horizontal(["a", "b"], "x", line) == expected


def test_horizontal_nests_sub_lists():
    assert expand_horizontal([[1, 2], [3]], "x", "m = [$x$, ]") == "m = [[1, 2], [3]]"


def test_horizontal_without_a_bracketed_placeholder_gives_none():
    assert expand_horizontal(["a"], "x", "n = $x$") is None


def test_find_idx_gives_spans_and_arithmetic():
    assert find_idx("a$idx$b$idx+2$$idx - 1$") == [
        ((1, 6), None, 0),
        ((7, 14), "+", 2),
        ((14, 23), "-", 1),
    ]


@pytest.mark.parametrize(
    "line, expected",
    [("$idx*2$", "4"), ("$idx / 2$", "1"), ("$idx + 1$", "3"), ("$IDX$ $idx*3$", "6 6")],
)
def test_idx_takes_the_same_arithmetic_as_other_placeholders(line, expected):
    assert expand_vertical([([2], "a")], line, "$x$", 0) == [expected]


def test_empty_items_expand_to_nothing_or_an_empty_list():
    assert vertical([], "y[$idx$] = $x$") == ""
    assert expand_horizontal([], "x", "a = [$x$, ]") == "a = []"


def test_upper_idx_without_idx_tokens_is_the_first_index():
    items = [([3], "a")]
    assert expand_vertical(items, "$IDX$ $x$", "$x$", 0) == ["3 a"]


def test_line_without_idx_skips_the_dimension_check():
    assert vertical([[1, 2]], "v = $x$") == "v = 1\nv = 2\n"


def test_multi_digit_offsets():
    assert vertical(["a"], "y[$idx+10$]$x$") == "y[10]a\n"


def test_horizontal_accepts_flattened_indexed_lists():
    flat = as_flattened([[1, 2], [3]])
    assert expand_horizontal(flat, "x", "m = ['$x$', ]") == "m = [['1', '2'], ['3']]"
