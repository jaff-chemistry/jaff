# ABOUTME: Tests for $$JAFF CALL: text from a @jaffgen_function in
# ABOUTME: jaffgen_template_functions.py, indented, REPLACEd, loaded lazily once

import textwrap

import pytest

from jaff import Network
from jaff.codegen import TemplateFunctions, jaffgen_function
from jaff.codegen.templating import TemplateParser
from jaff.errors import ParserError
from tests._paths import FIXTURES_DIR

FUNCTIONS = "jaffgen_template_functions.py"


@pytest.fixture(scope="module")
def net() -> Network:
    return Network(str(FIXTURES_DIR / "react_cie_hepp.jet"))


def functions(tmp_path, source: str) -> TemplateFunctions:
    """A loader for a functions file holding *source* (dedented)."""
    path = tmp_path / FUNCTIONS
    path.write_text(textwrap.dedent(source))

    return TemplateFunctions(path)


def render(net, tmp_path, text, funcs, name="t.cpp", echo=False) -> str:
    template = tmp_path / name
    template.write_text(text)

    return TemplateParser(
        net, template, echo_directives=echo, functions=funcs
    ).parse_file()


LINES = """
    from jaff.codegen import jaffgen_function

    @jaffgen_function
    def lines(net):
        return f"a = {len(net.species)};\\nb = 2;"
"""


# --------------------------------------------------------------------------- #
# Rendering                                                                    #
# --------------------------------------------------------------------------- #


def test_call_inserts_the_returned_text(net, tmp_path):
    out = render(net, tmp_path, "x\n$$JAFF CALL lines\ny\n", functions(tmp_path, LINES))

    assert out == "x\na = 5;\nb = 2;\ny\n"


def test_returned_lines_get_the_call_line_indentation(net, tmp_path):
    funcs = functions(
        tmp_path,
        """
        from jaff.codegen import jaffgen_function

        @jaffgen_function
        def block(net):
            return "if (x) {\\n    y();\\n\\n}\\n"
        """,
    )
    out = render(net, tmp_path, "\t  $$JAFF CALL block\n", funcs)

    assert out == "\t  if (x) {\n\t      y();\n\n\t  }\n"


def test_echo_comments_the_directive(net, tmp_path):
    funcs = functions(tmp_path, LINES)
    out = render(net, tmp_path, "    $$JAFF CALL lines\n", funcs, echo=True)

    assert out == "    // $$JAFF CALL lines\n    a = 5;\n    b = 2;\n"


def test_empty_return_renders_nothing(net, tmp_path):
    funcs = functions(
        tmp_path,
        """
        from jaff.codegen import jaffgen_function

        @jaffgen_function
        def nothing(net):
            return ""
        """,
    )

    assert render(net, tmp_path, "x\n$$JAFF CALL nothing\ny\n", funcs) == "x\ny\n"


def test_replace_applies_to_every_returned_line(net, tmp_path):
    text = '  $$JAFF CALL lines { REPLACE { "(\\w) =" -> "\\1 :=" } }\n'
    out = render(net, tmp_path, text, functions(tmp_path, LINES))

    assert out == "  a := 5;\n  b := 2;\n"


def test_replace_from_a_define(net, tmp_path):
    text = '$$JAFF DEFINE r { REPLACE { ";" -> "" } }\n$$JAFF CALL lines { USE r }\n'
    out = render(net, tmp_path, text, functions(tmp_path, LINES))

    assert out == "a = 5\nb = 2\n"


def test_decorator_name_argument_sets_the_call_name(net, tmp_path):
    funcs = functions(
        tmp_path,
        """
        from jaff.codegen import jaffgen_function

        @jaffgen_function("greeting")
        def _hello(net):
            return "hello"
        """,
    )

    assert render(net, tmp_path, "$$JAFF CALL greeting\n", funcs) == "hello\n"
    assert funcs.names() == ["greeting"]


def test_undecorated_functions_are_not_callable(net, tmp_path):
    funcs = functions(tmp_path, LINES + "\n    def helper(net):\n        return 'x'\n")

    assert funcs.names() == ["lines"]


def test_decorator_returns_the_function_unchanged():
    def f(net):
        return "x"

    assert jaffgen_function(f) is f
    assert jaffgen_function("other")(f) is f


# --------------------------------------------------------------------------- #
# Loading                                                                      #
# --------------------------------------------------------------------------- #


def test_file_is_not_executed_without_a_call(net, tmp_path):
    funcs = functions(tmp_path, "raise RuntimeError('should not run')\n")

    assert render(net, tmp_path, "$$JAFF SUB nspec\n$nspec$\n$$JAFF END\n", funcs) == (
        "5\n"
    )


def test_file_is_executed_once_for_many_templates(net, tmp_path):
    counter = tmp_path / "count.txt"
    funcs = functions(
        tmp_path,
        f"""
        from pathlib import Path
        from jaff.codegen import jaffgen_function

        path = Path({str(counter)!r})
        path.write_text(path.read_text() + "x" if path.exists() else "x")

        @jaffgen_function
        def one(net):
            return "1"
        """,
    )
    for name in ("a.cpp", "b.cpp"):
        render(net, tmp_path, "$$JAFF CALL one\n$$JAFF CALL one\n", funcs, name=name)

    assert counter.read_text() == "x"


def test_loaders_do_not_share_functions(tmp_path):
    (tmp_path / "1").mkdir()
    (tmp_path / "2").mkdir()
    first = functions(tmp_path / "1", LINES)
    second = functions(
        tmp_path / "2",
        """
        from jaff.codegen import jaffgen_function

        @jaffgen_function
        def other(net):
            return ""
        """,
    )

    assert first.names() == ["lines"]
    assert second.names() == ["other"]


# --------------------------------------------------------------------------- #
# Errors                                                                       #
# --------------------------------------------------------------------------- #


def error(net, tmp_path, text, funcs) -> ParserError:
    with pytest.raises(ParserError) as info:
        render(net, tmp_path, text, funcs)

    return info.value


def test_unknown_name_lists_the_known_ones(net, tmp_path):
    e = error(net, tmp_path, "x\n$$JAFF CALL nope\n", functions(tmp_path, LINES))

    assert "CALL has no function nope (known: lines)" in e.message
    assert e.nline == 2


def test_call_without_a_functions_file(net, tmp_path):
    e = error(net, tmp_path, "$$JAFF CALL lines\n", None)

    assert FUNCTIONS in e.message


def test_unknown_name_is_found_before_anything_renders(net, tmp_path):
    funcs = functions(
        tmp_path,
        """
        from jaff.codegen import jaffgen_function

        @jaffgen_function
        def boom(net):
            raise AssertionError("rendered")
        """,
    )
    e = error(net, tmp_path, "$$JAFF CALL boom\n$$JAFF CALL nope\n", funcs)

    assert "nope" in e.message


def test_file_that_fails_to_load(net, tmp_path):
    funcs = functions(tmp_path, "raise RuntimeError('broken file')\n")
    e = error(net, tmp_path, "$$JAFF CALL lines\n", funcs)

    assert FUNCTIONS in e.message and "broken file" in e.message
    assert isinstance(e.__cause__, RuntimeError)


def test_duplicate_registered_name(net, tmp_path):
    funcs = functions(
        tmp_path,
        """
        from jaff.codegen import jaffgen_function

        @jaffgen_function
        def twice(net):
            return "1"

        @jaffgen_function("twice")
        def again(net):
            return "2"
        """,
    )
    e = error(net, tmp_path, "$$JAFF CALL twice\n", funcs)

    assert "twice" in e.message and "more than once" in e.message


def test_function_that_raises(net, tmp_path):
    funcs = functions(
        tmp_path,
        """
        from jaff.codegen import jaffgen_function

        @jaffgen_function
        def boom(net):
            raise KeyError("missing")
        """,
    )
    e = error(net, tmp_path, "$$JAFF CALL boom\n", funcs)

    assert "CALL failed: KeyError" in e.message
    assert isinstance(e.__cause__, KeyError)


def test_non_string_return(net, tmp_path):
    funcs = functions(
        tmp_path,
        """
        from jaff.codegen import jaffgen_function

        @jaffgen_function
        def number(net):
            return 3
        """,
    )
    e = error(net, tmp_path, "$$JAFF CALL number\n", funcs)

    assert "CALL number returned int, expected str" in e.message


def test_option_group_is_rejected(net, tmp_path):
    e = error(net, tmp_path, "$$JAFF CALL lines { CSE { use -> true } }\n", None)

    assert "CALL takes no option group CSE" in e.message


def test_call_has_no_body(net, tmp_path):
    e = error(
        net, tmp_path, "$$JAFF CALL lines\n$$JAFF END\n", functions(tmp_path, LINES)
    )

    assert "END without a directive" in e.message
