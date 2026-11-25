# ABOUTME: Runs every docstring Examples section in src/jaff (one test per docstring), in
# ABOUTME: source order per module, in a scratch folder holding the documented input files

from __future__ import annotations

import builtins
import doctest
import importlib
import io
import logging
import pkgutil
import re

import matplotlib
import pytest

import jaff
from tests._paths import REPO_ROOT

matplotlib.use("Agg")

#: The real print: tests/conftest.py silences it, but examples print their output.
_PRINT = builtins.print

#: Developer scripts, not documented API.
_SKIP_PACKAGES = ("jaff._utils",)

#: Input files the examples describe in prose ("Given ``grid.csv``: ...").
_INPUT_FILES = {
    "grid.csv": "T,n,k\n10,1,0.1\n10,2,0.2\n20,1,0.3\n20,2,0.4\n",
    "nspec.cpp": "$$JAFF SUB nspec\nconst int n = $nspec$;\n$$JAFF END\n",
}

#: Repo folders the examples name by relative path (``networks/demos/...``).
_REPO_LINKS = ("networks",)

_FLAGS = doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE

#: A number in example output: an integer, a decimal or one in exponent notation.
_NUMBER = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")

#: Relative tolerance of numbers in example output: the last digits of a float differ
#: between platforms' maths libraries (macOS vs glibc).  No absolute tolerance, so tiny
#: values (cross sections ~1e-19) are compared as strictly as large ones.
_REL = 1e-12


class _ApproxChecker(doctest.OutputChecker):
    """Doctest output comparison with the numbers compared by ``pytest.approx``.

    Output that matches as text passes as usual.  Otherwise the text between the
    numbers must match exactly and every number must equal the expected one within
    :data:`_REL`, so documented output stays the real output of one platform.
    """

    def check_output(self, want: str, got: str, optionflags: int) -> bool:
        if super().check_output(want, got, optionflags):
            return True

        wanted, found = _NUMBER.findall(want), _NUMBER.findall(got)
        if len(wanted) != len(found):
            return False
        if any(
            float(g) != pytest.approx(float(w), rel=_REL, abs=0)
            for w, g in zip(wanted, found)
        ):
            return False

        # The numbers agree: compare the text with the expected numbers put in.
        expected = iter(wanted)

        return super().check_output(
            want, _NUMBER.sub(lambda _: next(expected), got), optionflags
        )


#: Examples known to fail, with why; empty this list rather than add to it.
XFAIL: dict[str, str] = {}


def _modules() -> list[str]:
    names = ["jaff"]
    for info in pkgutil.walk_packages(jaff.__path__, "jaff."):
        if not info.name.startswith(_SKIP_PACKAGES):
            names.append(info.name)

    return sorted(names)


def _examples() -> list[tuple[str, doctest.DocTest]]:
    finder = doctest.DocTestFinder(exclude_empty=True)
    found = []
    for name in _modules():
        module = importlib.import_module(name)
        tests = [t for t in finder.find(module) if t.examples]
        # Source order: a module's later examples may read files earlier ones wrote.
        found += [(name, t) for t in sorted(tests, key=lambda t: t.lineno or 0)]

    return found


_EXAMPLES = _examples()


@pytest.fixture(scope="session")
def _module_dirs(tmp_path_factory):
    """One scratch folder per module, shared by that module's examples."""
    dirs = {}

    def get(module: str):
        if module not in dirs:
            path = tmp_path_factory.mktemp(module.replace(".", "_"))
            for name, text in _INPUT_FILES.items():
                (path / name).write_text(text, encoding="utf-8")
            for name in _REPO_LINKS:
                (path / name).symlink_to(REPO_ROOT / name)
            dirs[module] = path

        return dirs[module]

    return get


def _params():
    for module, test in _EXAMPLES:
        reason = XFAIL.get(test.name)
        marks = [pytest.mark.xfail(reason=reason, strict=True)] if reason else []
        yield pytest.param(module, test, id=test.name, marks=marks)


@pytest.mark.parametrize(("module", "test"), list(_params()))
def test_docstring_example(module, test, _module_dirs, monkeypatch):
    monkeypatch.chdir(_module_dirs(module))
    monkeypatch.setattr(builtins, "print", _PRINT)
    # The banner and log lines are not part of the documented output.
    monkeypatch.setattr("jaff.core.network.network.print_motd", lambda *a, **k: None)
    logging.disable(logging.CRITICAL)
    out = io.StringIO()
    try:
        runner = doctest.DocTestRunner(checker=_ApproxChecker(), optionflags=_FLAGS)
        runner.run(test, out=out.write, clear_globs=True)
    finally:
        logging.disable(logging.NOTSET)
        matplotlib.pyplot.close("all")

    if runner.failures:
        pytest.fail(out.getvalue(), pytrace=False)


def test_examples_are_found():
    names = {test.name for _, test in _EXAMPLES}

    assert "jaff.drivers.csv.csv_to_df" in names
    assert len(names) > 100
    assert set(XFAIL) <= names, "XFAIL names an example that no longer exists"


class TestApproxChecker:
    @pytest.mark.parametrize(
        ("want", "got"),
        [
            ("1.2403886126668026e-19\n", "1.240388612666803e-19\n"),
            ("(1.95e-18, 15.735906475383057)\n", "(1.95e-18, 15.735906475383056)\n"),
            ("5.85154082319921e-8*photden[1]\n", "5.851540823199211e-8*photden[1]\n"),
        ],
    )
    def test_last_digit_differences_pass(self, want, got):
        assert _ApproxChecker().check_output(want, got, _FLAGS)

    @pytest.mark.parametrize(
        ("want", "got"),
        [
            ("1.2403886126668026e-19\n", "1.2403886e-19\n"),  # beyond the tolerance
            ("1.0e-19\n", "5.0e-19\n"),  # tiny values have no absolute tolerance
            ("2.0*photden[1]\n", "2.0*photden[2]\n"),  # an index is a number too
            ("2.0*photden[1]\n", "2.0*nden[1]\n"),  # text must match exactly
            ("(1.0, 2.0)\n", "(1.0, 2.0, 3.0)\n"),
        ],
    )
    def test_other_differences_fail(self, want, got):
        assert not _ApproxChecker().check_output(want, got, _FLAGS)
