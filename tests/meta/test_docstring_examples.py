# ABOUTME: Runs every docstring Examples section in src/jaff (one test per docstring), in
# ABOUTME: source order per module, in a scratch folder holding the documented input files

from __future__ import annotations

import builtins
import doctest
import importlib
import io
import logging
import pkgutil

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
        runner = doctest.DocTestRunner(optionflags=_FLAGS)
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
