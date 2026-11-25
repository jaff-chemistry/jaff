# ABOUTME: jaffgen finds jaffgen_template_functions.py among all its input files,
# ABOUTME: never writes it (nor jaffgen.toml) and shares one loader across the run

from pathlib import Path
from types import SimpleNamespace

import pytest

from jaff.cli.jaffgen._engine import JaffGen, run_config
from tests._paths import FIXTURES_DIR
from tests.helpers.jaffgen_run import run_folder, run_jaffgen

NETWORK = FIXTURES_DIR / "react_cie_hepp.jet"
FUNCTIONS = "jaffgen_template_functions.py"

SOURCE = """\
from jaff.codegen import jaffgen_function

@jaffgen_function
def nspec(net):
    return f"n = {len(net.species)}"
"""


def args(tmp_path: Path, **kw) -> SimpleNamespace:
    return SimpleNamespace(network=str(NETWORK), outdir=str(tmp_path / "out"), **kw)


def write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)

    return path


def test_functions_file_in_the_input_directory(tmp_path):
    indir = tmp_path / "in"
    write(indir / "sub" / FUNCTIONS, SOURCE)
    write(indir / "a.py", "$$JAFF CALL nspec\n")
    write(indir / "b.py", "x = 1\n    $$JAFF CALL nspec\n")

    gen = run_jaffgen(args(tmp_path, indir=str(indir), echo=False))
    out = run_folder(tmp_path / "out")

    assert (out / "a.py").read_text() == "n = 5\n"
    assert (out / "b.py").read_text() == "x = 1\n    n = 5\n"
    assert gen.functions is not None
    assert gen.functions.path == (indir / "sub" / FUNCTIONS).resolve()


def test_functions_file_from_files(tmp_path):
    funcs = write(tmp_path / "src" / FUNCTIONS, SOURCE)
    template = write(tmp_path / "src" / "a.py", "$$JAFF CALL nspec\n")

    run_jaffgen(args(tmp_path, files=[str(funcs), str(template)], echo=False))

    assert (run_folder(tmp_path / "out") / "a.py").read_text() == "n = 5\n"


def test_functions_file_and_jaffgen_toml_are_not_written(tmp_path):
    indir = tmp_path / "in"
    write(indir / FUNCTIONS, SOURCE)
    # In a subfolder: the run folder's own jaffgen.toml is jaffgen's config record.
    write(indir / "sub" / "jaffgen.toml", "[jaffgen]\n")
    write(indir / "a.py", "x = 1\n")

    run_jaffgen(args(tmp_path, indir=str(indir)))
    out = run_folder(tmp_path / "out")
    written = {p.relative_to(out).as_posix() for p in out.rglob("*") if p.is_file()}

    assert FUNCTIONS not in written
    assert "sub/jaffgen.toml" not in written
    assert "a.py" in written


def test_no_functions_file(tmp_path):
    indir = tmp_path / "in"
    write(indir / "a.py", "x = 1\n")

    gen = JaffGen(run_config(args(tmp_path, indir=str(indir))))

    assert gen.functions is None


def test_more_than_one_functions_file_is_an_error(tmp_path):
    indir = tmp_path / "in"
    first = write(indir / FUNCTIONS, SOURCE)
    second = write(tmp_path / "extra" / FUNCTIONS, SOURCE)
    write(indir / "a.py", "x = 1\n")

    with pytest.raises(ValueError, match="more than one") as info:
        JaffGen(run_config(args(tmp_path, indir=str(indir), files=[str(second)])))

    assert str(first.resolve()) in str(info.value)
    assert str(second.resolve()) in str(info.value)


def test_functions_file_is_not_run_at_collection(tmp_path):
    indir = tmp_path / "in"
    write(indir / FUNCTIONS, "raise RuntimeError('ran')\n")
    write(indir / "a.py", "x = 1\n")

    run_jaffgen(args(tmp_path, indir=str(indir)))

    assert (run_folder(tmp_path / "out") / "a.py").read_text() == "x = 1\n"
