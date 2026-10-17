# ABOUTME: Every marimo notebook in examples/ runs headless and exits 0
# ABOUTME: A cell that raises (including a failed shell command) makes the run fail

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

import pytest

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"

NOTEBOOKS: List[str] = sorted(p.name for p in EXAMPLES.glob("*.py"))

# External programs a notebook shells out to. It is skipped if any is missing, and
# off Linux altogether: Windows/macOS toolchains differ (e.g. no plain gfortran,
# multi-config build/Debug/*.exe layouts).
REQUIRES: Dict[str, List[str]] = {
    "04_codegen_jaffgen.py": ["gfortran", "make"],
    "05_kokkos_cpp.py": ["cmake", "clang++"],
}


def _params() -> List[Any]:
    params = []
    for name in NOTEBOOKS:
        tools = REQUIRES.get(name, [])
        missing = [tool for tool in tools if shutil.which(tool) is None]
        if tools and not sys.platform.startswith("linux"):
            marks = pytest.mark.skip(reason="builds external code; Linux only")
        elif missing:
            marks = pytest.mark.skip(reason=f"needs {', '.join(missing)}")
        else:
            marks = ()
        params.append(pytest.param(name, marks=marks, id=name))

    return params


@pytest.mark.parametrize("name", _params())
def test_notebook_runs(name: str, tmp_path: Path):
    proc = subprocess.run(
        [sys.executable, str(EXAMPLES / name)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        env={**os.environ, "MPLBACKEND": "Agg"},
    )

    assert proc.returncode == 0, proc.stdout + proc.stderr


def test_requirements_name_existing_notebooks():
    """A renamed notebook must not silently lose its tool requirements."""
    stale = set(REQUIRES) - set(NOTEBOOKS)

    assert not stale
