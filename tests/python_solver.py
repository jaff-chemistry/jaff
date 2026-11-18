# ABOUTME: Shared test helper: generate the python_solve_ivp template for a network
# ABOUTME: and import its chemistry and solver modules afresh
"""Generate and import the standalone Python solver (``python_solve_ivp``)."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from types import SimpleNamespace

from tests.jaffgen_run import run_jaffgen

MODULES = ("chemistry", "solver")


def generate(network: str, outdir: Path, template: str = "python_solve_ivp") -> Path:
    """Run jaffgen with *template* (default ``python_solve_ivp``); return the run folder."""
    args = SimpleNamespace(
        network=network,
        config=None,
        label=None,
        funcfile=None,
        duplicate_policy=None,
        expand_nuclei=None,
        strict=None,
        network_config=None,
        outdir=str(outdir),
        indir=None,
        files=None,
        template=template,
    )

    return run_jaffgen(args).output_dir


def load(folder: Path) -> SimpleNamespace:
    """Import the generated ``chemistry`` and ``solver`` modules of *folder* afresh.

    Both are dropped from ``sys.modules`` before and after, so solvers generated
    for different networks can be loaded in one test session.
    """
    for name in MODULES:
        sys.modules.pop(name, None)
    sys.path.insert(0, str(folder))
    try:
        chemistry = importlib.import_module("chemistry")
        solver = importlib.import_module("solver")
    finally:
        sys.path.remove(str(folder))
        for name in MODULES:
            sys.modules.pop(name, None)

    return SimpleNamespace(chemistry=chemistry, solver=solver)
