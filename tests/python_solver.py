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


def reference(
    network: str,
    outdir: Path,
    densities: dict[str, float],
    temperature: float,
    t_end: float,
    parameters: dict[str, float],
) -> tuple[float, dict[str, float]]:
    """(T, densities) at *t_end* from the Python template (Radau, rtol 1e-12).

    The referee for other templates' runs: the gas temperature evolves.
    """
    modules = load(generate(network, outdir))
    result = modules.solver.solve(
        densities,
        temperature,
        t_end,
        parameters,
        evolve_temperature=True,
        t_eval=[t_end],
        rtol=1e-12,
        atol=1e-20,
    )

    return float(result.temperature[-1]), {k: float(v[-1]) for k, v in result.n.items()}
