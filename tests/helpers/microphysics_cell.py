# ABOUTME: Shared helpers for tests of generated microphysics networks: build the one-zone
# ABOUTME: Microphysics test (tests/fixtures/microphysics_cell) against a run folder, run it

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

from tests._paths import FIXTURES_DIR

CELL = FIXTURES_DIR / "microphysics_cell"

#: Environment variables naming the Microphysics and AMReX checkouts
HOMES = ("MICROPHYSICS_HOME", "AMREX_HOME")


def available() -> bool:
    """Whether Microphysics and AMReX are named in the environment and make exists."""
    return (
        all(os.environ.get(name) for name in HOMES) and shutil.which("make") is not None
    )


def build(network: Path, workdir: Path) -> subprocess.CompletedProcess:
    """Copy the one-zone test into *workdir* and build it with the generated *network*.

    The network folder is a jaffgen ``microphysics`` run; its ``eos/`` is the EOS.
    """
    shutil.copytree(CELL, workdir, dirs_exist_ok=True)
    command = [
        "make",
        f"-j{os.cpu_count() or 1}",
        *(f"{name}={os.environ[name]}" for name in HOMES),
        f"NETWORK_CUSTOM_PATH={network}",
        f"EOS_PATH={network / 'eos'}",
    ]

    return subprocess.run(command, cwd=workdir, capture_output=True, text=True)


def run(
    workdir: Path,
    densities: dict[str, float],
    parameters: dict[str, float],
    *settings: str,
) -> subprocess.CompletedProcess:
    """Run the built test from *densities* with the network *parameters*.

    *settings* are extra ``key=value`` inputs (``unit_test.tmax=1e10``, ...).
    """
    initial = workdir / "initial.txt"
    initial.write_text(
        "".join(f"{name} {value!r}\n" for name, value in densities.items()),
        encoding="utf-8",
    )
    (exe,) = workdir.glob("main*.ex")
    command = [
        str(exe),
        "inputs",
        *(f"network.{name}={value!r}" for name, value in parameters.items()),
        *settings,
    ]

    return subprocess.run(command, cwd=workdir, capture_output=True, text=True)


def final_state(stdout: str) -> tuple[bool, float, dict[str, float]]:
    """(success, T, densities by name) from the test's output."""
    values = dict(line.split() for line in stdout.splitlines() if len(line.split()) == 2)
    success = values.pop("success") == "1"
    temperature = float(values.pop("T"))

    return success, temperature, {name: float(v) for name, v in values.items()}
