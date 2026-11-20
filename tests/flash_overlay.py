# ABOUTME: Shared helpers for tests of generated FLASH networks: a temporary overlay of the
# ABOUTME: FLASH checkout with the generated units added, its build, and one-zone runs

"""Build and run a jaffgen ``flash`` run folder without touching the FLASH checkout.

FLASH's ``setup`` takes its home from the path of ``bin/setup.py`` (``abspath``,
symlinks kept), so an overlay folder works as a FLASH home: ``bin/`` holds
symlinks to FLASH's scripts, ``lib/``, ``sites/`` and ``tools/`` are symlinks,
and ``source/`` is a symlink copy of FLASH's source with the generated folders
copied in (replacing any of the same name) and the network added to the
KromeChemistry composition's EXCLUSIVE list.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from collections.abc import Mapping
from pathlib import Path

import numpy as np

COMPOSITION_CONFIG = Path("source/Simulation/SimulationComposition/KromeChemistry/Config")


def available() -> bool:
    """Whether ``FLASH_HOME`` names a FLASH checkout and its build tools exist."""
    home = os.environ.get("FLASH_HOME")
    tools = ("make", "mpif90", "cp", "python3")

    return bool(home) and Path(home, "setup").is_file() and all(map(shutil.which, tools))


def overlay(run: Path, dest: Path, label: str) -> Path:
    """A FLASH home in *dest*: the checkout in ``FLASH_HOME`` plus *run*'s units."""
    home = Path(os.environ["FLASH_HOME"])
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "bin").mkdir()
    for entry in (home / "bin").iterdir():
        (dest / "bin" / entry.name).symlink_to(entry)
    for name in ("lib", "sites", "tools"):
        if (home / name).exists():
            (dest / name).symlink_to(home / name)
    shutil.copy(home / "setup", dest / "setup")
    subprocess.run(["cp", "-rs", str(home / "source"), str(dest / "source")], check=True)

    # each generated folder replaces the overlay's folder of the same path
    generated = sorted({f.parent for f in (run / "source").rglob("*") if f.is_file()})
    for folder in generated:
        target = dest / folder.relative_to(run)
        if target.exists():
            shutil.rmtree(target)
        shutil.copytree(folder, target)
    _register(dest / COMPOSITION_CONFIG, label)

    return dest


def _register(config: Path, label: str) -> None:
    """Add *label* to the EXCLUSIVE list of the composition Config (a real file)."""
    lines = config.read_text().splitlines()
    config.unlink()
    lines = [
        f"{line} {label}"
        if line.startswith("EXCLUSIVE") and label not in line.split()
        else line
        for line in lines
    ]
    config.write_text("\n".join(lines) + "\n")


def build(home: Path, label: str) -> subprocess.CompletedProcess:
    """Set up the one-zone simulation *label* in ``home/obj`` and build ``flash4``.

    ``+eos_multigamma`` makes the hydro unit take the multigamma EOS (the species'
    adiabatic indices) instead of its default single-gamma one.
    """
    setup = subprocess.run(
        ["./setup", label, "-objdir=obj", "-auto", "-1d", "-nxb=1", "-maxblocks=1"]
        + ["+ug", "+eos_multigamma"],
        cwd=home,
        capture_output=True,
        text=True,
    )
    if setup.returncode != 0:
        return setup

    return subprocess.run(
        ["make", f"-j{os.cpu_count() or 1}"],
        cwd=home / "obj",
        capture_output=True,
        text=True,
    )


def run(
    obj: Path,
    densities: dict[str, float],
    settings: Mapping[str, object],
    workdir: Path,
) -> subprocess.CompletedProcess:
    """Run ``flash4`` in *workdir* from *densities* with flash.par *settings* changed.

    *settings* maps runtime parameter names to values (strings are quoted, bools
    written as Fortran logicals); the rest of the generated flash.par stays.
    """
    workdir.mkdir(parents=True, exist_ok=True)
    (workdir / "initial.txt").write_text(
        "".join(f"{name} {value!r}\n" for name, value in densities.items())
    )
    par = (obj / "flash.par").read_text().splitlines()
    keys = {k.lower() for k in settings}
    lines = [line for line in par if line.split("=")[0].strip().lower() not in keys]
    lines += [f"{key} = {_par_value(value)}" for key, value in settings.items()]
    (workdir / "flash.par").write_text("\n".join(lines) + "\n")

    return subprocess.run(
        [str(obj / "flash4")], cwd=workdir, capture_output=True, text=True
    )


def _par_value(value: object) -> str:
    if isinstance(value, bool):
        return ".true." if value else ".false."
    if isinstance(value, str):
        return f'"{value}"'

    return repr(value)


def final_state(workdir: Path) -> tuple[float, float, dict[str, float]]:
    """(t, T, densities by name) of the last line of ``solution.dat``."""
    lines = (workdir / "solution.dat").read_text().splitlines()
    names = lines[0].split()[3:]
    last = np.array([float(v) for v in lines[-1].split()])

    return float(last[0]), float(last[1]), dict(zip(names, map(float, last[2:])))
