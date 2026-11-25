# ABOUTME: Shared test helper: run jaffgen from a CLI-like namespace of options,
# ABOUTME: as the typer command does (run_config, JaffGen(config), run())

from pathlib import Path
from types import SimpleNamespace

from jaff.cli import JaffGen
from jaff.cli.jaffgen._engine import run_config


def run_jaffgen(args: SimpleNamespace) -> JaffGen:
    """Load the run's config from *args*, build the run, write its output; return it."""
    gen = JaffGen(run_config(args))
    gen.run()

    return gen


def run_folder(outdir: Path) -> Path:
    """The one run folder jaffgen wrote in *outdir* (``<network>_<template>_<time>``)."""
    (folder,) = [p for p in outdir.iterdir() if p.is_dir()]

    return folder
