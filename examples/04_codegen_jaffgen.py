# ABOUTME: Example notebook: generate Python and Fortran solvers with the jaffgen CLI,
# ABOUTME: then run them; every shell command must exit 0 (marimo, runs headless)

import marimo

__generated_with = "0.24.2"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # 04 · Generating solver code with `jaffgen`

    The previous notebooks built a symbolic model of a network. To actually
    evolve the chemistry, that model has to become source code for a solver:
    rate coefficients, reaction fluxes and the rate equations
    $\mathrm{d}n_i/\mathrm{d}t$, written in the language your simulation uses.

    `jaffgen` is JAFF's command-line code generator. It takes a network and a set
    of **templates** (ordinary source files with a few JAFF directives in them)
    and fills the directives in with code for that network. This notebook
    generates a Python solver and a Fortran solver for `demo2`, the two-reaction
    CO network from `03_comparing_networks.py`, and runs both.

    Every command below goes through a small `run()` helper (bottom of the
    notebook) that shows the command's output and **fails the cell if the exit
    code is not 0**, so a broken step stops the notebook instead of passing
    silently.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Template collections

    `jaffgen --template NAME` selects a predefined collection shipped with JAFF.
    A collection is one folder, `templates/NAME/` at the repository root,
    holding two kinds
    of file:

    - files with JAFF directives are **rendered** for your network (rates,
      fluxes, ODEs, species indices);
    - files without directives come out **as they are**: hand-written drivers,
      build files and solver libraries that do not depend on the network.

    The table lists every collection and its files:
    """)
    return


@app.cell
def _(TEMPLATES_DIR, mo):
    _rows = []
    for _gen in sorted(TEMPLATES_DIR.iterdir()):
        if not _gen.is_dir():
            continue
        # jaffgen.toml is the collection's config, not an output file.
        _files = sorted(
            p for p in _gen.iterdir() if p.is_file() and p.name != "jaffgen.toml"
        )
        _rendered = [p.name for p in _files if "$$JAFF" in p.read_text(errors="ignore")]
        _rows.append(
            {
                "collection": _gen.name,
                "rendered (JAFF directives)": ", ".join(_rendered),
                "as-is (no directives)": ", ".join(
                    p.name for p in _files if p.name not in _rendered
                ),
            }
        )
    mo.ui.table(_rows, selection=None, show_column_summaries=False)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## What a template looks like

    Here is `chemistry.py` from the `python_solve_ivp` collection. It is plain
    Python apart from the blocks between a `$$JAFF` directive and `$$JAFF END`:

    - `$$JAFF REPEAT idx, specie IN species` repeats the enclosed line once per
      species, binding `idx` to its index and `specie` to its name; `$idx$` and
      `$specie$` mark where those values are substituted;
    - `$$JAFF REPEAT ... IN parameters` lists the network's free parameters, so
      the generated functions read exactly the ones the network uses;
    - `$$JAFF SUB` substitutes single values, such as the number of species or
      the EOS energy;
    - the `rhses` and `jacobian` blocks print the right-hand side and Jacobian
      (with `THERMAL { MODE -> dtdt }`, including the gas temperature), sharing
      common subexpressions through the `cses` blocks.

    The language is taken from the file extension, so the same directives print
    Fortran in a `.f90` file.
    """)
    return


@app.cell
def _(TEMPLATES_DIR, mo):
    mo.md(
        "```python\n"
        + (TEMPLATES_DIR / "python_solve_ivp" / "chemistry.py").read_text()
        + "\n```"
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Generating the Python solver

    The network to generate for is `demo2`:
    """)
    return


@app.cell
def _(DEMO2, mo):
    mo.md("```text\n" + DEMO2.read_text() + "\n```")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    Generate the `python_solve_ivp` collection for it into a scratch directory.
    Each run writes into a new folder there, named
    `<network>_<template>_<date>-<time>`, so earlier runs are kept:
    """)
    return


@app.cell
def _(DEMO2, JAFFGEN, run, scratch):
    run(
        [
            JAFFGEN,
            "--network",
            DEMO2,
            "--template",
            "python_solve_ivp",
            "--outdir",
            scratch / "python_solve_ivp",
        ]
    )
    (py_out,) = (scratch / "python_solve_ivp").iterdir()
    return (py_out,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    `jaffgen` wrote one file per template. `chemistry.py` is rendered for the
    network: species map, parameter names, EOS energy, right-hand side and
    Jacobian. `solver.py` (`solve()`, SciPy's Radau with the analytic Jacobian)
    and `main.py` (the driver) are copied as they are. `jaffgen` also writes
    `config.toml`, a record of the settings used for the run, which can be
    passed back with `--config` to repeat it.
    """)
    return


@app.cell
def _(file_table, py_out):
    file_table(py_out)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    The rendered `chemistry.py` shows the directives at work. Each block now
    holds the code for `demo2`, and the directives themselves are kept as
    comments (`--echo`, on by default) so you can see where generated code came
    from:
    """)
    return


@app.cell
def _(mo, py_out):
    mo.md("```python\n" + (py_out / "chemistry.py").read_text() + "\n```")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ### Running it

    The generated directory is a self-contained program that needs only NumPy
    and SciPy. The driver starts from the inputs block at the top of `main.py`;
    as generated, every density is 0 and every parameter a placeholder 1, so
    edit that block for a real run. Here it just shows the program runs and
    saves `solution.dat`.
    """)
    return


@app.cell
def _(py_out, run, sys):
    run([sys.executable, "main.py"], cwd=py_out)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    The generated modules can also be imported and used directly. The next cell
    calls `solve()` from the generated `solver.py`, starting from pure CO at
    100 K, with densities given by species name and the one parameter `demo2`
    uses, $A_V$. CO is ionised by starlight (`CO -> CO+ + e-`, attenuated
    by $A_V = 1$) while CO⁺ recombines with electrons (`CO+ + e- -> CO`). The
    two balance once ionisation and recombination rates are equal, and CO⁺ and
    e⁻ settle at the same equilibrium density.
    """)
    return


@app.cell
def _(PLOT_SIZE, import_generated, np, plt, py_out):
    _gen = import_generated(py_out, ["chemistry", "solver"])
    _chemistry, _solver = _gen["chemistry"], _gen["solver"]

    _year = 365.0 * 24 * 3600
    _t_eval = np.geomspace(1e2, 1e6, 100) * _year
    _result = _solver.solve(
        {"CO": 1e2},  # cm^-3
        100.0,  # K
        _t_eval[-1],
        {"av": 1.0},  # every name in _chemistry.PARAMETERS
        t_eval=_t_eval,
    )

    _fig, _ax = plt.subplots(figsize=PLOT_SIZE)
    for _name in ("CO", "CO+", "e-"):
        _ax.loglog(_result.t / _year, _result.n[_name], label=_name)
    _ax.set_xlabel("time (yr)")
    _ax.set_ylabel("number density (cm$^{-3}$)")
    _ax.legend()
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Generating the Fortran solver

    Switching language means switching collection; the network stays the same.
    `fortran_dlsodes` renders `chemistry.f90` (species, parameters, EOS energy,
    right-hand side and Jacobian) and `user_functions.f90` (stubs for functions
    the network calls without defining; `demo2` has none), and copies `solver.f90`
    (DLSODES with the analytic Jacobian), the driver `main.f90`, a `Makefile` and
    the ODEPACK sources.
    """)
    return


@app.cell
def _(DEMO2, JAFFGEN, run, scratch):
    run(
        [
            JAFFGEN,
            "--network",
            DEMO2,
            "--template",
            "fortran_dlsodes",
            "--outdir",
            scratch / "fortran_dlsodes",
        ]
    )
    (f90_out,) = (scratch / "fortran_dlsodes").iterdir()
    return (f90_out,)


@app.cell
def _(f90_out, file_table):
    file_table(f90_out)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    The same directives now print Fortran: double-precision `d` exponents and
    `merge` for conditions. The parameters are module variables the driver sets,
    and `implicit double precision (x)` types the common-subexpression
    temporaries:
    """)
    return


@app.cell
def _(f90_out, mo):
    mo.md("```fortran\n" + (f90_out / "chemistry.f90").read_text() + "\n```")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ### Building and running it

    The driver starts from the inputs block at the top of `main.f90`: one line
    per species (all 0 as generated), the parameters (placeholder 1), the
    temperature, the end time and whether to evolve the temperature. Here we start
    from 100 cm⁻³ of CO, as in the Python run, by editing that block. `make`
    compiles everything with `gfortran` (or Intel's `ifx` when it is installed)
    into an executable called `test`, which integrates for $10^6$ years, prints
    the final state and saves the evolution to `solution.dat`; it exits with a
    nonzero code if DLSODES reports a failure.
    """)
    return


@app.cell
def _(f90_out, re, run):
    _main = f90_out / "main.f90"
    _main.write_text(
        re.sub(r"= 0d0(  ! CO)$", r"= 1d2\1", _main.read_text(), flags=re.MULTILINE)
    )
    run(["make"], cwd=f90_out)
    return


@app.cell
def _(f90_out, run):
    run(["./test"], cwd=f90_out)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    After $10^6$ years CO⁺ and e⁻ have reached the same equilibrium as in the
    Python run, with the temperature held at 100 K.

    **Functions the network does not define:** photo-reactions written as
    `PHOTO` lines (like those in `demo1`) are emitted as `photorates(...)` calls,
    as are table interpolators. Neither collection implements them: the Python
    `chemistry.py` looks them up at call time, and the Fortran `user_functions.f90`
    has a stub for each to fill in. This notebook uses `demo2`, which needs none.

    **Next:** `05_kokkos_cpp.py` generates a C++ solver for a larger network and
    builds it with CMake.
    """)
    return


@app.cell(hide_code=True)
def _():
    import importlib
    import logging
    import os
    import re
    import subprocess
    import sys
    import tempfile
    from pathlib import Path
    from typing import Dict, List

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np

    from jaff.config import TEMPLATES_DIR

    # JAFF resets its logger to INFO on every JaffLogger(); disabling WARNING and
    # below is the only switch that sticks.
    logging.disable(logging.WARNING)

    DEMO2 = mo.notebook_dir().parent / "networks" / "demos" / "demo2.jet"
    # The jaffgen console script installed next to this Python interpreter.
    JAFFGEN = Path(sys.executable).parent / "jaffgen"
    PLOT_SIZE = (5.5, 3.5)

    # Returned from this cell so it stays alive: a TemporaryDirectory deletes its
    # folder when garbage-collected, and _-prefixed names are dropped after the cell.
    scratch_dir = tempfile.TemporaryDirectory()
    scratch = Path(scratch_dir.name)

    def run(cmd: List, cwd: Path | None = None, env: Dict | None = None):
        """Run *cmd*, show its output, and fail the cell on a nonzero exit code."""
        cmd = [str(c) for c in cmd]
        proc = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True)
        # Show installed executables by name (jaffgen, python), not full path.
        program = Path(cmd[0]).name if Path(cmd[0]).is_absolute() else cmd[0]
        shown = " ".join([program, *cmd[1:]])
        assert proc.returncode == 0, (
            f"`{shown}` exited with {proc.returncode}:\n{proc.stdout}{proc.stderr}"
        )

        return mo.md(f"```console\n$ {shown}\n{proc.stdout}{proc.stderr}\n```")

    def file_table(folder: Path):
        """List the files in *folder* with their sizes."""
        return mo.ui.table(
            [
                {"file": p.name, "bytes": p.stat().st_size}
                for p in sorted(folder.iterdir())
                if p.is_file()
            ],
            selection=None,
            show_column_summaries=False,
        )

    def import_generated(folder: Path, names: List[str]) -> Dict:
        """Import generated modules from *folder*, replacing stale copies."""
        sys.path.insert(0, str(folder))
        try:
            for name in names:
                sys.modules.pop(name, None)

            return {name: importlib.import_module(name) for name in names}
        finally:
            sys.path.remove(str(folder))

    return (
        DEMO2,
        JAFFGEN,
        PLOT_SIZE,
        TEMPLATES_DIR,
        file_table,
        import_generated,
        mo,
        np,
        os,
        plt,
        re,
        run,
        scratch,
        scratch_dir,
        sys,
    )


if __name__ == "__main__":
    app.run()
