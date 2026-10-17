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
    A collection has two parts:

    - files under `templates/generator/NAME/` contain JAFF directives and are
      **rendered** for your network (rates, fluxes, ODEs, species indices);
    - files under `templates/preprocessor/NAME/` that are not also in
      `generator/` are **added as they are**: hand-written drivers, build files
      and solver libraries that do not depend on the network.

    The table lists every collection and its files:
    """)
    return


@app.cell
def _(TEMPLATES_DIR, mo):
    _rows = []
    for _gen in sorted((TEMPLATES_DIR / "generator").iterdir()):
        if not _gen.is_dir():
            continue
        _rendered = sorted(p.name for p in _gen.iterdir() if p.is_file())
        _static_dir = TEMPLATES_DIR / "preprocessor" / _gen.name
        _static = (
            sorted(p.name for p in _static_dir.iterdir() if p.is_file())
            if _static_dir.is_dir()
            else []
        )
        _rows.append(
            {
                "collection": _gen.name,
                "rendered (generator/)": ", ".join(_rendered),
                "added as-is (preprocessor/)": ", ".join(
                    f for f in _static if f not in _rendered
                ),
            }
        )
    mo.ui.table(_rows, selection=None, show_column_summaries=False)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## What a template looks like

    Here is `rates.py` from the `python_solve_ivp` collection. It is plain Python
    apart from the block between `$$JAFF REPEAT` and `$$JAFF END`:

    - `$$JAFF REPEAT idx, rate IN rates` repeats the enclosed lines once per
      reaction, binding `idx` to the reaction's index and `rate` to its rate
      coefficient printed as Python;
    - `$idx$` and `$rate$` mark where those values are substituted.

    Other directives work the same way: `$$JAFF SUB` substitutes single values
    (such as the number of species), and the iterable after `IN` can be species,
    fluxes, ODE right-hand sides, Jacobian entries and more. The language is
    taken from the file extension, so the same directive prints Fortran in a
    `.f90` file.
    """)
    return


@app.cell
def _(TEMPLATES_DIR, mo):
    mo.md(
        "```python\n"
        + (TEMPLATES_DIR / "generator" / "python_solve_ivp" / "rates.py").read_text()
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
    Generate the `python_solve_ivp` collection for it into a scratch directory:
    """)
    return


@app.cell
def _(DEMO2, JAFFGEN, run, scratch):
    py_out = scratch / "python_solve_ivp"
    run(
        [
            JAFFGEN,
            "--network",
            DEMO2,
            "--template",
            "python_solve_ivp",
            "--outdir",
            py_out,
        ]
    )
    return (py_out,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    `jaffgen` wrote one file per template. `commons.py` holds the species
    indices and sizes, `rates.py` the rate coefficients, `fluxes.py` the
    reaction fluxes and `ode.py` the right-hand side of the rate equations;
    `main.py` is the hand-written driver that integrates them with SciPy's
    `solve_ivp`.
    """)
    return


@app.cell
def _(file_table, py_out):
    file_table(py_out)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    The rendered `rates.py` shows the directive at work. The block now holds one
    line per reaction of `demo2`, and the directive itself is kept as a comment
    (`--echo`, on by default) so you can see where generated code came from:
    """)
    return


@app.cell
def _(mo, py_out):
    mo.md("```python\n" + (py_out / "rates.py").read_text() + "\n```")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ### Running it

    The generated directory is a self-contained program. Running its driver
    integrates the chemistry from pure CO at 100 K for $10^6$ years. The driver
    calls `plt.show()`, so it is run with the non-interactive `Agg` backend;
    what matters here is that it exits 0.
    """)
    return


@app.cell
def _(os, py_out, run, sys):
    run(
        [sys.executable, "main.py"],
        cwd=py_out,
        env={**os.environ, "MPLBACKEND": "Agg"},
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    The generated modules can also be imported and used directly. The next cell
    imports `get_ode` from the generated `ode.py` and integrates it the same way
    the driver does. CO is ionised by starlight (`CO -> CO+ + e-`, attenuated
    by $A_V = 1$) while CO⁺ recombines with electrons (`CO+ + e- -> CO`). The
    two balance once ionisation and recombination rates are equal, and CO⁺ and
    e⁻ settle at the same equilibrium density.
    """)
    return


@app.cell
def _(PLOT_SIZE, import_generated, np, plt, py_out):
    from scipy.integrate import solve_ivp

    _gen = import_generated(py_out, ["commons", "rates", "fluxes", "ode"])
    _commons, _ode = _gen["commons"], _gen["ode"]

    _y0 = np.zeros(_commons.nvars)
    _y0[_commons.idx_co] = 1e2  # cm^-3
    _y0[_commons.idx_tgas] = 1e2  # K
    _year = 365.0 * 24 * 3600
    _t_eval = np.logspace(2, 6, 100) * _year
    _crate, _av = 1e-17, 1.0

    _sol = solve_ivp(
        lambda t, y: _ode.get_ode(y, y[_commons.idx_tgas], _crate, _av),
        (0, _t_eval[-1]),
        _y0,
        t_eval=_t_eval,
        method="LSODA",
    )
    assert _sol.success, _sol.message

    _fig, _ax = plt.subplots(figsize=PLOT_SIZE)
    for _name, _idx in [
        ("CO", _commons.idx_co),
        ("CO+", _commons.idx_coj),
        ("e-", _commons.idx_e),
    ]:
        _ax.loglog(_sol.t / _year, _sol.y[_idx], label=_name)
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
    `fortran_dlsodes` produces Fortran modules for the same rates, fluxes and
    ODEs, plus a driver (`main.f90`), a `Makefile` and the ODEPACK sources of the
    stiff solver DLSODES.
    """)
    return


@app.cell
def _(DEMO2, JAFFGEN, run, scratch):
    f90_out = scratch / "fortran_dlsodes"
    run(
        [
            JAFFGEN,
            "--network",
            DEMO2,
            "--template",
            "fortran_dlsodes",
            "--outdir",
            f90_out,
        ]
    )
    return (f90_out,)


@app.cell
def _(f90_out, file_table):
    file_table(f90_out)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    The same `rates` directive now prints Fortran: double-precision `d`
    exponents, `exp` and `**`, and 1-based indices. Fortran is case-insensitive,
    so species indices are generated in lower case (`idx_co`, `idx_coj`) and the
    electron is mapped to `idx_e`:
    """)
    return


@app.cell
def _(f90_out, mo):
    mo.md(
        "```fortran\n"
        + (f90_out / "commons.f90").read_text()
        + "\n"
        + (f90_out / "reactions.f90").read_text()
        + "\n```"
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ### Building and running it

    `make` compiles everything with `gfortran` (or Intel's `ifx` when it is
    installed) into an executable called `test`. The driver integrates for one
    second from the same initial state and prints the final time, temperature
    and densities; it exits with a nonzero code if DLSODES reports a failure.
    """)
    return


@app.cell
def _(f90_out, run):
    run(["make"], cwd=f90_out)
    return


@app.cell
def _(f90_out, run):
    run(["./test"], cwd=f90_out)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    After one second almost nothing has happened yet: each CO ionisation makes
    one CO⁺ and one electron, so both equal
    $k\,n_\mathrm{CO}\,t = 10^{-12}e^{-1.222}\times100\times1 \approx 2.95\times10^{-11}\ \mathrm{cm^{-3}}$,
    while CO is still at $100\ \mathrm{cm^{-3}}$ and the temperature is fixed at
    100 K.

    **A current limitation:** both collections expect rates that can be
    evaluated directly. Photo-reactions written as `PHOTO` lines (like those in
    `demo1`) are emitted as `photorates(...)` calls, which neither collection
    defines yet, so this notebook uses `demo2`.

    **Next:** `05_kokkos_cpp.py` generates a C++ solver for a larger network and
    builds it with CMake.
    """)
    return


@app.cell(hide_code=True)
def _():
    import importlib
    import logging
    import os
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
        run,
        scratch,
        scratch_dir,
        sys,
    )


if __name__ == "__main__":
    app.run()
