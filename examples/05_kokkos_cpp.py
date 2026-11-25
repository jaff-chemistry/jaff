# ABOUTME: Example notebook: generate a C++ chemistry solver for the COthin network,
# ABOUTME: build it with CMake and run it; every command must exit 0 (marimo)

import marimo

__generated_with = "0.24.2"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    # 05 · A C++ solver built with CMake

    `04_codegen_jaffgen.py` generated small Python and Fortran solvers. This
    notebook does the same for C++ on a realistic network, `COthin`: 37 species
    and 238 reactions of carbon, oxygen, hydrogen, helium and silicon chemistry,
    including reactions given over several temperature ranges.

    The `kokkos_ode` collection writes a C++20 solver that integrates the rate
    equations with the VODE (variable-order BDF) method from the
    [header-only-integrators](https://github.com/BenWibking/header-only-integrators)
    library. The build steps need `cmake`, `clang++` and network access, because
    CMake downloads the integrators library on first configure.

    As before, every command runs through `run()`, which fails the cell on a
    nonzero exit code.
    """)
    return


@app.cell
def _(COTHIN, load_network, mo):
    net = load_network(COTHIN)
    mo.ui.table(
        [
            {"attribute": "label", "value": net.label},
            {"attribute": "species", "value": str(net.species.count)},
            {"attribute": "reactions", "value": str(net.reactions.count)},
            {"attribute": "elements", "value": ", ".join(net.elements.symbols())},
        ],
        selection=None,
        show_column_summaries=False,
    )
    return (net,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Generating the code
    """)
    return


@app.cell
def _(COTHIN, JAFFGEN, run, scratch):
    run(
        [
            JAFFGEN,
            "--network",
            COTHIN,
            "--template",
            "kokkos_ode",
            "--outdir",
            scratch / "kokkos_ode",
        ]
    )
    # jaffgen writes into a new <network>_<template>_<date>-<time> folder.
    (out,) = (scratch / "kokkos_ode").iterdir()
    return (out,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Five files come out:

    - `chemistry_ode.hpp` defines a `ChemistryODE` struct: the species names,
      the network's parameters (`ChemistryODE::params`), `rhs()` (rate
      coefficients, then $\mathrm{d}n_i/\mathrm{d}t$ and $\mathrm{d}T/\mathrm{d}t$)
      and `jacobian()` (the exact derivatives, which stiff solvers need). It is
      large because every rate and Jacobian entry is written out explicitly, with
      common subexpressions factored into temporaries;
    - `user_functions.hpp` holds stubs for functions the network calls but does
      not define (COthin has none, so it is empty);
    - `chemistry_ode.cpp` is the driver: a block of inputs to edit, then the
      integration, a `solution.dat` file and a conservation check;
    - `CMakeLists.txt` fetches the integrators library and builds the driver;
    - `config.toml` records the settings `jaffgen` was run with.
    """)
    return


@app.cell
def _(file_table, out):
    file_table(out)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    The state is the species densities followed by the gas temperature. The
    network's parameters (here the visual extinction `av`, the cosmic-ray
    ionisation rate `crate` and the dust temperature `tdust`) are members of
    `ChemistryParameters`, NaN until set:
    """)
    return


@app.cell
def _(mo, out):
    _hpp = (out / "chemistry_ode.hpp").read_text()
    _start = _hpp.index("struct ChemistryParameters")
    _stop = _hpp.index("// Index of a species")
    mo.md("```cpp\n" + _hpp[_start:_stop].rstrip() + "\n```")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Setting the inputs

    The driver starts with a block of inputs: one line per species (initial
    density, 0 by default), one per parameter (a placeholder 1.0), the initial
    temperature, the end time, whether to evolve the temperature, and the
    tolerances. We fill it in with test values that exercise all reactions rather
    than a physical cloud: every species at $10^{-10}$ and H at $10^{-6}$
    cm$^{-3}$, typical parameters, the temperature fixed at 1000 K (some COthin
    rate fits overflow well below 200 K), and $10^{12}$ s (about 32 000 years) of
    evolution:
    """)
    return


@app.cell
def _(mo, net, np, out):
    initial = np.full(net.species.count, 1e-10)
    initial[0] = 1e-6
    _parameters = {"av": 1.0, "crate": 1.3e-17, "tdust": 20.0}

    _driver = out / "chemistry_ode.cpp"
    _text = _driver.read_text()
    for _i, _n in enumerate(initial):
        _text = _text.replace(f"  y[{_i}] = 0.0;", f"  y[{_i}] = {_n:.1e};", 1)
    for _name, _value in _parameters.items():
        _text = _text.replace(
            f"params.{_name} = 1.0;", f"params.{_name} = {_value:.6g};", 1
        )
    _text = _text.replace("temperature = 100.0;", "temperature = 1000.0;", 1)
    _text = _text.replace("t_end = 3.15576e13;", "t_end = 1e12;", 1)
    _driver.write_text(_text)

    _lines = _text.splitlines()
    _start = next(i for i, s in enumerate(_lines) if "inputs: edit here" in s)
    _stop = next(i for i, s in enumerate(_lines) if "const char *output" in s)
    mo.md("```cpp\n" + "\n".join(_lines[_start : _stop + 1]) + "\n```")
    return (initial,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Building

    Configure with CMake. The integrators library uses a C++ construct that GCC
    rejects, so the generated `CMakeLists.txt` picks `clang++` unless you name a
    compiler yourself (`-DCMAKE_CXX_COMPILER=...` or the `CXX` environment
    variable):
    """)
    return


@app.cell
def _(out, run):
    run(["cmake", "-S", ".", "-B", "build"], cwd=out)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    Then compile. Some rates are evaluated at the edge of their temperature range
    and underflow to zero, so clang may print `-Wliteral-range` warnings; they are
    harmless.
    """)
    return


@app.cell
def _(out, run):
    run(["cmake", "--build", "build", "-j"], cwd=out)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Running

    The driver takes no arguments: it stops if a parameter is unset, integrates
    with the analytic Jacobian and prints the final state:
    """)
    return


@app.cell
def _(out, run):
    run(["./build/chemistry_ode"], cwd=out)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    The end of the output is a conservation check. Chemistry moves atoms
    between species but never creates or destroys them, so the total amount of
    each element, and the net charge, must be the same before and after. Here they
    agree to round-off: the generated equations conserve every element exactly.

    `solution.dat` holds the state at $t = 0$ and at 100 log-spaced times up to
    the end, one row each: time, temperature, then every species. Its last row,
    labelled with the network's species names:
    """)
    return


@app.cell
def _(initial, mo, net, np, out):
    solution = np.loadtxt(out / "solution.dat")
    final = solution[-1, 2:]
    mo.ui.table(
        [
            {
                "index": i,
                "species": sp.name,
                "initial": f"{n0:.3e}",
                "final": f"{n1:.3e}",
            }
            for i, (sp, n0, n1) in enumerate(zip(net.species, initial, final))
        ],
        selection=None,
        page_size=40,
    )
    return (solution,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    The evolution of the most abundant species. Values far below $10^{-20}$ (the
    solver's absolute tolerance) are effectively zero, and tiny negative values
    are solver noise at that level, so both are clipped to $10^{-20}$ for
    plotting.
    """)
    return


@app.cell
def _(PLOT_SIZE, net, np, plt, solution):
    _floor = 1e-20
    _t, _n = solution[1:, 0], solution[1:, 2:]
    _fig, _ax = plt.subplots(figsize=(PLOT_SIZE[0] * 1.4, PLOT_SIZE[1] * 1.2))
    for _i in np.argsort(_n[-1])[::-1][:8]:
        _ax.loglog(_t, np.maximum(_n[:, _i], _floor), label=net.species.names()[_i])
    _ax.set_xlabel("time [s]")
    _ax.set_ylabel("number density [cm$^{-3}$]")
    _ax.legend(fontsize="small", ncols=2)
    _fig.tight_layout()
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    **Next:** `06_rate_tables.py` tabulates rate coefficients on a temperature
    grid instead of generating code for them.
    """)
    return


@app.cell(hide_code=True)
def _():
    import contextlib
    import io
    import logging
    import subprocess
    import sys
    import tempfile
    from pathlib import Path
    from typing import Dict, List

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np

    from jaff import Network

    # JAFF resets its logger to INFO on every JaffLogger(); disabling WARNING and
    # below is the only switch that sticks.
    logging.disable(logging.WARNING)

    COTHIN = mo.notebook_dir().parent / "networks" / "COthin" / "react_COthin.jet"
    # The jaffgen console script installed next to this Python interpreter.
    JAFFGEN = Path(sys.executable).parent / "jaffgen"
    PLOT_SIZE = (5.5, 3.5)

    # Returned from this cell so it stays alive: a TemporaryDirectory deletes its
    # folder when garbage-collected, and _-prefixed names are dropped after the cell.
    scratch_dir = tempfile.TemporaryDirectory()
    scratch = Path(scratch_dir.name)

    def load_network(path: Path, **kwargs) -> Network:
        """Load a network, hiding the banner JAFF prints on stdout."""
        with contextlib.redirect_stdout(io.StringIO()):
            return Network(str(path), **kwargs)

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

    return (
        COTHIN,
        JAFFGEN,
        PLOT_SIZE,
        file_table,
        load_network,
        mo,
        np,
        plt,
        run,
        scratch,
        scratch_dir,
    )


if __name__ == "__main__":
    app.run()
