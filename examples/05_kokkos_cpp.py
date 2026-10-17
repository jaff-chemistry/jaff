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

    The `kokkos_ode` collection writes a C++17 solver that integrates the rate
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
    out = scratch / "kokkos_ode"
    run([JAFFGEN, "--network", COTHIN, "--template", "kokkos_ode", "--outdir", out])
    return (out,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Three files come out:

    - `chemistry_ode.hpp` defines a `ChemistryODE` struct: species indices,
      `rhs()` (rate coefficients, then $\mathrm{d}n_i/\mathrm{d}t$) and
      `jacobian()` (the exact derivatives $\partial \dot n_i/\partial n_j$,
      which stiff solvers need). It is large because every rate and Jacobian
      entry is written out explicitly, with common subexpressions factored into
      temporaries;
    - `chemistry_ode.cpp` is the driver: it sets initial conditions, calls the
      integrator, prints the result and checks conservation;
    - `CMakeLists.txt` fetches the integrators library and builds the driver.
    """)
    return


@app.cell
def _(file_table, out):
    file_table(out)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    The start of the header shows the generated species indices (C++ identifiers
    cannot contain `+` or `-`, so ions get a `j` suffix and anions a `k`):
    """)
    return


@app.cell
def _(mo, out):
    _head = (out / "chemistry_ode.hpp").read_text().splitlines()[:40]
    mo.md("```cpp\n" + "\n".join(_head) + "\n```")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Building

    Configure with CMake, choosing `clang++` as the compiler. The integrators
    library currently uses a C++ construct that recent GCC versions reject, so
    clang is the safe choice.
    """)
    return


@app.cell
def _(out, run):
    run(["cmake", "-S", ".", "-B", "build", "-DCMAKE_CXX_COMPILER=clang++"], cwd=out)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    Then compile. Some rates are evaluated at the edge of their temperature range
    and underflow to zero, so clang prints `-Wliteral-range` warnings; they are
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

    The driver takes up to three positional arguments, `t_end [ic_file]
    [out_file]`, plus `--analytic-jac` to use the generated Jacobian instead of
    a finite-difference one. Without an initial-conditions file it starts every
    species at $10^{-10}$ and the first one (H) at $10^{-6}$, test values that
    exercise all reactions rather than a physical cloud.

    We write those same initial conditions to a file, so the run is explicit,
    integrate for $10^{12}$ s (about 32 000 years) with the analytic Jacobian,
    and ask for the final densities in an output file:
    """)
    return


@app.cell
def _(net, np, out, run):
    initial = np.full(net.species.count, 1e-10)
    initial[0] = 1e-6
    np.savetxt(out / "initial.txt", initial)

    run(
        [
            "./build/chemistry_ode",
            "1e12",
            "initial.txt",
            "final.txt",
            "--analytic-jac",
        ],
        cwd=out,
    )
    return (initial,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    The end of the output is a conservation check. Chemistry moves atoms
    between species but never creates or destroys them, so the total amount of
    each element, and the net charge, must be the same before and after. Here they
    agree to about $10^{-14}$, which is round-off: the generated equations
    conserve every element exactly.

    The final densities, read back from `final.txt` and labelled with the
    network's species names:
    """)
    return


@app.cell
def _(initial, mo, net, np, out):
    final = np.loadtxt(out / "final.txt")
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
    return (final,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    The bar chart compares the two. Species whose final value is far below
    $10^{-20}$ (the solver's absolute tolerance) are effectively zero, and tiny
    negative values are solver noise at that level, so both are clipped to
    $10^{-20}$ for plotting.
    """)
    return


@app.cell
def _(PLOT_SIZE, final, initial, net, np, plt):
    _floor = 1e-20
    _x = np.arange(net.species.count)
    _fig, _ax = plt.subplots(figsize=(PLOT_SIZE[0] * 1.8, PLOT_SIZE[1]))
    _ax.bar(_x - 0.2, np.maximum(initial, _floor), width=0.4, label="initial")
    _ax.bar(_x + 0.2, np.maximum(final, _floor), width=0.4, label="final")
    _ax.set_yscale("log")
    _ax.set_xticks(_x, net.species.names(), rotation=90)
    _ax.set_ylabel("number density")
    _ax.legend()
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
