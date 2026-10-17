# ABOUTME: Example notebook: tabulate rate coefficients on an adaptive temperature grid
# ABOUTME: with Network.write_table, inspect the HDF5 file and check accuracy (marimo)

import marimo

__generated_with = "0.24.2"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # 06 · Rate tables

    Generated solver code evaluates every rate coefficient from its formula at
    every step: powers, exponentials, piecewise branches. For large networks on
    GPUs that is expensive, and some simulation codes (e.g. quokka) instead
    **look rates up** in a precomputed table and interpolate.

    A JAFF rate table stores $k_i(T)$ for every reaction on one temperature grid,
    uniformly spaced in $\log T$. At run time, $\log k$ is interpolated linearly
    in $\log T$ between neighbouring grid points. This notebook builds such a
    table for `uclchem_small_gas` (562 reactions), looks inside the file, and
    checks how accurate the interpolation is.
    """)
    return


@app.cell
def _(UCLCHEM, load_network, mo):
    net = load_network(UCLCHEM)
    mo.ui.table(
        [
            {"attribute": "label", "value": net.label},
            {"attribute": "species", "value": str(net.species.count)},
            {"attribute": "reactions", "value": str(net.reactions.count)},
        ],
        selection=None,
        show_column_summaries=False,
    )
    return (net,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Which reactions can be tabulated

    A one-dimensional table can only hold rates that depend on temperature alone.
    Before sorting reactions, JAFF sets the visual extinction `av` to 0 and the
    cosmic-ray ionisation rate `crate` to 1 (so cosmic-ray rates are tabulated
    per unit `crate` and scaled at run time). Then:

    - rates that depend only on `tgas` are **tabulated**;
    - rates that still depend on something else (grain properties, densities,
      radiation) cannot be, and are left out;
    - rates that are constant over the temperature range are left out too,
      since a table adds nothing; pass `include_all=True` to keep them.
    """)
    return


@app.cell
def _(mo, net, sympy, tabulable):
    _tgas = sympy.Symbol("tgas")
    _counts = {"depends only on temperature": 0, "constant": 0, "depends on more": 0}
    for _rxn in net.reactions:
        _symbols = tabulable(_rxn.rate).free_symbols
        if not _symbols:
            _counts["constant"] += 1
        elif _symbols == {_tgas}:
            _counts["depends only on temperature"] += 1
        else:
            _counts["depends on more"] += 1
    mo.ui.table(
        [
            {"rate (with av = 0, crate = 1)": k, "reactions": v}
            for k, v in _counts.items()
        ],
        selection=None,
        show_column_summaries=False,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Building the table

    `Network.write_table` builds the grid adaptively. It starts from `nT`
    log-spaced temperatures and repeatedly inserts a point halfway (in
    $\log T$) between every pair. At each new point it compares the exact rate
    with the value interpolated from its neighbours, and stops once the largest
    relative error over all reactions is below `err_tol`. `nT_max` caps the
    grid size so a network that cannot reach `err_tol` stops with a warning
    instead of refining forever.

    We ask for 10 K to 10⁴ K with 5 % accuracy. The file extension picks the
    format: `.hdf5` here, `.txt` for a plain-text table.
    """)
    return


@app.cell
def _(net, scratch, time):
    table_file = scratch / "uclchem_rates.hdf5"
    _start = time.perf_counter()
    net.write_table(table_file, T_min=10.0, T_max=1e4, err_tol=0.05)
    build_seconds = time.perf_counter() - _start
    return build_seconds, table_file


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Inside the file

    Everything lives in one HDF5 group, `reaction_coeff`, in the layout quokka
    reads:

    - attributes describe the grid: the input variable (`temperature`, in K),
      the spacing (`log`) and the end points `xlo`, `xhi`; the grid itself is not
      stored, since it is fully determined by these and the number of points;
    - `data` holds one row per tabulated reaction, one column per temperature;
    - `output_names` labels each row with the reaction type, reactants and
      products; `output_units` gives each row's units.
    """)
    return


@app.cell
def _(build_seconds, h5py, mo, table_file):
    with h5py.File(table_file) as _f:
        _group = _f["reaction_coeff"]
        _attrs = {k: _group.attrs[k] for k in _group.attrs}
        _shape = _group["data"].shape
    mo.ui.table(
        [{"item": f"attribute {k}", "value": str(v)} for k, v in _attrs.items()]
        + [
            {"item": "data shape (reactions, temperatures)", "value": str(_shape)},
            {"item": "file size", "value": f"{table_file.stat().st_size / 1e6:.1f} MB"},
            {"item": "build time", "value": f"{build_seconds:.1f} s"},
        ],
        selection=None,
        show_column_summaries=False,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    The first few rows and their labels:
    """)
    return


@app.cell
def _(h5py, mo, np, table_file):
    with h5py.File(table_file) as _f:
        table_names = [n.decode() for n in _f["reaction_coeff/output_names"][:]]
        table_units = [u.decode() for u in _f["reaction_coeff/output_units"][:]]
        table_data = _f["reaction_coeff/data"][:]
        table_lo = float(_f["reaction_coeff"].attrs["xlo"][0])
        table_hi = float(_f["reaction_coeff"].attrs["xhi"][0])
    table_temps = np.logspace(np.log10(table_lo), np.log10(table_hi), table_data.shape[1])

    mo.ui.table(
        [
            {"row": i, "output_name": name, "unit": unit, "k(10 K)": f"{row[0]:.3e}"}
            for i, (name, unit, row) in enumerate(
                zip(table_names, table_units, table_data)
            )
        ],
        selection=None,
        page_size=8,
    )
    return table_data, table_names, table_temps


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## How accurate is it?

    The table is only exact at its grid points. To check it between them, the
    next cell takes one reaction, interpolates $\log k$ at the midpoint of
    every grid interval and compares it with the exact rate there. Pick a
    reaction; the default is one with a piecewise rate.
    """)
    return


@app.cell
def _(mo, net, table_names, table_row_name):
    _by_name = {table_row_name(rxn): rxn for rxn in net.reactions}
    _options = {name: name for name in table_names if name in _by_name}
    reaction_picker = mo.ui.dropdown(
        options=_options,
        value="2_body rate coefficient: {'H2': 1, 'C': 1} --> {'CH': 1, 'H': 1}",
        label="Reaction",
    )
    reaction_by_name = _by_name
    reaction_picker
    return reaction_by_name, reaction_picker


@app.cell
def _(
    PLOT_SIZE,
    lambdify,
    np,
    plt,
    reaction_by_name,
    reaction_picker,
    sympy,
    table_data,
    table_names,
    table_temps,
    tabulable,
):
    _rxn = reaction_by_name[reaction_picker.value]
    _row = table_data[table_names.index(reaction_picker.value)]
    _exact = lambdify(sympy.Symbol("tgas"), tabulable(_rxn.rate), "numpy")

    _mid = np.sqrt(table_temps[1:] * table_temps[:-1])
    _interp = np.exp(0.5 * (np.log(_row[1:]) + np.log(_row[:-1])))
    _truth = np.array([float(_exact(t)) for t in _mid])
    _rel_err = np.abs(_interp - _truth) / np.abs(_truth)

    _fig, (_top, _bottom) = plt.subplots(
        2, 1, figsize=(PLOT_SIZE[0], PLOT_SIZE[1] * 1.6), sharex=True
    )
    _top.loglog(table_temps, _row, label="table")
    _top.set_ylabel("k")
    _top.set_title(_rxn.verbatim)
    _bottom.loglog(_mid, _rel_err)
    _bottom.axhline(0.05, color="k", ls="--", lw=0.8, label="err_tol")
    _bottom.set_xlabel("T (K)")
    _bottom.set_ylabel("relative error")
    _bottom.legend()
    _fig.tight_layout()
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    For smooth rates the error between grid points is tiny. It peaks where the
    rate has a kink: a boundary between two temperature ranges, or a `clip`
    cutoff. For `H2 + C -> CH + H` it is the cutoff at 300 K: the rate is only
    given above 300 K, so below it is held at its 300 K value, and the
    $e^{-11700/T}$ curve meets that constant at an angle. A kink cannot be
    represented by linear interpolation, and
    halving the spacing only halves the error there, so these few points decide
    how many grid points the whole table needs.

    ## Accuracy versus size

    That first-order behaviour shows up directly when we rebuild the table for a
    few tolerances. Each halving of `err_tol` roughly doubles the grid, the file
    size and the build time.
    """)
    return


@app.cell
def _(h5py, mo, net, scratch, time):
    _rows = []
    for _tol in (0.2, 0.1, 0.05):
        _path = scratch / f"tol_{_tol}.hdf5"
        _start = time.perf_counter()
        net.write_table(_path, T_min=10.0, T_max=1e4, err_tol=_tol)
        _seconds = time.perf_counter() - _start
        with h5py.File(_path) as _f:
            _n = _f["reaction_coeff/data"].shape[1]
        _rows.append(
            {
                "err_tol": _tol,
                "grid points": _n,
                "file size (MB)": round(_path.stat().st_size / 1e6, 1),
                "build time (s)": round(_seconds, 1),
            }
        )
    mo.ui.table(_rows, selection=None, show_column_summaries=False)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    The default `err_tol=0.01` needs about 16 000 points for this network. For
    the full KIDA network (about 7 600 reactions), several rates jump by more than
    an order of magnitude across the 1 K gaps between their temperature ranges,
    and reaching 10 % takes roughly 64 000 points; `nT_max` (default 16 385)
    stops the refinement before that. Improving this is tracked in
    [issue #217](https://github.com/jaff-chemistry/jaff/issues/217).

    One more caveat: `output_units` currently says `cm^3 s^-1` for every row,
    which is only right for two-body reactions. Cosmic-ray rates are per second
    (per unit `crate`) and three-body rates are in cm⁶ s⁻¹.

    This is the last notebook in the series.
    """)
    return


@app.cell(hide_code=True)
def _():
    import contextlib
    import io
    import logging
    import tempfile
    import time
    from collections import Counter
    from pathlib import Path

    import h5py
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import sympy
    from sympy import lambdify

    from jaff import Network

    # JAFF resets its logger to INFO on every JaffLogger(); disabling WARNING and
    # below is the only switch that sticks.
    logging.disable(logging.WARNING)

    UCLCHEM = (
        mo.notebook_dir().parent
        / "networks"
        / "uclchem_small_gas"
        / "uclchem_small_gas_network.jet"
    )
    PLOT_SIZE = (5.5, 3.5)

    # Returned from this cell so it stays alive: a TemporaryDirectory deletes its
    # folder when garbage-collected, and _-prefixed names are dropped after the cell.
    scratch_dir = tempfile.TemporaryDirectory()
    scratch = Path(scratch_dir.name)

    def load_network(path: Path, **kwargs) -> Network:
        """Load a network, hiding the banner JAFF prints on stdout."""
        with contextlib.redirect_stdout(io.StringIO()):
            return Network(str(path), **kwargs)

    def tabulable(rate: sympy.Expr) -> sympy.Expr:
        """The rate as write_table sees it: av set to 0 and crate to 1."""
        return rate.subs({"av": 0.0, "crate": 1.0})

    def table_row_name(rxn) -> str:
        """The output_names label write_table gives *rxn*."""
        rtype = "2_body" if rxn.type == "unknown" else rxn.type
        reactants = dict(Counter(rxn.reactants.names()))
        products = dict(Counter(rxn.products.names()))

        return f"{rtype} rate coefficient: {reactants} --> {products}"

    return (
        PLOT_SIZE,
        UCLCHEM,
        h5py,
        lambdify,
        load_network,
        mo,
        np,
        plt,
        scratch,
        scratch_dir,
        sympy,
        table_row_name,
        tabulable,
        time,
    )


if __name__ == "__main__":
    app.run()
