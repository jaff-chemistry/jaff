# ABOUTME: Example notebook: load a chemical network and inspect its species, elements
# ABOUTME: and stoichiometry; first stop for new JAFF users (marimo, runs headless too)

import marimo

__generated_with = "0.24.2"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    # 01 · Network basics

    A chemical network is a list of reactions, each with a rate coefficient. JAFF
    reads such a list, figures out which species take part, checks that the
    chemistry is consistent, and turns everything into symbolic expressions that
    later notebooks use for plotting, tabulating and generating solver code.

    This notebook covers the first step: **loading a network and looking at what
    JAFF built from it.** We use `networks/demos/demo1.jet`, a small demo network
    written to show off the input formats JAFF understands.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## The input file

    Below is the raw `demo1.jet`. It is deliberately messy: one file mixes several
    community formats. JAFF detects the format **line by line**, so all of these
    can live together:

    - **PRIZMO** lines (`H+ + E -> H  [0, 300]  rate`): reactants, products, an
      optional temperature range in brackets, then the rate expression;
    - a **KIDA** fixed-width line (`N2 CR N N ...`);
    - a **UMIST/UDFA** colon-separated line (`326:CE:CO:N2+:...`);
    - **KROME** comma-separated lines, both the default layout and a custom one
      declared with `@format:`.

    The `VARIABLES{...}` block defines shorthands (`invt`, `t32`) that are
    substituted into the rate expressions while parsing. Rates may use
    *primitive variables* such as `Tgas` (gas temperature), `crate` (cosmic-ray
    ionisation rate) and `av` (visual extinction); these stay symbolic.
    """)
    return


@app.cell
def _(DEMO1, mo):
    mo.md("```text\n" + DEMO1.read_text() + "\n```")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Loading the network

    Everything starts with one constructor call, `Network(path)`. It runs a fixed
    pipeline:

    1. **Read** the file and detect each line's format;
    2. **Parse & assemble**: every new species name becomes a `Specie`, every line
       becomes a `Reaction`;
    3. **Standardise symbols**: densities in rate expressions are rewritten onto
       one vector, `nden[i, 0]`, where `i` is the species index;
    4. **Validate**: mass and charge conservation, sinks/sources, duplicates;
    5. **Finalise**: build stoichiometry matrices and the set of elements.

    JAFF prints a banner when it loads; the helper cell at the bottom of this
    notebook hides it so the output stays readable.
    """)
    return


@app.cell
def _(DEMO1, load_network):
    net = load_network(DEMO1)
    return (net,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    The result is a `Network` object. Its basic metadata tells us what was built.
    The label defaults to the file name without its extension.
    """)
    return


@app.cell
def _(mo, net):
    mo.ui.table(
        [
            {"attribute": "net.label", "value": net.label},
            {"attribute": "net.filename", "value": net.filename.name},
            {"attribute": "net.species.count", "value": net.species.count},
            {"attribute": "net.reactions.count", "value": net.reactions.count},
            {
                "attribute": "net.elements.symbols()",
                "value": ", ".join(net.elements.symbols()),
            },
        ],
        selection=None,
        show_column_summaries=False,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Species

    `net.species` is an ordered catalogue of `Specie` objects. The order matters:
    a species' `index` is its position in the density vector `nden`, and therefore
    its slot in every array the generated solver code uses.

    JAFF works out each species' mass and charge from its name alone. Note that
    the input file wrote electrons as `E`; JAFF normalises that to `e-`. Masses
    are in grams, charges in units of the elementary charge.
    """)
    return


@app.cell
def _(mo, net):
    species_rows = [
        {
            "index": sp.index,
            "name": sp.name,
            "mass [g]": f"{sp.mass:.4e}",
            "charge": sp.charge,
            "elements": ", ".join(sp.elements.symbols()),
        }
        for sp in net.species
    ]
    mo.ui.table(species_rows, selection=None, page_size=20)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    Each species can also render itself as LaTeX via `Specie.latex()`, which is
    handy for plot labels and reports.
    """)
    return


@app.cell
def _(mo, net):
    mo.ui.table(
        [{"name": sp.name, "LaTeX": mo.md(f"${sp.latex()}$")} for sp in net.species],
        selection=None,
        page_size=20,
        show_column_summaries=False,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Querying the catalogues

    The catalogues offer helpers for common questions, so you rarely need to
    loop by hand. For example, `names()` lists every species,
    `charged("name")` keeps only the ions and electrons, and
    `elements.symbols()` lists the chemical elements present.
    """)
    return


@app.cell
def _(mo, net):
    mo.ui.table(
        [
            {"call": "net.species.names()", "result": str(net.species.names())},
            {
                "call": 'net.species.charged("name")',
                "result": str(net.species.charged("name")),
            },
            {"call": "net.elements.symbols()", "result": str(net.elements.symbols())},
        ],
        selection=None,
        show_column_summaries=False,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Stoichiometry

    JAFF stores the network's structure in two integer matrices of shape
    `(n_reactions, n_species)`:

    - `reactant_matrix[i, j]`: how many of species *j* reaction *i* consumes;
    - `product_matrix[i, j]`: how many it produces.

    Their difference, `product_matrix - reactant_matrix`, is the net change in
    each species per reaction. It is the backbone of the rate equations
    $\mathrm{d}n_j/\mathrm{d}t = \sum_i (P_{ij} - R_{ij})\, k_i \prod n_\text{reactants}$
    that the code generators emit. The table below shows it with one row per
    reaction (zeros are left blank so the pattern stands out).
    """)
    return


@app.cell
def _(mo, net):
    net_change = net.product_matrix - net.reactant_matrix
    names = net.species.names()
    stoich_rows = [
        {"reaction": rxn.verbatim}
        | {name: (int(v) if v else "") for name, v in zip(names, row)}
        for rxn, row in zip(net.reactions, net_change)
    ]
    mo.ui.table(stoich_rows, selection=None, page_size=20)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    Every reaction should conserve mass and charge, so each row's changes,
    weighted by species mass (or charge), must sum to zero. JAFF already checked
    this during loading (phase 4); here is the charge check done by hand for every
    reaction:
    """)
    return


@app.cell
def _(net):
    import numpy as np

    charges = np.array([sp.charge for sp in net.species])
    charge_change = (net.product_matrix - net.reactant_matrix) @ charges
    assert not charge_change.any(), "a reaction does not conserve charge"
    charge_change
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    All zeros: every reaction conserves charge.

    **Next:** `02_reactions_and_rates.py` looks inside the reactions themselves:
    their rate expressions, temperature ranges and cross-sections.
    """)
    return


@app.cell(hide_code=True)
def _():
    import contextlib
    import io
    import logging
    from pathlib import Path

    import marimo as mo

    from jaff import Network

    # JAFF resets its logger to INFO on every JaffLogger(); disabling WARNING and
    # below is the only switch that sticks.
    logging.disable(logging.WARNING)

    DEMO1 = mo.notebook_dir().parent / "networks" / "demos" / "demo1.jet"

    def load_network(path: Path, **kwargs) -> Network:
        """Load a network, hiding the banner JAFF prints on stdout."""
        with contextlib.redirect_stdout(io.StringIO()):
            return Network(str(path), **kwargs)

    return DEMO1, load_network, mo


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
