# ABOUTME: Example notebook: inspect reactions, symbolic rate coefficients, code output,
# ABOUTME: multi-range (piecewise) rates and photo cross-sections (marimo, runs headless)

import marimo

__generated_with = "0.24.2"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    # 02 · Reactions and rates

    `01_network_basics.py` looked at the species a network contains. This notebook
    looks at the reactions that connect them: how JAFF identifies and classifies
    each reaction, why every rate coefficient is kept as a **symbolic expression**
    rather than a number, how that expression is printed in other programming
    languages, and how JAFF stitches a reaction given over several temperature
    ranges into one rate. It ends with photo-reactions, whose rates come from
    cross-section data.

    We reuse `networks/demos/demo1.jet`.
    """)
    return


@app.cell
def _(DEMO1, load_network):
    net = load_network(DEMO1)
    return (net,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Every reaction at a glance

    `net.reactions` is an ordered catalogue of `Reaction` objects. The table lists
    the key attributes of each one:

    - `verbatim`: the human-readable equation;
    - `type`: the category the parser assigned while reading the line;
    - `tmin` / `tmax`: the temperature range (K) over which the rate is valid,
      `None` meaning unbounded;
    - `reactants` / `products`: species catalogues, so they repeat a species when
      it appears twice (`N + N -> N2`).
    """)
    return


@app.cell
def _(mo, net):
    mo.ui.table(
        [
            {
                "verbatim": rxn.verbatim,
                "type": rxn.type,
                "tmin": rxn.tmin,
                "tmax": rxn.tmax,
                "reactants": ", ".join(rxn.reactants.names()),
                "products": ", ".join(rxn.products.names()),
            }
            for rxn in net.reactions
        ],
        selection=None,
        page_size=20,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ### Pseudo-species and reaction types

    Some reactants start with an underscore: `_PHOTON` and `_CR`. These are
    **pseudo-species**. They record what drives the reaction (a photon, a cosmic
    ray) and make it a distinct reaction from, say, collisional ionisation, but
    they are not real species: they never appear in `net.species` and get no
    density.

    The `type` comes from the parser, not from evaluating anything. A `_PHOTON`
    reactant makes a reaction `photo`, a `_CR` reactant makes it `cosmic_ray`, and
    three or more real reactants make it `3_body`. For PRIZMO lines without a
    pseudo-species the parser then looks at the rate: one using `crate` is
    `cosmic_ray`, one using `av` (visual extinction, i.e. shielded starlight) is
    `photo`. That is why `H -> H+ + e-` (rate `0.46*crate`) is a cosmic-ray
    reaction and `C -> C+ + e-` (rate depends on `av`) is a photo-reaction.
    Everything else is `unknown`, which simply means "ordinary two-body gas-phase
    chemistry".
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## One reaction up close

    Indexing the catalogue gives a single `Reaction`. We take the first one,
    radiative recombination of hydrogen, `H+ + e- -> H`.

    A reaction has two identities:

    - `serialized` sorts the species **names** on each side and joins them
      (`.` between species, `__` between the two sides). JAFF uses it for `==`
      and hashing, so two reactions are equal when they have the same reactants
      and products, whatever order the file wrote them in.
    - `serialized_exploded` does the same with each species' atom-level form, so
      isomers (e.g. `HCO+` and `HOC+`) collapse together. It is used to spot
      isomer twins.

    It also has two indices: `index` is the row number in the source file (it can
    have gaps), while `catalogue_index` is the dense position in `net.reactions`
    that generated code uses as an array subscript.
    """)
    return


@app.cell
def _(mo, net):
    rxn = net.reactions[0]
    mo.ui.table(
        [
            {"attribute": "verbatim", "value": rxn.verbatim},
            {"attribute": "serialized", "value": rxn.serialized},
            {"attribute": "serialized_exploded", "value": rxn.serialized_exploded},
            {"attribute": "index", "value": str(rxn.index)},
            {"attribute": "catalogue_index", "value": str(rxn.catalogue_index)},
            {"attribute": "tmin, tmax", "value": f"{rxn.tmin}, {rxn.tmax}"},
            {
                "attribute": "rxn == net.reactions[0]",
                "value": str(rxn == net.reactions[0]),
            },
            {
                "attribute": "rxn == net.reactions[1]",
                "value": str(rxn == net.reactions[1]),
            },
        ],
        selection=None,
        show_column_summaries=False,
    )
    return (rxn,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## The rate is symbolic

    `rxn.rate` is a [SymPy](https://www.sympy.org) expression in the gas
    temperature `tgas`, not a number. Keeping it symbolic means one rate can be
    differentiated exactly, printed in any language, or evaluated on a
    temperature grid later, without re-parsing.

    The file line was `H+ + E -> H  [0, 300]  3.61e-12*t32**(-0.75)`. Three things
    happened on the way in:

    1. the `VARIABLES{}` shorthand `t32 = Tgas / 3d2` was substituted and the
       constants folded, so the rate is $2.6\times10^{-10}\,T^{-0.75}$;
    2. the lower bound `0` means "no lower bound", so `tmin` is `None`;
    3. above the upper bound (300 K) the reaction's temperature cutoff applies.
       The default, `clip`, holds the rate at its value at 300 K, so the rate
       becomes a `Piecewise`: the fit below 300 K, a constant above it.
    """)
    return


@app.cell
def _(rxn):
    rxn.rate
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Because it is a SymPy object, the usual SymPy tools work on it directly. For
    example, the exact temperature derivative $\mathrm{d}k/\mathrm{d}T$, which
    solvers need for the Jacobian when temperature is evolved:
    """)
    return


@app.cell
def _(rxn):
    import sympy

    sympy.diff(rxn.rate, sympy.Symbol("tgas"))
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Printing a rate in other languages

    `Reaction.get_code(lang=...)` prints the rate as source code for a target
    language. Notice how the `Piecewise` becomes a ternary in C/C++, `merge` in
    Fortran and an `if`/`else` expression in Python and Rust. This is the same
    printer the code generators in `04_codegen_jaffgen.py` use.

    Pick a reaction and a language:
    """)
    return


@app.cell
def _(mo, net):
    reaction_picker = mo.ui.dropdown(
        options={rxn.verbatim: i for i, rxn in enumerate(net.reactions)},
        value=net.reactions[0].verbatim,
        label="Reaction",
    )
    language_picker = mo.ui.dropdown(
        options=["c", "cpp", "fortran", "python", "rust", "julia", "r"],
        value="fortran",
        label="Language",
    )
    mo.hstack([reaction_picker, language_picker], justify="start")
    return language_picker, reaction_picker


@app.cell
def _(language_picker, mo, net, reaction_picker):
    _picked = net.reactions[reaction_picker.value]
    mo.md(
        f"```{language_picker.value}\n{_picked.get_code(lang=language_picker.value)}\n```"
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Flux expressions

    A rate coefficient $k$ alone is not how fast a reaction proceeds; the flux is
    $k$ times the densities of the reactants. `get_flux_expression()` writes that
    product the way generated code indexes its arrays: `k[i]` for the rate
    coefficient of reaction `i` and `y[idx_X]` for the density of species `X`
    (ions get a `j` suffix, so `H+` is `idx_hj`). The names are configurable
    through `rate_variable`, `species_variable` and `brackets`.
    """)
    return


@app.cell
def _(mo, net):
    mo.ui.table(
        [
            {
                "reaction": rxn.verbatim,
                "flux": rxn.get_flux_expression(idx=rxn.catalogue_index),
            }
            for rxn in net.reactions
        ],
        selection=None,
        page_size=20,
        show_column_summaries=False,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Plotting rates

    `jaff.plotting.plot_rates` evaluates rate expressions on a temperature grid
    and plots them. It accepts reactions, SymPy expressions, or lists of either.

    The demo network has four temperature-dependent rates. Three are electron
    recombinations, which are fast at low temperature and slow down as the gas
    heats. Note the kink in the H⁺ recombination rate at 300 K, where `clip`
    takes over.
    """)
    return


@app.cell
def _(PLOT_STYLE, net, plot_rates):
    recombinations = [
        net.reactions.get("H+ + e- -> H"),
        net.reactions.get("C+ + e- -> C"),
        net.reactions.get("N2+ + e- -> N + N"),
    ]
    _fig, _ = plot_rates(
        recombinations,
        tmin=10,
        tmax=1e4,
        labels=[rxn.verbatim for rxn in recombinations],
        style=PLOT_STYLE,
        show=False,
        close=False,
    )
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    The fourth, collisional dissociation of H₂ by electrons, behaves completely
    differently. Its rate $5.6\times10^{-11}\,T^{1/2}\,e^{-102124/T}$ carries a
    Boltzmann factor for the 102 124 K (about 8.8 eV) energy barrier: at 10 K it
    is around $10^{-316}$, about 300 orders of magnitude below the
    recombinations. Plotted on the same axes it would squash them flat, so it gets
    its own plot over the range where it matters. It climbs by about 18
    orders of magnitude between 2 000 K and 10 000 K.
    """)
    return


@app.cell
def _(PLOT_STYLE, net, plot_rates):
    _dissociation = net.reactions.get("H2 + e- -> H + H + e-")
    _fig, _ = plot_rates(
        _dissociation,
        tmin=2e3,
        tmax=3e4,
        labels=[_dissociation.verbatim],
        style=PLOT_STYLE,
        show=False,
        close=False,
    )
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Rates over several temperature ranges

    Rate fits are often only accurate over a limited temperature span, so
    databases list the **same reaction several times**, each line with its own
    range. JAFF keeps each line as a `RateSegment` (in `rxn.rate_segments`) and
    merges them into one piecewise rate:

    - inside each range, that range's own fit is used;
    - in a gap between two ranges, the rate is interpolated between the two
      boundary values, linearly in $\log k$ (so rates spanning decades stay
      smooth);
    - outside the outermost ranges the temperature cutoff decides: `clip` holds
      the boundary value, `extrapolate` lets the end fits continue.

    Only one range may be open below (`tmin` is `None`) and only one open above,
    and ranges must not overlap; otherwise loading fails with an error.

    To see this in isolation, the next cell writes a two-line network for
    $\mathrm{H + H \to H_2}$ with ranges 10–100 K and 200–1000 K, in two folders.
    The second folder also gets a `jaff.toml` that switches the cutoff to
    `extrapolate`; `Network` picks up a `jaff.toml` sitting next to the network
    file automatically.
    """)
    return


@app.cell
def _(load_network):
    import tempfile
    from pathlib import Path

    TWO_RANGE = (
        "H + H -> H2  [10, 100]    1e-17*(Tgas/1e2)**0.5\n"
        "H + H -> H2  [200, 1000]  3e-17*(Tgas/3e2)\n"
    )
    EXTRAPOLATE_TOML = '[network.rates]\nT_cutoff = "extrapolate"\n'

    scratch = tempfile.TemporaryDirectory()
    two_range = {}
    for _cutoff in ("clip", "extrapolate"):
        _folder = Path(scratch.name) / _cutoff
        _folder.mkdir()
        (_folder / "two_range.jet").write_text(TWO_RANGE)
        if _cutoff == "extrapolate":
            (_folder / "jaff.toml").write_text(EXTRAPOLATE_TOML)
        two_range[_cutoff] = load_network(_folder / "two_range.jet").reactions[0]
    return (two_range,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    Both networks hold the same two segments:
    """)
    return


@app.cell
def _(mo, two_range):
    mo.ui.table(
        [
            {"tmin": seg.tmin, "tmax": seg.tmax, "rate": str(seg.rate)}
            for seg in two_range["clip"].rate_segments
        ],
        selection=None,
        show_column_summaries=False,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    The merged rate differs only outside 10–1000 K. With `clip` the `Piecewise`
    gains a constant branch at each end; with `extrapolate` the end fits run on.
    In both, the middle branch bridges the 100–200 K gap.
    """)
    return


@app.cell
def _(mo, two_range):
    mo.vstack(
        [
            mo.md("**clip**"),
            two_range["clip"].rate,
            mo.md("**extrapolate**"),
            two_range["extrapolate"].rate,
        ]
    )
    return


@app.cell
def _(PLOT_STYLE, plot_rates, two_range):
    _fig, _ = plot_rates(
        [two_range["clip"], two_range["extrapolate"]],
        tmin=1,
        tmax=1e4,
        labels=["clip", "extrapolate"],
        title="H + H -> H2 over two ranges",
        style=PLOT_STYLE,
        show=False,
        close=False,
    )
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Photo-reactions and cross-sections

    A photo-reaction's rate depends on the radiation field, which the network
    file does not know. JAFF therefore stores it as a placeholder,
    `photorates(index, threshold, ...)`, and keeps the reaction's
    photo-absorption and photodecay **cross-sections** σ(E) instead. When code is
    generated with a radiation field enabled, the placeholder is replaced by an
    integral of σ(E) over the photon spectrum in each radiation band.

    The demo network's two photo-reactions show the placeholder:
    """)
    return


@app.cell
def _(mo, net):
    mo.ui.table(
        [
            {"reaction": rxn.verbatim, "rate": str(rxn.rate)}
            for rxn in net.reactions
            if "_PHOTON" in rxn.reactants.names()
        ],
        selection=None,
        show_column_summaries=False,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    `jaff.plotting.plot_xsecs` draws the cross-sections. Hydrogen photoionisation
    switches on at the 13.6 eV ionisation threshold and falls off steeply above
    it; its data (from the NORAD database) give only the ionisation (photodecay)
    cross-section, so the panel has one curve. For
    CH₂ photodissociation the two differ: photo-absorption counts every absorbed
    photon, photodecay only those that actually break the molecule. Pass
    `energy_unit="nm"` to plot against wavelength instead of photon energy.
    """)
    return


@app.cell
def _(PLOT_STYLE, net, plot_xsecs):
    _fig, _ = plot_xsecs(
        [
            net.reactions.get("H + _PHOTON -> H+ + e-", type="photo"),
            net.reactions.get("CH2 + _PHOTON -> CH + H", type="photo"),
        ],
        grid=(1, 2),
        processes=["photodecay", "photo_absorption"],
        style=PLOT_STYLE,
        show=False,
        close=False,
    )
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    **Next:** `03_comparing_networks.py` compares two networks species by species
    and reaction by reaction.
    """)
    return


@app.cell(hide_code=True)
def _():
    import contextlib
    import io
    import logging
    from pathlib import Path as _Path

    import marimo as mo

    from jaff import Network
    from jaff.plotting import plot_rates, plot_xsecs

    # JAFF resets its logger to INFO on every JaffLogger(); disabling WARNING and
    # below is the only switch that sticks.
    logging.disable(logging.WARNING)

    DEMO1 = mo.notebook_dir().parent / "networks" / "demos" / "demo1.jet"

    # The default "paper" style draws journal-column-sized panels; enlarge them
    # for reading on screen.
    PLOT_STYLE = {"panel_size": (4.0, 3.0)}

    def load_network(path: _Path, **kwargs) -> Network:
        """Load a network, hiding the banner JAFF prints on stdout."""
        with contextlib.redirect_stdout(io.StringIO()):
            return Network(str(path), **kwargs)

    return DEMO1, PLOT_STYLE, load_network, mo, plot_rates, plot_xsecs


if __name__ == "__main__":
    app.run()
