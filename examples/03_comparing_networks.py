# ABOUTME: Example notebook: compare two networks species-by-species and
# ABOUTME: reaction-by-reaction, including rates of shared reactions (marimo)

import marimo

__generated_with = "0.24.2"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    # 03 · Comparing networks

    You often need to know how two networks differ: when switching to a newer
    release of a rate database, when checking a reduced network against the full
    one it was cut from, or when merging chemistry from two sources. This
    notebook compares `demo1` (from the previous notebooks) with `demo2`, a
    two-reaction network, and shows that matching reactions is not the same as
    matching rates.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## The two networks

    `demo2.jet` contains just two PRIZMO lines. The `label` argument of
    `Network` gives a network a readable name; without it the label would be the
    file name.
    """)
    return


@app.cell
def _(DEMOS, mo):
    mo.md("```text\n" + (DEMOS / "demo2.jet").read_text() + "\n```")
    return


@app.cell
def _(DEMOS, load_network):
    net1 = load_network(DEMOS / "demo1.jet")
    net2 = load_network(DEMOS / "demo2.jet", label="very small network")
    return net1, net2


@app.cell
def _(mo, net1, net2):
    mo.ui.table(
        [
            {
                "label": net.label,
                "species": net.species.count,
                "reactions": net.reactions.count,
            }
            for net in (net1, net2)
        ],
        selection=None,
        show_column_summaries=False,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## What counts as "the same"

    JAFF matches species and reactions by their **serialized** identity, never by
    position in the file:

    - a species' `serialized` form is its canonical name, so `E` and `e-` are the
      same species;
    - a reaction's `serialized` form sorts the species names on each side, so the
      order the file lists them in does not matter. `demo1` writes
      `CO+ + E -> CO` while `demo2` writes `E + CO+ -> CO`; both serialize to
      the same string and are equal.

    The verbatim form, by contrast, keeps the file's order: `demo2` stores this
    reaction as `e- + CO+ -> CO`. `net.reactions.get()` accepts either form, but
    only the serialized one finds a reaction regardless of how a file wrote it,
    so it is the right key for looking up a reaction across networks.

    The serialized form ignores the rate and the temperature range. Two
    networks can share a reaction but disagree on how fast it goes; we come back
    to this at the end.
    """)
    return


@app.cell
def _(mo, net1, net2):
    _r1 = net1.reactions.get("CO+ + e- -> CO")
    _r2 = net2.reactions.get("e- + CO+ -> CO")
    mo.ui.table(
        [
            {
                "network": net1.label,
                "original line": _r1.original_string,
                "serialized": _r1.serialized,
            },
            {
                "network": net2.label,
                "original line": _r2.original_string,
                "serialized": _r2.serialized,
            },
            {"network": "equal?", "original line": "", "serialized": str(_r1 == _r2)},
        ],
        selection=None,
        show_column_summaries=False,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Species

    Comparing two networks is then a set operation on serialized forms: the
    intersection is what they share, the differences are what each one lacks.
    The table lists every species that appears in either network, with one
    column per network marking whether it contains that species.
    """)
    return


@app.cell
def _(compare, mo, net1, net2):
    species_cmp = compare(
        {sp.serialized: sp.name for sp in net1.species},
        {sp.serialized: sp.name for sp in net2.species},
        net1.label,
        net2.label,
    )
    mo.ui.table(species_cmp, selection=None, page_size=20)
    return (species_cmp,)


@app.cell
def _(mo, net1, net2, species_cmp, summary):
    mo.md(summary(species_cmp, net1.label, net2.label, "species"))
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    `demo2` is a strict subset of `demo1`: everything it contains also appears in
    `demo1`, which has eleven species `demo2` lacks.

    ## Reactions

    The same comparison for reactions:
    """)
    return


@app.cell
def _(compare, mo, net1, net2):
    reaction_cmp = compare(
        {rxn.serialized: rxn.verbatim for rxn in net1.reactions},
        {rxn.serialized: rxn.verbatim for rxn in net2.reactions},
        net1.label,
        net2.label,
    )
    mo.ui.table(reaction_cmp, selection=None, page_size=20)
    return (reaction_cmp,)


@app.cell
def _(mo, net1, net2, reaction_cmp, summary):
    mo.md(summary(reaction_cmp, net1.label, net2.label, "reactions"))
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Shared reactions, different rates

    Both networks contain `CO+ + e- -> CO` and `CO -> CO+ + e-`, so a
    species-level comparison calls them identical. Their rates tell a different
    story. `demo1` uses constant rates; `demo2` makes the recombination depend on
    temperature and the ionisation depend on visual extinction $A_V$ (it is
    driven by starlight, which dust attenuates).

    The table evaluates both at $T = 100$ K and $A_V = 1$. The recombination
    rates differ by more than two orders of magnitude, a difference that changes
    the CO⁺ abundance a model predicts. When comparing networks, always check
    the rates of the reactions they share.
    """)
    return


@app.cell
def _(mo, net1, net2, reaction_cmp):
    CONDITIONS = {"tgas": 100.0, "av": 1.0}

    _rows = []
    for _row in reaction_cmp:
        if not (_row[net1.label] and _row[net2.label]):
            continue
        _r1 = net1.reactions.get(_row["serialized"])
        _r2 = net2.reactions.get(_row["serialized"])
        _k1 = float(_r1.rate.subs(CONDITIONS))
        _k2 = float(_r2.rate.subs(CONDITIONS))
        _rows.append(
            {
                "reaction": _row["name"],
                f"rate in {net1.label}": str(_r1.rate),
                f"rate in {net2.label}": str(_r2.rate),
                f"k ({net1.label})": f"{_k1:.3e}",
                f"k ({net2.label})": f"{_k2:.3e}",
                "ratio": f"{_k2 / _k1:.3g}",
            }
        )
    mo.ui.table(_rows, selection=None, show_column_summaries=False)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Built-in helpers

    `Network.compare_species(other)` and `Network.compare_reactions(other)` run
    the same set comparison and report it on the console: the lists of species or
    reactions only in one network, only in the other, and in both, followed by
    the counts. They are meant for quick checks in a terminal session; they
    return nothing, so in a notebook the tables above are more useful.

    **Next:** `04_codegen_jaffgen.py` turns a network into solver source code
    with the `jaffgen` command-line tool.
    """)
    return


@app.cell(hide_code=True)
def _():
    import contextlib
    import io
    import logging
    from pathlib import Path
    from typing import Dict, List

    import marimo as mo

    from jaff import Network

    # JAFF resets its logger to INFO on every JaffLogger(); disabling WARNING and
    # below is the only switch that sticks.
    logging.disable(logging.WARNING)

    DEMOS = mo.notebook_dir().parent / "networks" / "demos"

    def load_network(path: Path, **kwargs) -> Network:
        """Load a network, hiding the banner JAFF prints on stdout."""
        with contextlib.redirect_stdout(io.StringIO()):
            return Network(str(path), **kwargs)

    def compare(
        first: Dict[str, str], second: Dict[str, str], label1: str, label2: str
    ) -> List[Dict]:
        """One row per serialized key in either mapping, flagging membership."""
        names = {**second, **first}
        return [
            {
                "name": names[key],
                "serialized": key,
                label1: key in first,
                label2: key in second,
            }
            for key in sorted(names, key=names.get)
        ]

    def summary(rows: List[Dict], label1: str, label2: str, noun: str) -> str:
        """Count shared and one-sided rows of a :func:`compare` table."""
        both = sum(r[label1] and r[label2] for r in rows)
        only1 = sum(r[label1] and not r[label2] for r in rows)
        only2 = sum(r[label2] and not r[label1] for r in rows)
        return (
            f"**{both}** {noun} in both · **{only1}** only in `{label1}` · "
            f"**{only2}** only in `{label2}`"
        )

    return DEMOS, compare, load_network, mo, summary


if __name__ == "__main__":
    app.run()
