---
tags:
    - User-guide
    - Network
---

# Visualising a network

`net.visualize` draws a network's species and reactions as static matplotlib
figures. It shows structure: which species each reaction consumes and
produces. The graph and species views can also scale their arrows by rate
coefficient (see [Rate-scaled arrows](#rate-scaled-arrows)). There are three
views:

```python
from jaff import Network

net = Network("networks/demos/water_co.jet")
net.visualize.graph(save="water_co_graph.pdf")
net.visualize.stoichiometry(save="water_co_matrix.pdf")
net.visualize.species("H2O", save="water.pdf")
```

All return `(fig, ax)`; see
[`NetworkVisualizer`](../../reference/plotting/networkvisualizer/index.md) for
the full reference.

---

## The graph

`graph()` lays the network out as a layered graph:

- **Species** stand in columns by number of atoms, lightest on the left
  ("0–1 ATOMS", "2 ATOMS", ...): reading left to right follows molecules being
  built up. Electrons and other atom-free species share the first column.
- **Reactions** are small squares labelled with their position in
  `net.reactions` and coloured by reaction type, with the equation written
  above. A reaction between species of one column sits beside that column; one
  linking columns sits between them. Crowded gutters split into several
  sub-columns, reactions reaching further right placed further right.
- **Edges**: thin grey lines lead from each reactant into its reaction; coloured
  arrows lead out to the products. An edge that has to pass other columns runs
  in the gap under a row, never through an equation.

Each species sits level with the reactions it takes part in, so most edges are
short.

## The stoichiometry matrix

`stoichiometry()` draws the net change `product_matrix - reactant_matrix`: one
row per reaction (its equation, then its number in a circle coloured by type)
and one column per species, grouped by number of atoms. A cell shows how many of
the species the reaction consumes or produces; colour tells which (in the
default theme consumed cells are filled, produced cells outlined). A species
that appears on both sides with no net change gets no cell.

## One species

`species(name)` draws the reactions that form and destroy one species, which
stays readable in networks of any size:

- the species sits in the middle, outlined in bold;
- on its left, the reactions **forming** it (numbered squares coloured by type,
  equation above), each linked to its other reactants further left;
- on its right, the reactions **destroying** it, each linked to its other
  products further right.

A reaction forms the species when it produces more of it than it consumes, and
destroys it when it consumes more; one with the species on both sides and no
net change (a catalyst) is left out. Co-products of forming reactions and
co-reactants of destroying ones appear only in the equations. The legend lists
the reaction types shown, in the colours of the whole network.

---

## Rate-scaled arrows

`graph()` and `species(name)` can draw each reaction's product arrows with a
width set by its rate coefficient, so the fast channels stand out. Pass the
values the rate coefficients need as `rate_values`:

```python
import numpy as np

net = Network("GOW")
nden = np.zeros(len(net.species))
nden[net.species["H"].index] = 1.0
nden[net.species["H2"].index] = 100.0

net.visualize.species(
    "CO",
    rate_values={
        "tgas": 20.0, "av": 10.0, "chi": 1.0, "crate": 2e-16,
        "ncol_CO": 1e18, "ncol_H2": 1e22, "nden": nden,
    },
    rate_functions={"fshield_co_interp2d": 7.3e-5},
    save="gow_co.pdf",
)
```

![CO in GOW: the reactions forming and destroying CO, product arrows scaled by rate coefficient](../../assets/figures/working-with-networks/gow_co_rates_light.png#only-light)
![CO in GOW: the reactions forming and destroying CO, product arrows scaled by rate coefficient](../../assets/figures/working-with-networks/gow_co_rates_dark.png#only-dark)

- **Values.** Every symbol in the rate coefficients drawn needs a value:
  `tgas`, `av`, `crate`, column densities, ... Densities `nden[i]` take one
  array under `"nden"`, in species order. An unknown function (one the `.jfunc`
  file calls but does not define, such as GOW's CO shielding table
  `fshield_co_interp2d`) needs an entry in `rate_functions`: a constant, or a
  callable taking the function's numeric arguments. A missing name raises a
  `ValueError` listing every missing one. The species view evaluates only the
  reactions it draws, so it needs only their values.
- **Widths.** Arrow widths go logarithmically from the smallest rate drawn
  (thinnest) to the largest (thickest), between the two widths of `arrow_width`
  (`(0.4, 3.5)` points by default). Rates more than 12 decades below the largest
  get the thinnest width, so one negligible rate (above, the shielded CO
  photodissociation) does not squash the scale. The grey reactant lines keep
  their fixed width.
- **Key.** A row below the legend shows up to three reference widths, at powers
  of ten within the range drawn.

Rate coefficients have units that depend on the reaction order (s⁻¹ for one
reactant, cm³ s⁻¹ for two, ...), so widths compare reactions of the same order
best. The same values evaluate the coefficients directly, see
[Evaluating a rate coefficient](reactions.md#evaluating-a-rate-coefficient).

---

## Options

All three take the same keywords (`species` takes the species name first):

| keyword | default   | meaning                                                                     |
| ------- | --------- | --------------------------------------------------------------------------- |
| `theme` | `"vivid"` | colour theme: a name below or a `NetworkTheme`                              |
| `save`  | `None`    | file to save the figure to; the format follows the suffix                   |
| `show`  | `True`    | show the figure (skipped on non-interactive backends)                       |
| `close` | `True`    | close the figure in pyplot afterwards; the returned figure still works      |
| `tight` | `False`   | crop the figure to its drawn content, with no padding around it             |

`graph` and `species` also take `rate_values`, `rate_functions` and
`arrow_width` (see [Rate-scaled arrows](#rate-scaled-arrows)).

### Themes

| theme     | reaction types                         | matrix cells                         |
| --------- | -------------------------------------- | ------------------------------------ |
| `vivid`   | saturated colours (default)            | consumed filled / produced outlined  |
| `okabe`   | Okabe–Ito (colour-blind safe)          | red / blue                           |
| `neutral` | Okabe–Ito on squares and arrowheads    | consumed filled / produced outlined  |
| `muted`   | soft, desaturated colours              | purple / orange                      |
| `print`   | greys only                             | consumed filled / produced outlined  |
| `dark`    | bright colours on a dark background    | light filled / outlined              |

Reaction types get colours in order of first appearance. Reactions of type
`unknown`, and any types beyond a theme's palette (`vivid` has ten colours),
share a grey "other" colour.

```python
net.visualize.graph(theme="print", save="water_co_graph.pdf", show=False)
```

---

## Size

Spacing is fixed, so labels never shrink: the figure grows with the network.
A legend wider than the drawing wraps onto more rows, in balanced columns, so
it never makes the figure wider.
`water_co.jet` (20 reactions) is about 10 × 5 in, GOW (50 reactions) about
11 × 10 in, and a network of a few hundred reactions is a poster. Save large
networks as PDF or SVG and zoom in; raster formats beyond matplotlib's
65 536-pixel limit fail to save.
