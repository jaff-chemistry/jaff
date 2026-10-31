---
tags:
    - User-guide
    - Network
---

# Visualising a network

`net.visualize` draws a network's species and reactions as static matplotlib
figures. It shows structure only: which species each reaction consumes and
produces, not rates or physics. There are two views:

```python
from jaff import Network

net = Network("networks/demos/water_co.jet")
net.visualize.graph(save="water_co_graph.pdf")
net.visualize.stoichiometry(save="water_co_matrix.pdf")
```

Both return `(fig, ax)`; see
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

---

## Options

Both methods take the same keywords:

| keyword | default   | meaning                                                                     |
| ------- | --------- | --------------------------------------------------------------------------- |
| `theme` | `"vivid"` | colour theme: a name below or a `NetworkTheme`                              |
| `save`  | `None`    | file to save the figure to; the format follows the suffix                   |
| `show`  | `True`    | show the figure (skipped on non-interactive backends)                       |
| `close` | `True`    | close the figure in pyplot afterwards; the returned figure still works      |

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
`water_co.jet` (20 reactions) is about 10 × 5 in, GOW (50 reactions) about
11 × 10 in, and a network of a few hundred reactions is a poster. Save large
networks as PDF or SVG and zoom in; raster formats beyond matplotlib's
65 536-pixel limit fail to save.
