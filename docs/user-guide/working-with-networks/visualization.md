---
tags:
    - User-guide
    - Network
---

# Visualising a network

`net.visualize` draws a network's species and reactions as static matplotlib
figures. It shows structure: which species each reaction consumes and
produces. The graph and species views can also colour reactions by rate
coefficient (see [Rate colouring](#rate-colouring)). There are three views:

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
- **Reactions** are small markers labelled with their position in
  `net.reactions`, with the equation written above. The reaction type sets both
  the marker's colour and its shape (square, circle, diamond, triangle,
  hexagon, pentagon, inverted triangle, in order of first appearance), so the
  type stays readable in print and for colour-blind readers. A reaction between species of one column sits beside that column; one
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
- on its left, the reactions **forming** it (numbered markers, shaped and
  coloured by type, equation above), each linked to its other reactants
  further left;
- on its right, the reactions **destroying** it, each linked to its other
  products further right.

A reaction forms the species when it produces more of it than it consumes, and
destroys it when it consumes more; one with the species on both sides and no
net change (a catalyst) is left out. Co-products of forming reactions and
co-reactants of destroying ones appear only in the equations. The legend lists
the reaction types shown, in the colours of the whole network.

---

## Rate colouring

`graph()` and `species(name)` can colour each reaction by its rate coefficient
*k* or its rate 𝓡 (*k* times the reactant densities), so the fast channels stand
out. The values come from [analysis mode](analysis.md): enable it,
then pass `colour_by="rate_coefficient"` or `colour_by="rate"`:

```python
net = Network("GOW")
net.analysis.enable(
    parameters={
        "tgas": {"value": 20.0, "label": r"$T_\mathrm{gas} = {}$ K"},
        "av": {"value": 10.0, "label": r"$A_V$"},
        "chi": 1.0,
        "crate": {"value": 2e-16, "label": r"$\zeta_\mathrm{CR} = {}\ \mathrm{s^{-1}}$"},
        "ncol_CO": 1e18,
        "ncol_H2": 1e22,
    },
    densities={"H": {"value": 1.0, "label": True}, "H2": {"value": 100.0, "label": True}, ...},
    functions={"fshield_co_interp2d": 7.3e-5},
)
net.visualize.species("CO", colour_by="rate_coefficient", save="gow_co.pdf")
```

![CO in GOW: the reactions forming and destroying CO, coloured by rate coefficient](../../assets/figures/working-with-networks/gow_co_rates_light.png#only-light)
![CO in GOW: the reactions forming and destroying CO, coloured by rate coefficient](../../assets/figures/working-with-networks/gow_co_rates_dark.png#only-dark)

- **Values.** Every symbol in the expressions drawn needs a value in
  `net.analysis`: parameters (`tgas`, `av`, `crate`, column densities, ...),
  the densities they use (all reactant densities for 𝓡), and a value for each
  unknown function (one the `.jfunc` file calls but does not define, such as
  GOW's CO shielding table `fshield_co_interp2d`). A missing value raises a
  `ValueError` listing every missing one; `colour_by` with analysis mode off
  raises too. The species view evaluates only the reactions it draws, so it
  needs only their values.
- **Colours.** Each reaction's marker and product arrows take a colour from the
  theme's colormap, logarithmically from the smallest value drawn (lightest) to
  the largest (darkest; brightest in `dark`). Values more than 12 decades below
  the largest get the lightest colour, so one negligible rate (above, the
  shielded CO photodissociation) does not squash the scale. The grey reactant
  lines are unchanged. A rounded colour bar below the legend shows the scale,
  with up to six powers of ten, titled by the quantity.
- **Types.** The colour now shows the rate, so the reaction type is read from
  the marker shape alone; the legend draws the shapes empty.
- **Values legend.** Every analysis entry with a label is listed in a box under
  the colour bar, grouped as parameters, number densities and functions. A `{}`
  in the label takes the value (3 significant figures, `m × 10ⁿ` outside
  0.01–1000); a label without it reads `label = value`; `label=True` gives a
  default such as `n(H₂) = 100`. Entries without a label are used but not
  listed, and with no labelled entry there is no box. A function shows its
  value for the parameters.
- **Colormap.** Each theme has its own (see [Themes](#themes)); pass
  `rate_cmap="magma"` (any matplotlib colormap or its name) to use another, as
  given.

Rate coefficients have units that depend on the reaction order (s⁻¹ for one
reactant, cm³ s⁻¹ for two, ...), so *k* colours compare reactions of the same
order best; rates 𝓡 (cm⁻³ s⁻¹) compare every reaction. The same values evaluate
the coefficients and rates directly, see
[Evaluating rate coefficients and rates](reactions.md#evaluating-rate-coefficients-and-rates).

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

`graph` and `species` also take `colour_by` and `rate_cmap` (see
[Rate colouring](#rate-colouring)).

### Themes

| theme     | reaction types                         | matrix cells                         | rate colormap       |
| --------- | -------------------------------------- | ------------------------------------ | ------------------- |
| `vivid`   | saturated colours (default)            | consumed filled / produced outlined  | `YlGnBu`            |
| `okabe`   | Okabe–Ito (colour-blind safe)          | red / blue                           | `cividis` (reversed) |
| `neutral` | Okabe–Ito on markers and arrowheads    | consumed filled / produced outlined  | `YlGnBu`            |
| `muted`   | soft, desaturated colours              | purple / orange                      | `PuBu`              |
| `print`   | greys only                             | consumed filled / produced outlined  | `Greys`             |
| `dark`    | bright colours on a dark background    | light filled / outlined              | `viridis`           |

A theme's `rates` field (a `RateColormap`: a matplotlib colormap name and the
part of it used, e.g. `RateColormap("YlGnBu", 0.3, 1.0)`) sets its rate
colormap; the light end is trimmed so slow reactions stay visible.

Reaction types get colours and marker shapes in order of first appearance. At
most seven types get their own: reactions of type `unknown`, and any types
beyond the seventh or beyond a theme's palette, share a grey "other" colour and
an octagon.

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
