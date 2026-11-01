---
tags:
    - User-guide
    - Elements
---

# Elements

Elements are not something you declare — JAFF **derives** them. When a
[`Network`](network.md) loads, it walks every species' atom list
([`Specie.exploded`](species.md)), collects the distinct atomic symbols, and
exposes them through `Network.elements`. This page is about what
`Network.elements` contains and the composition matrices it builds.

As with [species](species.md) and [reactions](reactions.md), there are **two
layers**:

- a single **`Element`** — one entry from the periodic table, with its mass,
  atomic number, and isotope counts;
- the **`Elements` catalogue** — the sorted, de-duplicated set of elements the
  network actually uses, plus the stoichiometry matrices.

```python
from jaff import Network

net = Network("networks/h_photoionization/h_photo.jet")

net.elements.count  # 1
net.elements.symbols()  # ['H']
```

<!-- prettier-ignore -->
!!! note "Why only `H`?"
    The hydrogen network contains `H`, `H+`, and `e-`, yet a single element
    comes out. Element extraction keeps only **alphabetic** atom tokens, so the
    charge markers (`+`, `-`) and the electron (`e-`) are dropped — they are not
    chemical elements. What remains is the set of real periodic-table symbols.

---

## The two layers

An `Element` knows about _itself_ — "I am oxygen, 8 protons, atomic weight
15.999". It carries no notion of which species it appears in.

The `Elements` catalogue knows about _the set_ — which elements the network
uses, in a fixed sorted order, and how each species is composed from them (the
[composition matrices](#composition-matrices)).

Both classes are **flyweights**: identical instances are reused rather than
rebuilt. Constructing the same element, or an element set over the same species,
hands back the very same object:

```python
net.elements["H"] is net.elements["H"]  # True — same cached Element
```

This keeps element objects unique across the whole program, so identity
(`is`) and equality agree.

---

## The individual `Element`

Every `Element` is a record of the **periodic table**,
`jaff.chemistry.periodic_table()`: JAFF's bundled element table, read once and
shared by every network. It is reference data, not part of a network:

```python
from jaff.chemistry import periodic_table
from jaff.core import Element

periodic_table()  # the whole table: a read-only {symbol: Element} mapping
periodic_table("C")  # one element → Element('C')
Element("C")  # the same object
periodic_table("C").mass  # 1.994473e-23 (grams)
```

Besides the chemical elements, the table holds pseudo-elements such as `e-` and
`GRAIN` (only `name` and `mass` set).

### Attributes

| Attribute     | Type    | Description                                     |
| ------------- | ------- | ----------------------------------------------- |
| `symbol`      | `str`   | Periodic-table symbol (`"H"`, `"C"`, `"O"`)     |
| `name`        | `str`   | Full element name, capitalised (`"Hydrogen"`)   |
| `mass`        | `float` | Standard atomic weight in **grams (CGS)**       |
| `atomic_mass` | `float` | Standard atomic weight in **atomic mass units** |
| `protons`     | `int`   | Atomic number                                   |
| `neutrons`    | `int`   | Neutron count of the most common isotope        |
| `electrons`   | `int`   | Electron count of the neutral atom              |

Note the two mass fields: `mass` is the physical mass in grams (like
[`Specie.mass`](species.md)), while `atomic_mass` is the dimensionless atomic
weight.

```python
h = net.elements["H"]

h.symbol  # 'H'
h.name  # 'Hydrogen'      ← capitalised
h.mass  # 1.673773e-24    ← grams (CGS)
h.atomic_mass  # 1.008           ← amu
h.protons  # 1
h.neutrons  # 0
h.electrons  # 1
```

<!-- prettier-ignore -->
!!! warning "An `Element` has no `index`"
    Unlike a `Specie` or `Reaction`, an `Element` carries no position field. Its
    place in the matrices is its rank in the sorted symbol list — recover it with
    `net.elements.symbols().index("H")`, not `element.index`.

### Comparing and printing an `Element`

An `Element`'s identity for comparison is its **symbol**. Equality, ordering,
and hashing all reduce to that string:

```python
net1 = Network("networks/uclchem_small_gas/uclchem_small_gas_network.jet")
C, H, O = net1.elements["C"], net1.elements["H"], net1.elements["O"]

C == H  # False — different symbols
C == net1.elements["C"]  # True

C < H  # True — '<' / '>' compare symbols lexicographically
H > C  # True
sorted(
    [O, H, C]
)  # [ElementObject(symbol='C'), ElementObject(symbol='H'), ElementObject(symbol='O')]
```

Because the symbol is the hash key (and instances are flyweights), elements work
cleanly in `set`s and as `dict` keys.

<!-- prettier-ignore -->
!!! warning "Comparisons only work element-to-element"
    Comparing an `Element` against a bare string raises `TypeError`:
    ```python
    C == "C"  # TypeError: '==' not supported … 'Element' and 'C'
    C.symbol == "C"  # True
    ```

Printing an element gives its symbol; `repr` wraps it for debugging:

```python
str(net.elements["H"])  # 'H'                        ← __str__ is the symbol
repr(net.elements["H"])  # "ElementObject(symbol='H')"
```

---

## The `Elements` catalogue

`net.elements` is **sorted alphabetically by symbol** and de-duplicated. That
fixed order is what pins the row order of the composition matrices, so it never
shifts between runs.

A single species also has its own element set, reachable via
[`Specie.elements`](species.md):

```python
net1.species["H2O"].elements  # Catalogue(['H', 'O'])
```

### How many elements

`net.elements.count` is the number of unique elements, and the catalogue is
sized, so `len(net.elements)` returns the same value — `count` is the cached
attribute, `len()` the Pythonic spelling:

```python
net.elements.count  # 1
len(net.elements)  # 1   — identical
```

### Finding an element

```python
net1.elements["O"]  # by symbol            → Element
net1.elements[0]  # by index into the set → Element('C')
net1.elements.from_symbol("O")  # same as ["O"]
net1.elements.from_name("Oxygen")  # by full name (capitalised)

"C" in net1.elements  # True — membership is by symbol
```

Lookup is by **symbol** (or full name via `from_name`); there is no `ne`
flag here, because the electron is never an element in the first place.

### Bulk accessors

The hydrogen network has a single element, so switch to a richer one to see
these. Each returns a `Vector` aligned to the sorted symbol order.

```python
net1 = Network("networks/uclchem_small_gas/uclchem_small_gas_network.jet")

net1.elements.symbols()  # ['C', 'H', 'He', 'Mg', 'O']
net1.elements.names()  # ['Carbon', 'Hydrogen', 'Helium', 'Magnesium', 'Oxygen']
net1.elements.atomic_masses()  # [12.011, 1.008, 4.003, 24.305, 15.999]
net1.elements.protons()  # [6, 1, 2, 12, 8]
net1.elements.masses()  # standard atomic weights in grams (CGS)
net1.elements.neutrons()  # neutron counts
net1.elements.electrons()  # electron counts
```

---

## Composition matrices

The point of the `Elements` catalogue is to describe **how every species is
built from elements**. Two matrices do this, both with shape
`(n_elements, n_species)` — element rows, species columns — and both returned as
plain `list[list[int]]` (not NumPy arrays), cached after the first call.

<!-- prettier-ignore -->
!!! warning "Row = element, column = species"
    Index as `matrix[element_row][species_col]`. The row order matches
    `net.elements.symbols()`; the column order matches `net.species` (so the
    column is just `specie.index`). Since an `Element` has no `index`, get the
    row from the symbol list.

### Density matrix

`density_matrix()[i][j]` is the **number of atoms** of element _i_ in species
_j_ — the stoichiometric composition.

```python
net1 = Network("networks/uclchem_small_gas/uclchem_small_gas_network.jet")
density = net1.elements.density_matrix()

syms = list(net1.elements.symbols())  # ['C', 'H', 'He', 'Mg', 'O']
h_row = syms.index("H")  # 1
o_row = syms.index("O")  # 4
h2o = net1.species["H2O"].index  # column for H2O

density[h_row][h2o]  # 2   — two H atoms in H2O
density[o_row][h2o]  # 1   — one O atom
```

### Truth matrix

`truth_matrix()[i][j]` is the **presence mask** — `1` if element _i_ appears in
species _j_, else `0` — i.e. `density_matrix()` clipped to 0/1. Same shape and
indexing as the density matrix.

```python
truth = net1.elements.truth_matrix()
syms = list(net1.elements.symbols())
truth[syms.index("H")][net1.species["H2O"].index]  # 1  — H2O contains H
```

---

## Common patterns

### Print the elemental composition of each species

```python
net1 = Network("networks/uclchem_small_gas/uclchem_small_gas_network.jet")
density = net1.elements.density_matrix()
syms = list(net1.elements.symbols())

print(f"{'Species':<10}", "  ".join(f"{s:>3}" for s in syms))
print("-" * (10 + 5 * len(syms)))
for sp in net1.species:
    col = [density[row][sp.index] for row in range(len(syms))]
    counts = "  ".join(f"{n:>3}" for n in col)
    print(f"{sp.name:<10}  {counts}")
```

### Find every species containing a given element

```python
def species_with_element(net, symbol):
    density = net.elements.density_matrix()
    row = list(net.elements.symbols()).index(symbol)
    return [s for s in net.species if density[row][s.index] > 0]


carbon_species = species_with_element(net1, "C")
print([s.name for s in carbon_species])
```
