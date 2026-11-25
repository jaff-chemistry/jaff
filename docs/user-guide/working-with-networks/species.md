---
tags:
    - User-guide
    - Species
---

# Species

When a [`Network`](network.md) is loaded, every chemical entity it mentions is
turned into a structured object. They are
produced by the parser and can be accessed through `net.species`. This page is
about what `net.species` actually contains and how to use it.

Special pseudo-species — names starting with `_`, such as `_PHOTON` or `_CR` —
are also built as `Specie` objects (with `index` `-1`), but they only appear in
the reactions that mention them; they are **not** added to `net.species`. That is
why `h_photo`, whose photo-ionization reaction reads `H + _PHOTON -> H+ + e-`,
has three species, not four.

There are **two layers**, and keeping them apart is the key to the whole API:

- A single **`Specie`** — one atom, molecule, or ion, with its mass, charge,
  composition, and the identifiers JAFF needs for code generation;
- the **`Species` Catalogue** — the ordered, multi-keyed collection that holds
  every `Specie` in the network and lets you query them in bulk.

`net.species` _is_ the `Catalogue`. Indexing into it (`net.species["H"]`) hands
you one `Specie` object. Throughout this page the running example is the hydrogen
photo-ionization network, which is small enough to print in full yet still
contains a neutral atom, a cation, and the electron:

```python
from jaff import Network

net = Network("networks/h_photoionization/h_photo.jet")

net.species.count  # 3
net.species.names()  # ['H', 'H+', 'e-']
```

---

## The two layers

A `Specie` knows about _itself_ — "I am H⁺, my mass is 1.67e-24 g, my charge is
+1". It has no idea what reactions it takes part in or where it sits relative
to other species; that context belongs to the catalogue.

The `Species` catalogue knows about _the set_ — ordering (each specie's
`index`), how to find a specie three different ways, and how to project any
per-species attribute into a flat array for the solver and code generator. It
holds `Specie` objects but adds no chemistry of its own.

```text
net.species                      ← the Species catalogue (the set)
  ├── net.species[0]   → Specie  ← one entity (H)
  ├── net.species[1]   → Specie  ← one entity (H+)
  └── net.species[2]   → Specie  ← one entity (e-)
```

Everything below is split along this seam: first the individual specie, then
the catalogue, then the one specie that bends both — the electron.

---

## The individual `Specie`

Parsing a name like `"H2O+"` is the single most important thing a `Specie`
does. From that one string it derives the composition, the mass, the charge,
and the identifiers used downstream. Once built, a `Specie` is just a bundle of
those derived facts.

### Attributes

| Attribute    | Type        | Description                                                                  |
| ------------ | ----------- | ---------------------------------------------------------------------------- |
| `name`       | `str`       | Chemical name exactly as written in the network file (`"H"`, `"H+"`, `"e-"`) |
| `mass`       | `float`     | Total mass in **grams (CGS)**, summed over the constituent atoms             |
| `charge`     | `int`       | Net charge in elementary-charge units (`0` neutral, `>0` cation, `<0` anion) |
| `index`      | `int`       | Zero-based position of this specie inside `net.species`                      |
| `exploded`   | `list[str]` | Sorted list of atomic symbols with repetition; the charge token is included  |
| `serialized` | `str`       | Canonical identity string: `exploded` sorted and joined with `/`             |
| `fidx`       | `str`       | Code-safe flat identifier used in generated C/Fortran/Python source          |
| `elements`   | `Elements`  | Lazily built [`Elements`](elements.md) collection for this single specie     |
| `is_special` | `bool`      | `True` for a pseudo-species (name starts with `_`, e.g. `_PHOTON`, `_CR`)    |
| `is_core`    | `bool`      | `not is_special` — a real species taking part in the kinetics                |

<!-- prettier-ignore -->
!!! tip "Specie latex representation"
    `latex` is a **method**, not an attribute — call it (`specie.latex()`).

```python
h = net.species["H"]

h.name  # 'H'
h.mass  # 1.673773e-24   ← grams, not amu
h.charge  # 0
h.index  # 0
h.exploded  # ['H']
h.serialized  # 'H'
h.fidx  # 'idx_h'
h.latex()  # '{\\rm H}'
```

!!! warning "Mass is in grams (CGS), not atomic mass units"
`specie.mass` is the physical mass in grams — `H` is `1.673773e-24`, not
`1.008`. JAFF works in CGS so the value drops straight into rate and energy
expressions. If you want the atomic weight in amu, read it from the element
data instead.

### How a name becomes composition

The interesting attributes are the ones the parser _derives_. The same H⁺
specie shows all of them at once:

```python
hp = net.species["H+"]

hp.charge  # 1          ← counted from trailing '+' signs
hp.exploded  # ['+', 'H'] ← atoms (sorted), charge carried as a token
hp.serialized  # '+/H'      ← canonical identity
hp.fidx  # 'idx_hj'   ← '+' → 'j' so it's a legal identifier
hp.latex()  # '{\\rm H^{+}}'
```

- **`exploded`** is the formula flattened into individual atoms and sorted, so
  `H2O` becomes `['H', 'H', 'O']`. The net charge rides along as a `'+'` or
  `'-'` token, which is what lets the serialized form encode charge too.
- **`serialized`** is the canonical identity. Because it is _sorted_, structural
  isomers collapse onto the same string — `HCO+` and `HOC+` both serialize to
  `'+/C/H/O'`. The catalogue uses this as an alternate lookup key (see below),
  and the network uses it to detect isomers and duplicate reactions.
- **`charge`** is read only from `+`/`-` characters at the _end_ of the name.
- **`fidx`** is the name made safe for generated source: lower-cased, `'+' → 'j'`,
  `'-' → 'k'`, so `H2O+` becomes `idx_h2oj`. This is the symbol the code
  generator emits to index this specie.
- **`latex()`** renders the name for plots and tables — subscripted counts,
  superscripted charge, roman element font. Pass `dollars=True` to wrap it in
  `$…$` math delimiters.

### Comparing and printing a `Specie`

A `Specie`'s identity for comparison is its **name**. Equality, ordering, and
hashing all reduce to the name string:

```python
H, Hp, e = net.species["H"], net.species["H+"], net.species["e-"]

H == Hp  # False — different names
H == net.species[0]  # True  — same specie, compared by name

H < Hp  # True  — '<' / '>' compare names lexicographically
Hp > H  # True

sorted([Hp, e, H])  # [SpecieObject('H'), SpecieObject('H+'), SpecieObject('e-')]
```

Because the name is also the hash key, species can be used directly in `set`s
and as `dict` keys. Sorting is by raw string order, so capitals sort before
lowercase (`'H'` < `'H+'` < `'e-'`).

<!-- prettier-ignore -->
!!! warning "Comparisons only work specie-to-specie"
    Comparing a `Specie` against a bare string raises `TypeError` — use the name
    explicitly instead:
    ```python
    H == "H"  # TypeError: '==' not supported … 'Specie' and 'H'
    H.name == "H"  # True
    ```

Printing a specie gives back its chemical name; `repr` wraps it for debugging:

```python
str(net.species["H+"])  # 'H+'             ← __str__ is just the name
repr(net.species["H+"])  # "SpecieObject('H+')"
print(net.species["H+"])  # H+
```

`str(specie)` returning the plain name is what lets you drop a `Specie` straight
into f-strings and joins without touching `.name`.

---

## The electron, a special citizen

The electron is a species like any other in the catalogue, but it breaks enough
rules to deserve its own section — and the catalogue API is built around that.

**It must be named `e-`.** The parser rejects `E`, `E-`, and (in any case)
`e`, `el`, `els`, `electrons` and the misspelling `eletron` with a fatal error.
Other spellings are not caught by that check: `electron`, for example, is
parsed as an ordinary formula and fails with `KeyError: 'e'`. Either way there
is exactly one spelling for the electron everywhere in JAFF.

**It does not decompose into atoms.** Where `H+` explodes to `['+', 'H']`, the
electron stays whole:

```python
e = net.species["e-"]

e.charge  # -1
e.mass  # 9.10938371e-28 ← the electron mass, in grams
e.exploded  # ['e-']         ← not split, not an atom token
e.serialized  # 'e-'
e.fidx  # 'idx_e'        ← special-cased, not 'idx_ek'
```

**It is often handled separately by solvers**, so almost every catalogue
accessor takes an `ne` ("no electron") flag. Set `ne=True` and the electron
drops out of the returned vector while every other specie keeps its place:

```python
net.species.names()  # ['H', 'H+', 'e-']
net.species.names(ne=True)  # ['H', 'H+']

net.species.charges()  # [0, 1, -1]
net.species.charges(ne=True)  # [0, 1]
```

When you need the electron's slot directly, ask the catalogue:

```python
net.species.e_idx()  # 2   (None if the network has no electron)
```

---

## The `Species` Catalogue

`net.species` is an ordered collection that you can query three ways and slice
into bulk arrays. The ordering is fixed at load time and matches every
`Specie.index`, the stoichiometry matrices, and the generated code.

### Three ways to find a specie

```python
net.species["H"]  # by name        → SpecieObject('H')
net.species["+/H"]  # by serialized  → SpecieObject('H+')
net.species[1]  # by index       → SpecieObject('H+')
net.species[-1]  # negative index → SpecieObject('e-')
```

Lookup by serialized form is what makes isomer-insensitive queries possible: any
arrangement of the same atoms and charge resolves to the same specie.

### How many species

`net.species.count` is the number of species in the catalogue, and the catalogue
is sized — so `len(net.species)` returns the same value. `count` is a cached
attribute kept in step as species are added; `len()` is the Pythonic spelling.
Use whichever reads better:

```python
net.species.count  # 3
len(net.species)  # 3   — identical
```

### Iteration and membership

Iterating yields each `Specie` in index order. Membership accepts a name _or_ a
serialized string:

```python
for sp in net.species:
    print(f"{sp.index:>3}  {sp.name:<4}  {sp.mass:.3e} g  charge {sp.charge:+d}")

"H" in net.species  # True   (by name)
"+/H" in net.species  # True   (by serialized form)
```

### Bulk accessors

Each method returns a `Vector` — a flat, index-aligned sequence ready to feed
the solver or code generator. Most take the `ne` flag described in the electron
section; `normalized_names()`, `e_idx()` (and the `count` attribute) do not.

| Method               | Returns             | h_photo result                        |
| -------------------- | ------------------- | ------------------------------------- |
| `count`              | `int` (attribute)   | `3`                                   |
| `names()`            | `Vector[str]`       | `['H', 'H+', 'e-']`                   |
| `normalized_names()` | `Vector[str]`       | `['h', 'hj', 'ek']`                   |
| `masses()`           | `Vector[float]`     | `[1.674e-24, 1.673e-24, 9.109e-28]`   |
| `charges()`          | `Vector[int]`       | `[0, 1, -1]`                          |
| `charge_truths()`    | `Vector[int]`       | `[0, 1, 1]`                           |
| `serialized()`       | `Vector[str]`       | `['H', '+/H', 'e-']`                  |
| `exploded()`         | `Vector[list[str]]` | `[['H'], ['+', 'H'], ['e-']]`         |
| `latex()`            | `Vector[str]`       | `['${\\rm H}$', '${\\rm H^{+}}$', …]` |
| `elements()`         | `Vector[Elements]`  | one `Elements` collection per specie  |
| `e_idx()`            | `int or None`       | `2`                                   |

A few are easy to misread:

- **`normalized_names()`** makes each name a legal lowercase identifier:
  `'+' → 'j'`, `'-' → 'k'`. So `H+` becomes `'hj'` and `e-` becomes `'ek'`.
  Pass `pos`/`neg` to override the replacement strings. This uses the same
  `j`/`k` scheme as `fidx`, so the two agree for ordinary species; they differ
  only for the electron, which `normalized_names()` gives as `'ek'` whereas
  `get_fidx` special-cases it to `idx_e`.
- **`charge_truths()`** is a 0/1 mask, `1` where the specie is charged — useful
  for charge-conservation terms.
- **`latex()`** on the catalogue defaults to `dollars=True` (wrapped in `$…$`),
  whereas `Specie.latex()` defaults to `dollars=False`.

### Filtering by charge

`neutral()` and `charged()` return the matching `Specie` objects, or — if you
pass an attribute name — that attribute projected out of each match. `charged()`
also honours the `ne` flag.

```python
net.species.neutral()  # [SpecieObject('H')]
net.species.neutral("name")  # ['H']

net.species.charged()  # [SpecieObject('H+'), SpecieObject('e-')]
net.species.charged("name")  # ['H+', 'e-']
net.species.charged("index")  # [1, 2]
net.species.charged("name", ne=True)  # ['H+']   (electron excluded)
```

The attribute name must be a real `Specie` attribute (`name`, `mass`, `charge`,
`index`, `exploded`, `serialized`, `fidx`, `elements`); anything else raises
`ValueError`.

---

## Common patterns

### Print a formatted summary

```python
print(f"{'Idx':>3}  {'Name':<5}  {'Mass (g)':>12}  {'Charge':>7}")
print("-" * 34)
for sp in net.species:
    print(f"{sp.index:>3}  {sp.name:<5}  {sp.mass:>12.3e}  {sp.charge:>+7d}")
```

### Find every ion

```python
ions = net.species.charged("name", ne=True)  # exclude the electron
print(f"{len(ions)} ions: {list(ions)}")
```

### Find hydrogen-bearing species

`exploded` makes "contains an H atom" an exact test rather than a substring
guess (`"H"` also appears in `"He"`):

```python
h_bearing = [sp for sp in net.species if "H" in sp.exploded]
```

### Mass-sorted list

```python
by_mass = sorted(net.species, key=lambda s: s.mass)
for s in by_mass:
    print(f"{s.name}: {s.mass:.3e} g")
```
