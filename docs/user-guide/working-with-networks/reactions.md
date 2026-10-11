---
tags:
    - User-guide
    - Reaction
---

# Reactions

When a [`Network`](network.md) is loaded, every line of the network file becomes
a structured object. They are produced by the parser and can be accessed through
`Network.reactions`. This page is about what `Network.reactions` actually
contains and how to use it.

As with [species](species.md), there are **two layers**:

- a single **`Reaction`** — one chemical transformation, carrying its reactants,
  products, symbolic rate, temperature bounds, and energy budget;
- the **`Reactions` Catalogue** — the ordered, doubly-keyed collection that holds
  every `Reaction` in the network and lets you query them in bulk.

`net.reactions` _is_ the Catalogue. Indexing into it hands you one `Reaction`.
The running example throughout is the hydrogen photo-ionization network — two
reactions, one of each interesting kind: a photo-ionization and its inverse
recombination.

```python
from jaff import Network

net = Network("networks/h_photoionization/h_photo.jet")

net.reactions.count            # 2
net.reactions.verbatim()       # ['H + _PHOTON -> H+ + e-', 'H+ + e- -> H']
```

<!-- prettier-ignore -->
!!! note "The `_PHOTON` pseudo-species"
    The photo-ionization reads `H + _PHOTON -> H+ + e-`, not `H -> H+ + e-`.
    `_PHOTON` is a **special pseudo-species** (its name starts with `_`) that
    marks the driving photon; cosmic-ray reactions carry `_CR` likewise. These
    pseudo-species give a photo/CR reaction its distinct identity and
    serialization, but they are **not** real species — they are excluded from
    the kinetics and the integrated ODE state, and never appear in
    `net.species`. See [reaction types](#reaction-types).

---

## The two layers

A `Reaction` knows about _itself_ — which species go in, which come out, how
fast it proceeds, and how much energy it releases. It does not know its place in
the network or how it couples to other reactions; that context belongs to the
catalogue.

The `Reactions` catalogue knows about _the set_ — ordering (each reaction's
`catalogue_index`), two ways to look a reaction up, and how to project any
per-reaction attribute into a flat array for the solver and code generator.

```text
net.reactions                       ← the Reactions catalogue (the set)
  ├── net.reactions[0]  → Reaction  ← one transformation (H + _PHOTON -> H+ + e-)
  └── net.reactions[1]  → Reaction  ← one transformation (H+ + e- -> H)
```

Everything below follows that seam: the individual reaction first, then the
catalogue, with the symbolic rate and the two notions of "same reaction" as the
ideas to hold onto.

---

## The individual `Reaction`

### The rate is symbolic, not a number

The single most important thing about a `Reaction` is that its `rate` is a
**SymPy expression**, not a float. It is a formula in the gas temperature
`tgas` (and possibly other symbols like `crate`, `av`, `ntot`), frozen at load
time but not yet evaluated.

```python
rec = net.reactions[1]      # H+ + e- -> H  (recombination)
rec.rate                    # 1.65941781598291e-10/tgas**0.7
```

Because the rate is symbolic, a single reaction can be turned into many concrete
forms on demand — differentiated, compiled to six languages, or plotted:

```python
import sympy

rec.get_sympy()                          # 1.65941781598291e-10/tgas**0.7
sympy.diff(rec.get_sympy(), "tgas")      # analytic dk/dT, also symbolic
```

A photo-reaction's rate is instead an unevaluated `photorates(...)` call. The
reaction is classified as photochemical by the parser (via its `_PHOTON`
reactant), not by inspecting the rate — see [reaction types](#reaction-types):

```python
net.reactions[0].rate       # photorates(1, 13.6, 1.0e+99)
```

<!-- prettier-ignore -->
!!! note "Photo-reactions during code generation"
    Photo reaction rate during code generation are not kept as photorates if radiation is enabled. Actual photorates are calculated as a function of photo number/energy density for each radiation band

### Attributes

| Attribute             | Type            | Description                                                               |
| --------------------- | --------------- | ------------------------------------------------------------------------- |
| `reactants`           | `Species`       | Catalogue of reactant species (may repeat, e.g. for 3-body reactions)     |
| `products`            | `Species`       | Catalogue of product species                                              |
| `rate`                | `sympy.Basic`   | Symbolic rate-coefficient expression (units depend on reaction order)     |
| `tmin`                | `float or None` | Lower temperature bound of rate validity (K); `None` = unbounded          |
| `tmax`                | `float or None` | Upper temperature bound of rate validity (K); `None` = unbounded          |
| `dE`                  | `sympy.Basic`   | Energy released per reaction event (erg), from a `.jfunc` aux function    |
| `dRad`                | `sympy.Basic`   | Radiation energy emission per photon energy (eV) per event                |
| `verbatim`            | `str`           | Human-readable equation `"R1 + R2 -> P1 + P2"`                            |
| `index`               | `int`           | File-side source number (`source_index`): file-order row number, **gapped** when rows merge. Not an array subscript |
| `catalogue_index`     | `int`           | Dense `0..count-1` position inside `net.reactions`; the array subscript used by generated code |
| `serialized`          | `str`           | Canonical **name-level** identity (isomer-sensitive)                      |
| `serialized_exploded` | `str`           | Canonical **atom-level** identity (isomer-insensitive)                    |
| `type`                | `str`           | Reaction type concluded by the parser (verbatim): gas-phase `"photo"`/`"cosmic_ray"`/`"3_body"`/`"unknown"`, or a grain surface-mechanism type — see [Reaction types](#reaction-types) |
| `custom_rad_rate`     | `bool`          | `True` when the radiation rate came from a `.jfunc`, not cross-sections   |
| `xsecs_dict`          | `XsecsProps or None` | Photo cross-section data for the reaction's single decay channel: `photon_energy` (eV) plus `photo_absorption` and `photodecay` (cm²), `database` (source used) and `photodecay_expr` (Verner symbolic σ(E), else `None`); else `None` |

<!-- prettier-ignore -->
!!! tip "`reactants` and `products` are `Species` Catalogues"
    They are not plain name lists — they are full [`Species`](species.md)
    collections, so every per-species accessor works on them
    (`rxn.reactants.names()`, `rxn.products.charges()`, …). The same species may
    appear more than once, so they are built without the length check.

```python
rxn = net.reactions[0]

rxn.verbatim                 # 'H + _PHOTON -> H+ + e-'
rxn.reactants.names()        # ['H', '_PHOTON']
rxn.products.names()         # ['H+', 'e-']
rxn.tmin, rxn.tmax           # (None, None)
rxn.index                    # 0
```

### Conservation is checked at construction

Every `Reaction` validates mass and charge conservation when it is built — a
violation is logged as a warning (or aborts the load when the network is opened
with `errors=True`). You can re-run the checks yourself:

```python
rxn.check_mass()    # True  — mass balances within one electron mass
rxn.check_charge()  # True  — net charge identical on both sides
```

`check_mass` deliberately tolerates a one-electron-mass discrepancy
(`9.109e-28 g`), so an ionization that "loses" an electron still passes.

---

## Reaction identity: two serialized forms

Just as a [`Specie`](species.md) has a canonical `serialized` identity, so does
a reaction — but a reaction has **two**, and the difference is the whole point.

```python
rxn.serialized            # 'H._PHOTON__H+.e-'    ← name-level  (isomer-sensitive)
rxn.serialized_exploded   # 'H._PHOTON__+/H.e-'   ← atom-level  (isomer-insensitive)
```

Both sort the species on each side, join the species on a side with `.`, and
separate reactants from products with `__`. (The `.` joiner — not `_` — is used
because special pseudo-species names start with `_`.) They differ in _what_ they
sort:

- **`serialized`** uses species **names**. `HCO+` and `HOC+` are different here.
  This is the form used for `==`, hashing, and the catalogue's serialized
  lookup. Two reactions are equal when their `serialized` strings match.
- **`serialized_exploded`** uses each species' atom-level `serialized` form, so
  isomers collapse together. This is what `is_isomer_version` compares to decide
  whether two reactions are the same chemistry written with different isomer
  names.

```python
net.reactions[0].is_isomer_version(net.reactions[1])   # False
```

`is_isomer_version` returns `True` only when the atom-level forms match _and_ at
least one species name differs — i.e. genuine isomer twins, not identical
reactions.

### Comparing and printing a `Reaction`

A reaction's identity for comparison is its **name-level `serialized` form**.
Equality, ordering, and hashing all reduce to that string — so two reactions
written differently but with the same sorted reactants and products are equal,
regardless of `index`, rate, or temperature bounds:

```python
r0, r1 = net.reactions[0], net.reactions[1]

r0 == r1               # False — different serialized forms
r0 == net.reactions[0] # True  — same reaction

r0 < r1                # compares serialized strings, not index order
sorted([r1, r0])       # [ReactionObject(H+ + e- -> H), ReactionObject(H + _PHOTON -> H+ + e-)]
```

<!-- prettier-ignore -->
!!! warning "Ordering is by serialized string, not by index"
    `<` / `>` sort on the `serialized` form, which need not match catalogue
    order. Above, `r1` sorts before `r0` even though its `index` is larger.
    Comparing against a non-`Reaction` (e.g. a string) raises `TypeError`.

Because equality is by `serialized`, reactions can be used in `set`s and as
`dict` keys, and `==` is the basis for the duplicate-reaction check the network
runs at load time.

Printing a reaction gives its human-readable equation; `repr` wraps it:

```python
str(net.reactions[0])    # 'H + _PHOTON -> H+ + e-'            ← __str__ is the verbatim
repr(net.reactions[0])   # 'ReactionObject(H + _PHOTON -> H+ + e-)'
print(net.reactions[0])  # H + _PHOTON -> H+ + e-
```

`str(reaction)` returning the verbatim equation is what lets you drop a
`Reaction` straight into f-strings and log lines.

---

## Reaction types

The `type` attribute holds the type **concluded by the network-format parser**
as it read the file — it does not inspect the rate expression. The parser decides
structurally (a `_PHOTON` reactant → photo, a `_CR`/`_CRP`/`_CRPHOT` reactant →
cosmic-ray, three or more real reactants → 3-body), so the type survives even
custom rates. Gas-grain formats (e.g. UCLCHEM) additionally classify by the
mechanism keyword carried in a reactant slot (`FREEZE`, `THERM`, `DESCR`, …) into
a surface-mechanism type — see the
[UCLCHEM format notes](../designing-networks/network-formats.md#uclchem-format).

Generic gas-phase types:

| Type           | Meaning                                          |
| -------------- | ------------------------------------------------ |
| `"photo"`      | Radiation-driven (photodissociation/ionisation)  |
| `"cosmic_ray"` | Cosmic-ray driven                                |
| `"3_body"`     | Three-body reaction                              |
| `"unknown"`    | Unclassified                                     |

Grain/ice surface-mechanism types (gas-grain networks):

| Type                                | Mechanism                          |
| ----------------------------------- | ---------------------------------- |
| `"freeze"`                          | Freeze-out (accretion) onto grains |
| `"desorption_thermal"`              | Thermal desorption                 |
| `"desorption_cr"`                   | Cosmic-ray desorption              |
| `"desorption_uvcr"`                 | CR-induced UV photodesorption      |
| `"desorption_h2"`                   | H₂-formation desorption            |
| `"eley_rideal"` / `"eley_rideal_desorption"` | Eley-Rideal surface reaction (+ reactive desorption) |
| `"langmuir_hinshelwood"` / `"langmuir_hinshelwood_desorption"` | Langmuir-Hinshelwood diffusion reaction (+ reactive desorption) |
| `"h2_formation"`                    | Parameterised H₂ formation         |
| `"bulk_swap"` / `"surface_swap"`    | Bulk↔surface ice exchange          |

The set is open-ended: any string a parser supplies is stored verbatim.

```python
net.reactions[0].type   # 'photo'
net.reactions[1].type   # 'unknown'
net.reactions.types()     # ['photo', 'unknown']
```

### Temperature ranges and piecewise rates

A reaction is valid over `[tmin, tmax]`. When the **same** reaction (same species
*and* same `type`) appears in the network file across several **disjoint**
temperature ranges — a common way to fit a rate accurately over a wide
temperature span — JAFF keeps each range as a `RateSegment` and merges them into
a single **piecewise** rate expression on `net.reactions[i].rate`:

- each range contributes its own rate over that range;
- the gap between two adjacent ranges is bridged by **linear interpolation** of
  the two bounding rates;
- behaviour past the outermost bounds follows the reaction's
  [`T_cutoff`](jaff-toml.md#temperature-cutoffs) (`clip` holds the boundary
  value; `extrapolate` lets the end segments run on).

The individual pieces live on `rxn.rate_segments`; `rxn.tmin` / `rxn.tmax` span
the full merged range. Two rows for the *same* reaction, same `type`, over the
*same or overlapping* range is a genuine conflict and raises a `ParserError` —
deduplicate the source network.

!!! note "Different mechanisms are different reactions"
    Two file rows with the same species but a different `type` (e.g. thermal vs
    cosmic-ray desorption) are **not** merged — they are distinct reactions with
    their own rates. Only same-species **and** same-`type` rows over disjoint
    ranges are stitched into one piecewise rate.

### Photo-reactions, a special citizen

Photo-reactions are to `Reactions` what the electron is to `Species` — present
in the catalogue like any other, but special-cased throughout. A photo-reaction
carries a cross-section table instead of an analytic rate:

```python
photo = net.reactions[0]

photo.rate                              # photorates(1, 13.6, 1.0e+99)
photo.xsecs_dict.keys()                 # units, _equations, photon_energy,
                                        #   photo_absorption, photodecay,
                                        #   database, photodecay_expr
len(photo.xsecs_dict["photon_energy"])  # number of grid points (energies in eV)
photo.xsecs_dict["photodecay"]          # cross sections in cm^2 (or None)
```

The `photon_energy` grid is in eV and each process array is in cm² (or `None`
when that process has no data for the reaction). The `_equations` sub-dict
carries the boolean `pa` photo-absorption flag and `decay_type`
(`"ionization"` or `"dissociation"`), identifying the single decay channel held
in `photodecay`. `xsecs_dict["database"]` names the source (`"norad"`,
`"verner"` or `"leiden"`); for Verner, `photon_energy` / `photodecay` are
`None` and the symbolic σ(E) (eV, cm²) is in `xsecs_dict["photodecay_expr"]`.

The catalogue gives you dedicated ways to pick them out:

```python
net.reactions.photo_reactions()         # Vector[Reaction] — just the photo ones
net.reactions.photo_reaction_truths()   # [1, 0]  — 1 where photo
net.reactions.photo_reaction_indices()  # [0]     — their positions
```

When generating code, and radiation is disabled, a photo-reaction's rate is emitted with an `$IDX$`
placeholder that the code generator fills in with the real array index later:

```python
photo.get_code(lang="cxx")   # 'photorates($IDX$, 13.6000000000000, 1.0e+99)'
```

---

## The `Reactions` catalogue

`net.reactions` is ordered — the catalogue position (`Reaction.catalogue_index`)
matches the stoichiometry-matrix rows and the generated reaction arrays. Note
this is **not** `Reaction.index`: `index` is the file-side source number, which
is gapped when rows merge (`duplicate_policy` or a temperature-range merge). The
catalogue can be looked up two ways.

### Ways to find a reaction

```python
net.reactions[0]                       # by index        → Reaction
net.reactions[-1]                      # negative index  → last reaction
net.reactions["H + _PHOTON -> H+ + e-"]  # by verbatim string
net.reactions["H._PHOTON__H+.e-"]        # by serialized form
```

String lookup returns a **single** `Reaction`. Because several reactions can
share one serialized form or verbatim string when they differ only by
mechanism/`type` (e.g. thermal vs cosmic-ray desorption), an ambiguous string
key **raises** `KeyError`. Disambiguate with a `(serialized, type)` tuple:

```python
net.reactions["C_DUST__C", "desorption_thermal"]   # exact reaction
net.reactions["C_DUST__C", "desorption_cr"]         # the other mechanism
```

The typed helpers:

```python
net.reactions.from_verbatim("H + _PHOTON -> H+ + e-")      # Reaction or None
net.reactions.from_serialized("H._PHOTON__H+.e-")          # list[Reaction] (all types)
net.reactions.all("H._PHOTON__H+.e-")                      # Vector[Reaction], [] if none
net.reactions.get("H + _PHOTON -> H+ + e-", type="photo")  # Reaction or None; never raises
```

`from_serialized` always returns a **list** (every reaction with that serialized
form); `all` is the same but yields an empty `Vector` instead of raising when the
key is absent. `get` resolves a `(serialized, type)` pair and returns `None`
rather than raising on a miss or ambiguity.

### Iteration and count

`net.reactions.count` is the number of reactions, and the catalogue is sized, so
`len(net.reactions)` returns the same value — `count` is the cached attribute,
`len()` the Pythonic spelling. (The same holds for a single reaction's
`reactants` / `products`, which are `Species` catalogues: `len(rxn.reactants)`
equals `rxn.reactants.count`.)

```python
for rxn in net.reactions:
    print(f"{rxn.index:>3}  {rxn.verbatim:<16}  {rxn.type}")

net.reactions.count   # 2
len(net.reactions)    # 2   — identical
```

### Bulk accessors

Each returns a `Vector` aligned to catalogue order.

| Method                  | Returns                 | h_photo result                       |
| ----------------------- | ----------------------- | ------------------------------------ |
| `verbatim()`            | `Vector[str]`           | `['H + _PHOTON -> H+ + e-', 'H+ + e- -> H']` |
| `types()`              | `Vector[str]`           | `['photo', 'unknown']`               |
| `rates()`               | `Vector[Basic]`         | the two symbolic rate expressions    |
| `reactants()`           | `Vector[Species]`       | one `Species` catalogue per reaction |
| `products()`            | `Vector[Species]`       | one `Species` catalogue per reaction |
| `tmins()` / `tmaxes()`  | `Vector[float or None]` | `[None, None]` / `[None, None]`      |
| `dE()` / `dRad()`       | `Vector[Basic]`         | energy / radiation expressions       |
| `serialized()`          | `Vector[str]`           | `['H._PHOTON__H+.e-', 'H+.e-__H']`   |
| `serialized_exploded()` | `Vector[str]`           | atom-level serialized strings        |

<!-- prettier-ignore -->
!!! warning "`reactants()` / `products()` return `Species`, not name lists"
    Each element is a full [`Species`](species.md) catalogue, not a list of
    strings. Call `.names()` on it if you want the names.

### Filter by type

```python
net.reactions.photo_reactions()              # the photo subset
net.reactions.with_type("cosmic_ray")       # cosmic-ray reactions
net.reactions.with_type("unknown")          # everything unclassified
```

The type key is whatever string the parser stored — the gas-phase set
`"photo"`, `"cosmic_ray"`, `"3_body"`, `"unknown"` (there is no `"CR"`), plus any
grain surface-mechanism type (see [Reaction types](#reaction-types)).

---

## `Reaction` methods

### Species membership

`has_reactant` / `has_product` test that **all** given species are on that side;
`has_any_species` tests for **any** on either side. Each accepts a name, a
`Specie`, or a list.

```python
rec = net.reactions[1]            # H+ + e- -> H

rec.has_reactant("H+")            # True
rec.has_reactant(["H+", "e-"])    # True  — all present
rec.has_product("H")              # True
rec.has_any_species("e-")         # True  — on either side
```

### String representations

```python
rxn = net.reactions[0]

rxn.verbatim                # 'H + _PHOTON -> H+ + e-'   (also rxn.get_verbatim())
rxn.get_latex()             # '${\\rm H} + {\\rm _PHOTON}\\,\\to\\,{\\rm H^{+}} + {\\rm e^{-}}$'
rxn.serialize()             # 'H._PHOTON__H+.e-'
rxn.serialize_exploded()    # 'H._PHOTON__+/H.e-'
```

### Code generation

`get_code` renders the rate expression as source for a target language. Supported
keys: `python`, `c`, `cxx`, `fortran`, `rust`, `julia`, `r`.

```python
rec = net.reactions[1]

rec.get_code(lang="python")    # '1.65941781598291e-10*tgas**(-0.7)'
rec.get_code(lang="cxx")       # '...*std::pow(tgas, -0.7...)'
rec.get_code(lang="fortran")   # '1.65941781598291d-10*tgas**(-0.7d0)'
rec.get_code(lang="julia")     # '1.65941781598291e-10 * tgas .^ (-0.7)'
```

`get_flux_expression` builds the reaction flux `k[i] * y[idx_R1] * y[idx_R2]…`
from each reactant's `fidx`, with configurable variable names and brackets:

```python
rec.get_flux_expression(idx=1)
# 'k[1] * y[idx_hj] * y[idx_e]'

rec.get_flux_expression(idx=1, rate_variable="k",
                        species_variable="nden", brackets="()")
# 'k(1) * nden(idx_hj) * nden(idx_e)'
```

### Plotting

The `Reaction` plot methods are thin wrappers over the
[`jaff.plotting`](../../api/plotting/index.md) free functions `plot_rates` and
`plot_xsecs`, applying the house style. They return the `(fig, ax)` they drew
on, so plots can be composed or saved.

```python
rec.plot_rate_coefficient()         # rate vs temperature (log–log)

photo = net.reactions[0]
photo.plot_xsecs()                              # all processes, overlay, eV vs Mb
photo.plot_xsecs(processes="photodecay")        # one process only
photo.plot_xsecs(layout="subplots")             # one stacked panel per process
photo.plot_xsecs(energy_unit="nm", xsec_unit="cm^2")  # wavelength + cm² axes
photo.plot_xsecs(shade=True, show_bands=True)         # shade + band-averaged bars
photo.plot_xsecs(save=True, filename="h_xsec.pdf")    # write to disk
```

`plot_rate_coefficient` spans `[tmin, tmax]`, defaulting to `2.73 K` and `1e6 K`
when a bound is `None`. `plot_xsecs` is a no-op (returns `None`) for non-photo
reactions (those with `xsecs_dict is None`) or when no requested process has
data. `plot_rate_coefficient` returns `None` for photo reactions, whose rate
carries a symbolic radiation-density variable that cannot be evaluated against
temperature.

#### Comparing several reactions

To overlay multiple reactions on one axes, call the free functions (or the
`Reactions` catalogue methods) with a list:

```python
from jaff.plotting import plot_rates, plot_xsecs

plot_rates(list(net.reactions))                 # every rate on shared axes
net.reactions.plot_rates()                      # equivalent, via the catalogue

# The free functions accept any list of reactions (photo_reactions() returns one).
plot_xsecs(net.reactions.photo_reactions(), show_bands=True)
net.reactions.plot_xsecs()                      # all reactions; non-photo skipped
```

Each curve gets a distinct line width (thinner in front) and a legend entry.
See the [`jaff.plotting`](../../api/plotting/index.md) reference for the full
options, palettes, and theming.

---

## Common patterns

### Conservation audit

```python
bad = [r.verbatim for r in net.reactions
       if not (r.check_mass() and r.check_charge())]

if bad:
    print(f"Conservation failures ({len(bad)}):")
    for v in bad:
        print(f"  {v}")
else:
    print("All reactions conserve mass and charge.")
```

### Formation and destruction pathways

```python
def pathways(net, species_name):
    formed    = [r for r in net.reactions if r.has_product(species_name)]
    destroyed = [r for r in net.reactions if r.has_reactant(species_name)]
    print(f"{species_name}: {len(formed)} formation / {len(destroyed)} destruction")
    for r in formed:
        print(f"  +  {r.verbatim}")
    for r in destroyed:
        print(f"  -  {r.verbatim}")

pathways(net, "H+")
```

### Group reactions by type

```python
from collections import Counter

Counter(net.reactions.types())     # Counter({'photo': 1, 'unknown': 1})
```

### Export to CSV

```python
import csv

with open("reactions.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["index", "reaction", "type", "tmin", "tmax",
                "n_reactants", "n_products"])
    for rxn in net.reactions:
        w.writerow([rxn.index, rxn.verbatim, rxn.type,
                    rxn.tmin, rxn.tmax,
                    len(rxn.reactants), len(rxn.products)])
```
