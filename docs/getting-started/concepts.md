---
tags:
    - Introduction
icon: phosphor/stack
---

# Basic Concepts

JAFF does two things: it gives every common chemical-network format a
single in-memory representation, and it turns that representation into
source code in the language of your choice — with
first-class support for explicit photochemistry.

This page walks the objects you will meet, in the order you meet them. The
running example throughout is the hydrogen photo-ionization network, the same
one used across the rest of the user guide. It is small enough to print in full
yet still contains a neutral atom, a cation, the electron, and a genuine
photo-reaction:

```text
H + _PHOTON -> H+ + e-
H+ + e- -> H
```

Two reactions, three species, and rate coefficients that depend on temperature
or other parameters. A rate coefficient tells you how fast a reaction runs. For a reaction

$$ \alpha A + \beta B \rightarrow \gamma C $$

with stoichiometric coefficients $\alpha$, $\beta$, $\gamma$, the rate is

$$ r = k\ [A]^{\alpha} [B]^{\beta} $$

where $k$ is the rate coefficient and $[A]$, $[B]$ are the number densities of
$A$ and $B$. In astrophysics a reaction can fire for many reasons — thermal
collisions, cosmic-ray impacts, photons from a nearby source, or spontaneous
decay — and each cause carries its own rate law. Keeping those causes
distinguishable is one of the things JAFF is built to do.

---

## The shape of the library

There are **two halves**, and keeping them apart is the key to the whole API:

- the **model** — a loaded [`Network`](../user-guide/working-with-networks/network.md)
  and the typed objects it holds: [`Species`](../user-guide/working-with-networks/species.md),
  [`Reactions`](../user-guide/working-with-networks/reactions.md), and
  [`Elements`](../user-guide/working-with-networks/elements.md);
- the **code generator** — the
  [`jaffgen`](../user-guide/code-generation/index.md) CLI, which expands your
  `$$JAFF` templates against the model to emit rate, ODE, and Jacobian source.

```text
Network                          ← the loaded model
  ├── net.species   → Species    ← the chemical entities
  ├── net.reactions → Reactions  ← the transformations
  └── net.elements  → Elements   ← the atoms they are made of
        │
        ▼
jaffgen --network … --files …    ← expands $$JAFF templates against the model
  → generated/                   ← network-specific source out
```

Everything below follows that seam: first the model, then the generator.

---

## The `Network`

`Network` is the entry point for every JAFF workflow. One constructor call
reads a file, auto-detects its format, parses it into the typed catalogues, and
validates the chemistry.

```python
from jaff import Network

net = Network("networks/h_photoionization/h_photo.jet")

len(net.species)  # 3   — number of species
len(net.reactions)  # 2   — number of reactions
net.label  # 'h_photo'   — network identifier
```

A loaded `Network` carries:

- `species` — the [`Species`](../user-guide/working-with-networks/species.md)
  catalogue; look an entity up by name with `net.species["H"]`;
- `reactions` — the [`Reactions`](../user-guide/working-with-networks/reactions.md)
  catalogue;
- `elements` — the [`Elements`](../user-guide/working-with-networks/elements.md)
  catalogue, plus mass and composition information;
- `filename` — the path the network was read from;
- `label` — a human-readable identifier (defaults to the file stem).

<!-- prettier-ignore -->
!!! tip "Validate while you load"
    Chemistry checks run on every load — missing sinks and sources, isomer
    clashes, missing recombinations and mass- or charge-conservation violations are
    reported as warnings. Pass `strict=True` to stop at the first problem
    instead:
    ```python
    net = Network("mynetwork.dat", strict=True)
    ```

---

## Species

A **`Specie`** is one chemical entity — atom, molecule, or ion. From the single
name string in the file it derives the composition, mass, charge, and the
identifiers code generation needs.

```python
H = net.species["H"]

H.name  # 'H'
H.mass  # 1.673773e-24   ← grams, not amu
H.charge  # 0
H.index  # 0
H.latex()  # '{\\rm H}'      ← latex() is a method, not an attribute
```

The attributes worth knowing up front:

- `name` — the formula exactly as written in the file;
- `mass` — total mass in **grams (CGS)**, summed over the constituent atoms;
- `charge` — net charge in elementary-charge units (`0` neutral, `>0` cation,
  `<0` anion);
- `index` — zero-based position inside `net.species`, used everywhere arrays are
  indexed.

<!-- prettier-ignore -->
!!! warning "Mass is in grams (CGS), not atomic mass units"
    `specie.mass` is the physical mass in grams — `H` is `1.673773e-24`, not
    `1.008`. JAFF works in CGS so the value drops straight into rate and energy
    expressions. For the atomic weight in amu, read it from the element data
    instead.

The [Species page](../user-guide/working-with-networks/species.md) covers the
derived identifiers (`exploded`, `serialized`, `fidx`), the catalogue's bulk
accessors, and the electron's special handling.

---

## Reactions

A **`Reaction`** is one chemical transformation. It carries the species going
in, the species coming out, a symbolic rate, optional temperature bounds, and
an energy budget.

```python
rxn = net.reactions[0]

rxn.verbatim  # 'H + _PHOTON -> H+ + e-'   ← verbatim is an attribute
rxn.type  # 'photo'                     ← type is an attribute
rxn.rate_coefficient  # a SymPy expression for the rate coefficient
rxn.rate.expanded  # the coefficient times the reactant densities
rxn.reactants  # species consumed
rxn.products  # species created
```

- `reactants` / `products` — the species on each side, as `Species` catalogues;
- `rate_coefficient` — the rate coefficient as a **SymPy expression**, not a
  number; it still contains the temperature symbol so it can be differentiated
  and emitted as code;
- `rate` — `rate_coefficient` times the number density of each reactant (also a
  SymPy expression; pseudo-species such as `_PHOTON` carry no density);
- `type` — the classification (`photo`, `cosmic_ray`, `3_body`, `unknown`, or a
  grain surface-mechanism type such as `freeze`/`desorption_*` in gas-grain
  networks), concluded by the network-format parser as it reads the file;
- `tmin` / `tmax` — the temperature window the rate is valid over (`None` means
  unbounded);
- `verbatim` — the human-readable reaction string.

**Rate expressions.** Most thermal reactions use an Arrhenius-type law:

$$k(T) = \alpha \left(\frac{T}{300}\right)^\beta e^{-\gamma/T}$$

with $\alpha$ the pre-exponential factor, $\beta$ the temperature exponent,
$\gamma$ the activation parameter, and $T$ the temperature in Kelvin.
Photo-reactions instead carry a `photorates(...)` call and a `_PHOTON`
pseudo-reactant — the parser uses the latter to classify them as photo-reactions
(the `type` attribute no longer inspects the rate). The
[Reactions page](../user-guide/working-with-networks/reactions.md)
goes through every reaction type and the catalogue API.

---

## Network files

JAFF reads several community formats and detects which one it is looking at
automatically:

- **KIDA** — the KInetic Database for Astrochemistry;
- **UDFA** — the UMIST Database for Astrochemistry;
- **PRIZMO** — the PRIZMO astrochemical code;
- **KROME** — the KROME package for astrochemistry;
- **UCLCHEM** — the UCL Chemistry and Dust code.

Whatever the source, the loaded model is identical — that uniform
representation is the point of the parser.

---

## Code generation

JAFF's code generation is **template-driven**, and the
[`jaffgen`](../user-guide/code-generation/jaffgen.md) CLI is how you run it. You
write an ordinary source file — C, C++, Fortran, Python, Rust, Julia, or R — and
mark up the parts that depend on the network with `$$JAFF` directives. `jaffgen`
loads a network, expands every directive against it, and writes the result out.

A template is real source code with small generated regions. Every line that is
**not** a directive is copied verbatim; directive blocks are filled in:

```cpp
// rates.cpp — a template
$$JAFF SUB nreact
const int NREACT = $nreact$;
$$JAFF END

void compute_rates(double* k, double T) {
    $$JAFF REPEAT idx, rate IN rates
    k[$idx$] = $rate$;
    $$JAFF END
}
```

Run it against the network:

```bash
jaffgen --network networks/h_photoionization/h_photo.jet --files rates.cpp
```

The expanded file lands in `generated/`, keeping its name. `SUB` swaps a single
value; `REPEAT` loops a line once per item in a collection (here `rates`). The
directive lines themselves stay in the output as comments:

```cpp
// rates.cpp — a template
// $$JAFF SUB nreact
const int NREACT = 2;
// $$JAFF END

void compute_rates(double* k, double T) {
    // $$JAFF REPEAT idx, rate IN rates
    k[0] = photorates(1, 13.6, 1.0e+99);
    k[1] = …;
    // $$JAFF END
}
```

The directive language is small — `SUB`, `REPEAT`, `REDUCE`, `GET`, `HAS`,
`DEFINE`, `LANG`, `COMMENT` and `END`, plus options such as `REPLACE` and `USE`. The
[Template Syntax](../user-guide/code-generation/template-syntax.md) page is the
full reference. Three ideas explain what the collections emit: indexing, CSE,
and the ODE/Jacobian pair.

### Array indexing

Languages disagree on where arrays start and how they are bracketed. JAFF reads
the target language from each file's extension (or from a `$$JAFF LANG` line
in the file) and emits the right convention automatically:

| Language | Starting index | Example  |
| -------- | -------------- | -------- |
| C/C++    | 0              | `arr[0]` |
| Python   | 0              | `arr[0]` |
| Fortran  | 1              | `arr(1)` |

So the same template gives `k[0]` in C++ and `k(1)` in Fortran with no change.

### Common Subexpression Elimination (CSE)

The indexed collections (`rates`, `odes`, `jacobian`, …) can factor repeated
work into temporaries so each sub-expression is computed once. You opt in by
adding the `cse` variable to a `REPEAT`:

=== "Without CSE"

    ```cpp
    $$JAFF REPEAT idx, rate IN rates
    k[$idx$] = $rate$;
    $$JAFF END
    ```

=== "With CSE"

    ```cpp
    $$JAFF REPEAT idx, rate, cse IN rates
    const double x$idx$ = $cse$;
    k[$idx$] = $rate$;
    $$JAFF END
    ```

### ODEs and the Jacobian

A chemical network is a system of ODEs describing how concentrations change:

$$\frac{dy_i}{dt} = \sum_j \nu_{ij} R_j$$

with $y_i$ the concentration of species $i$, $R_j$ the rate of reaction $j$, and
$\nu_{ij}$ the stoichiometric coefficient of species $i$ in reaction $j$. The
`odes` collection emits exactly this right-hand side; implicit solvers also need
the **Jacobian** $J_{ij} = \partial f_i / \partial y_j$ (where $f_i = dy_i/dt$),
which the `jacobian` collection emits. Because every `rate` is a symbolic
expression, JAFF differentiates it exactly rather than numerically:

```cpp
$$JAFF REPEAT idx, ode IN odes
dydt[$idx$] = $ode$;
$$JAFF END

$$JAFF REPEAT idx, expr IN jacobian
J[$idx$][$idx$] = $expr$;
$$JAFF END
```

---

## Configuration files

Settings live in two TOML files with [one schema](../user-guide/configuration/index.md):
a `jaff.toml` next to a network file describes how that network is built
(temperature cutoffs, radiation, dust, the equation of state, …) and is picked
up by every `Network(...)` that loads it; a `jaffgen.toml` declares a
code-generation run once, so the command collapses to:

```bash
jaffgen --config jaffgen.toml
```

Command-line options (and `--set key.path=value`) still win over the files, so
a file is a baseline you can override per run. Each run writes its complete
config to `config.toml` in the output directory. Editors can
[complete and check](../user-guide/configuration/editor-setup.md) these files
as you type.

---

## Putting it together

A first analysis pass touches only the model — no code generation at all:

```python
from jaff import Network

net = Network("networks/h_photoionization/h_photo.jet", strict=True)

for sp in net.species:
    print(f"{sp.name}: {sp.mass:.3e} g  charge {sp.charge:+d}")

for rxn in net.reactions:
    print(f"{rxn.verbatim}  [{rxn.type}]")
```

The code-generation pass is a `jaffgen` invocation over your templates, exactly
as shown above. Everything else in the user guide is a deeper cut through one of
these objects.

---

## Common terms

| Term                 | Definition                                                      |
| -------------------- | --------------------------------------------------------------- |
| **Species**          | A chemical entity (atom, molecule, ion)                         |
| **Reaction**         | A chemical transformation between species                       |
| **Rate coefficient** | Function determining reaction speed                             |
| **Stoichiometry**    | Ratio of reactants to products                                  |
| **ODE**              | Ordinary Differential Equation describing concentration changes |
| **Jacobian**         | Matrix of partial derivatives of the ODEs                       |
| **CSE**              | Common Subexpression Elimination (an optimization)              |
| **Template**         | A source file with `$$JAFF` directives for code generation       |
| **Network**          | A collection of species and reactions                           |
| **Index offset**     | Starting index for generated arrays (`0` or `1`)                |

---

## Next steps

Now that the pieces have names:

1. [Working with Networks](../user-guide/working-with-networks/index.md) — inspect species, reactions, and elements in depth;
2. [Code Generation](../user-guide/code-generation/index.md) — run `jaffgen` over your templates;
3. [Template Syntax](../user-guide/code-generation/template-syntax.md) — write your own `$$JAFF` templates;
4. [API Reference](../reference/index.md) — the complete surface.
