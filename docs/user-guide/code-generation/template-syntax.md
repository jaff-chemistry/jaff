---
tags:
    - User-guide
    - Code-generation
    - Template
---

# Template Syntax

A JAFF template is an **ordinary source file** — C, C++, Fortran, Python, Rust,
Julia, R, or any custom format containing `$JAFF` directives. The engine
reads the file line by line and does exactly one of two things with each line:

- **A directive line** (or the content inside a directive block) is **expanded**
  against the loaded [`Network`](../working-with-networks/network.md);
- **Every other line** is copied to the output **verbatim**.

So a template is your real source code, with small generated regions created
by directives. The directive lines themselves are printed with the output too.

Templates are processed by the [`jaffgen`](jaffgen.md) CLI or the
[builder API](../advanced-code-generation/builder.md).

---

## Directive syntax

A directive is a **block** opened by a `$JAFF` line and closed by `$JAFF END`.
The line between them is the template; it is repeated/substituted per item:

```text
// $JAFF COMMAND arguments  $[modifiers]$
content with $token$ placeholders
// $JAFF END
```

It consists of three pieces:

- **Comment token** — taken automatically from the file extension: `//` for
  C/C++/Rust, `#` for Python/R, `!` for Fortran, plus `--` and `%`. A directive
  is only recognised when the line starts with the language's comment token
  immediately followed by `$JAFF`.
- **`$token$` placeholders** — substitution points in the body, written between
  dollar signs (`$nspec$`, `$specie$`, `$idx$`).
- **`$[modifiers]$`** — an optional trailing list (`SORT`, `REPLACE`,
  `RADIATION`, …); see [modifiers](#modifiers) and [REPLACE](#replace-directive).

---

## Commands

| Command  | Purpose                                                       |
| -------- | ------------------------------------------------------------- |
| `SUB`    | Substitute scalar network values (counts, label, …)           |
| `REPEAT` | Iterate a collection — one block per item, or an inline array |
| `GET`    | Retrieve one property of a named entity                       |
| `HAS`    | Test whether a species, reaction, or element exists           |
| `REDUCE` | Build a summed expression over a collection                   |
| `END`    | Close the current directive block                             |

Every command accepts a trailing `$[... REPLACE pattern replacement ...]$`
modifier (see [REPLACE Directive](#replace-directive)).

All examples below use the hydrogen photo-ionization network (`H`, `H+`, `e-`).

---

## SUB

Substitute one or more scalar network values.

```text
// $JAFF SUB var1, var2, ...
content with $var1$, $var2$ tokens
// $JAFF END
```

**Available tokens:**

| Token      | Type   | Description                                |
| ---------- | ------ | ------------------------------------------ |
| `nspec`    | `int`  | Number of species                          |
| `nreact`   | `int`  | Number of reactions                        |
| `nelem`    | `int`  | Number of elements                         |
| `nbands`   | `int`  | Number of radiation bands (`0` if none)    |
| `label`    | `str`  | Network label                              |
| `filename` | `str`  | Template file name                         |
| `filepath` | `Path` | Full template file path                    |
| `e_idx`    | `int`  | Index of the electron species              |
| `dedt`     | `str`  | Language-specific total `dE/dt` expression (see `DEDT_TYPE`) |
| `dtdt`     | `str`  | Language-specific total `dT/dt` expression |

Integer tokens accept inline arithmetic with `+`, `-`, `*`, `/`:

```cpp
// $JAFF SUB nspec
#define NSPEC    $nspec$
#define NSPEC_P1 $nspec+1$
// $JAFF END
```

**Output:**

```cpp
// $JAFF SUB nspec
#define NSPEC    3
#define NSPEC_P1 4
// $JAFF END
```

Several tokens can share one directive:

```cpp
// $JAFF SUB nspec, nreact
#define NSPEC    $nspec$
#define NREACT   $nreact$
// $JAFF END
```

---

## REPEAT

Iterate a collection. With `idx` in the variable list, the body is emitted once
per item (**vertical**); without it, the body's bracketed list is filled inline
(**horizontal**).

### Vertical mode

One output line per item. `$idx$` is the zero-based index; `$idx+N$` / `$idx-N$`
offset it. The number of `$idx$` tokens on a line must match the collection's
dimensionality (1 for lists, 2 for matrices).

```cpp
// $JAFF REPEAT idx, specie IN species
species_names[$idx$] = "$specie$";
// $JAFF END
```

Output:

```cpp
// $JAFF REPEAT idx, specie IN species
species_names[0] = "H";
species_names[1] = "H+";
species_names[2] = "e-";
// $JAFF END
```

For a 2D collection, use two `$idx$` tokens:

```cpp
// $JAFF REPEAT idx, element IN element_density_matrix
dens[$idx$][$idx$] = $element$;
// $JAFF END
```

### Horizontal mode

With no `idx`, the body must contain the variable inside brackets. The engine
fills the bracket with every item, reusing the surrounding bracket, optional
quotes, and the separator it finds:

```cpp
// $JAFF REPEAT specie_charge IN specie_charges
double charges[] = {$specie_charge$, };
// $JAFF END
```

Output:

```cpp
// $JAFF REPEAT specie_charge IN specie_charges
double charges[] = {0, 1, -1};
// $JAFF END
```

Everything between the bracket and the variable is treated as the separator (it
recognises `,`, `;`, `:` and whitespace), and quotes around the token are
preserved:

```cpp
// $JAFF REPEAT specie IN species
std::vector<std::string> species = {"$specie$", };
// $JAFF END
```

Output:

```cpp
std::vector<std::string> species = {"H", "H+", "e-"};
```

### CSE (common subexpression elimination)

`rates`, `odes`, `rhses`, `radodes`, and `jacobian` support a third variable, `cse`, that
holds repeated sub-expressions into named temporaries. The temporaries are
emitted only where a sub-expression is **actually shared** across items — a tiny
network with no repetition produces none.

```cpp
// $JAFF REPEAT idx, rate, cse IN rates
const double x$idx$ = $cse$;
k[$idx$] = $rate$;
// $JAFF END
```

On a network with shared sub-expressions (here the GOW network), the `$cse$`
line expands to the temporary cascade and the rate line references them:

```cpp
const double x0 = nden[0] + nden[1] + 2*nden[3] + 2*nden[4] + ...;
const double x1 = crate*(1.5*nden[0] + 2.3*nden[3])/x0;
const double x2 = 1.0/nden[0];
// ... x3 … xN ...
k[0] = x1;
k[1] = ...;
```

The temporaries are named after the identifier that contains `$idx$` on the
`$cse$` line: `x$idx$` gives `x0, x1, …` and `tmp$idx$_value` gives
`tmp0_value, tmp1_value, …`, with or without spaces around `=`. That line must
be the first line of the block, since the temporaries have to be declared before
the other lines use them. The identifier must start with a letter or underscore
before `$idx$`, and the index takes no offset, so `$idx$x`, `cse[$idx$]` and
`cse$idx+1$` are rejected with a `ParserError`.

### Modifiers

Modifiers go inside `$[...]$` at the end of the command line.

| Modifier    | Values       | Description                                     | Supported           |
| ----------- | ------------ | ----------------------------------------------- | ------------------- |
| `SORT`      | `TRUE/FALSE` | Sort items before expansion                     | All                 |
| `THERMAL`   | `none/dedt/dtdt` | Thermal ODE row: none, `dE/dt` or `dT/dt` (default `dedt` for `rhses`, `none` for `jacobian`) | `rhses`, `jacobian` |
| `DEDT_TYPE` | `volumetric/specific/per_particle/molar` | Internal-energy form used when `THERMAL dedt` (default `volumetric`) | `rhses`, `jacobian` |
| `RADIATION` | `TRUE/FALSE` | Include radiation ODE / Jacobian terms          | `rhses`, `jacobian` |
| `REPLACE`   | `pat repl`   | Regex replacement on the output                 | All                 |
| `POS`       | `str`        | `+` replacement string for normalized sign      | `species_with_normalized_sign` |
| `NEG`       | `str`        | `-` replacement string for normalized sign      | `species_with_normalized_sign` |

`REPLACE` rewrites generated output after expansion. This is handy for mapping JAFF's
standard symbols (`tgas`, `nden[…]`, `radeden[…]`, `photden[…]`, `rflux[…]`) onto your
code's own names:

```cpp
// $JAFF REPEAT idx, rate IN rates $[REPLACE tgas T]$
state.rates[$idx$] = $rate$;
// $JAFF END
```

Modifiers chain. A realistic Jacobian directive remaps several array accessors
and switches on radiation and the thermal (`dE/dt`) row at once:

```cpp
// $JAFF REPEAT idx, expr, cse IN jacobian $[REPLACE nden\[\s*(\d+)\s*\] state.xn[\1] REPLACE radeden\[\s*(\d+)\s*\] state.rn[2*\1] REPLACE photden\[\s*(\d+)\s*\] state.rn[2*\1] REPLACE rflux\[\s*(\d+)\s*\] state.rn[2*\1+1] RADIATION True THERMAL dedt]$
```

`POS` and `NEG` set the strings substituted for the `+` and `-` charge markers when
expanding `species_with_normalized_sign`. They default to `j` / `k` (so `H+` → `hj`,
`e-` → `ek`); override them to match another convention, e.g. `$[POS j NEG k]$`.

### Expression-generating collections

These produce indexed code expressions. Only the five marked support `cse`.

| Collection         | Variables              | Description                      | CSE |
| ------------------ | ---------------------- | -------------------------------- | --- |
| `rates`            | `idx, rate`            | Rate-coefficient expressions     | ✓   |
| `odes`             | `idx, ode`             | Full ODE expressions (dy/dt)     | ✓   |
| `rhses`            | `idx, rhs`             | RHS including dE/dt              | ✓   |
| `jacobian`         | `idx, expr`            | Jacobian matrix elements (2D)    | ✓   |
| `radodes`          | `idx, radode`          | Radiation moment ODEs            | ✓   |
| `flux_expressions` | `idx, flux_expression` | Flux = rate × reactant densities | ✗   |
| `ode_expressions`  | `idx, ode_expression`  | ODE terms without assignment     | ✗   |

### List-iterating collections

These iterate plain value lists.

**Species:**

| Collection                     | Variable                      |
| ------------------------------ | ----------------------------- |
| `species`                      | `specie`                      |
| `species_with_normalized_sign` | `specie_with_normalized_sign` |
| `specie_masses`                | `specie_mass`                 |
| `specie_charges`               | `specie_charge`               |
| `specie_masses_ne`             | `specie_mass_ne`              |
| `specie_charges_ne`            | `specie_charge_ne`            |
| `neutral_species`              | `neutral_specie`              |
| `charged_species`              | `charged_specie`              |
| `neutral_specie_indices`       | `neutral_specie_index`        |
| `charged_specie_indices`       | `charged_specie_index`        |
| `charged_specie_indices_ne`    | `charged_specie_index_ne`     |
| `neutral_specie_masses`        | `neutral_specie_mass`         |
| `charged_specie_masses`        | `charged_specie_mass`         |
| `charged_specie_masses_ne`     | `charged_specie_mass_ne`      |
| `charge_truths`                | `charge_truth`                |
| `charge_truths_ne`             | `charge_truth_ne`             |

**Reactions:**

| Collection               | Variable               |
| ------------------------ | ---------------------- |
| `reactions`              | `reaction`             |
| `reactants`              | `reactant`             |
| `products`               | `product`              |
| `photo_reactions`        | `photo_reaction`       |
| `photo_reaction_indices` | `photo_reaction_index` |
| `photo_reaction_truths`  | `photo_reaction_truth` |
| `tmins`                  | `tmin`                 |
| `tmaxes`                 | `tmax`                 |

**Elements:**

| Collection               | Variable  |
| ------------------------ | --------- |
| `elements`               | `element` |
| `element_density_matrix` | `element` |
| `element_truth_matrix`   | `element` |

**Radiation:**

| Collection             | Variable             |
| ---------------------- | -------------------- |
| `radiation_boundaries` | `radiation_boundary` |

Band edges, of length `nbands + 1`. Empty if the network has no radiation.

---

## GET

Retrieve one property of a named entity (species, reaction, or element).

```text
// $JAFF GET property FOR entity_name
content with $property$ token
// $JAFF END
```

**Available properties:**

| Property            | Entity            | Returns                             |
| ------------------- | ----------------- | ----------------------------------- |
| `specie_idx`        | species name      | `int` — index in the species array  |
| `specie_mass`       | species name      | `float` — mass in **grams (CGS)**   |
| `specie_charge`     | species name      | `int` — charge                      |
| `specie_latex`      | species name      | `str` — LaTeX name                  |
| `element_idx`       | element symbol    | `int` — index in the element array  |
| `reaction_idx`      | verbatim reaction | `int` — dense catalogue position in the reaction array (`Reaction.catalogue_index`, not the gapped source `index`) |
| `reaction_tmin`     | verbatim reaction | `float` — minimum temperature       |
| `reaction_tmax`     | verbatim reaction | `float` — maximum temperature       |
| `reaction_verbatim` | verbatim reaction | `str` — verbatim string             |

```cpp
// $JAFF GET specie_idx FOR H+
const int hplus_idx = $specie_idx$;
// $JAFF END

// $JAFF GET specie_mass FOR H
const double h_mass = $specie_mass$;   // grams (CGS)
// $JAFF END
```

Output:

```cpp
const int hplus_idx = 1;
const double h_mass = 1.673773e-24;   // grams (CGS)
```

<!-- prettier-ignore -->
!!! note "`specie_mass` is in grams, not amu"
    Like [`Specie.mass`](../working-with-networks/species.md), this value is the
    physical mass in CGS grams. For the atomic weight in amu, read it from the
    element data instead.

---

## HAS

Return `1` if an entity exists in the network, `0` otherwise.

```text
// $JAFF HAS identity entity_name
content with $identity$ token (expands to 0 or 1)
// $JAFF END
```

**Identities:** `specie`, `reaction`, `element`.

```cpp
// $JAFF HAS specie e-
const int HAS_ELECTRON = $specie$;
// $JAFF END

// $JAFF HAS element C
const int HAS_CARBON = $element$;
// $JAFF END
```

Output (the hydrogen network has the electron but no carbon):

```cpp
// $JAFF HAS specie e-
const int HAS_ELECTRON = 1;
// $JAFF END

// $JAFF HAS element C
const int HAS_CARBON = 0;
// $JAFF END
```

---

## REDUCE

Build a single summed expression over a collection by expanding the `$( ... )$`
region once per item and joining the pieces with `+`. Each piece and the whole
sum are wrapped in parentheses, so the result is safe to use as an operand
(`2 * $(...)$`, `1 / $(...)$`). A reduction over an empty collection expands
to `0`.

```text
// $JAFF REDUCE var IN collection
result = $( ... $var$ ... )$;
// $JAFF END
```

```cpp
// $JAFF REDUCE specie_charge IN specie_charges
double total_charge = $($specie_charge$)$;
// $JAFF END
```

Output:

```cpp
// $JAFF REDUCE specie_charge IN specie_charges
double total_charge = ((0) + (1) + (-1));
// $JAFF END
```

Multiple collections can be reduced together (they must be the same length) —
each `$var$` inside the `$()$` is indexed in lockstep. For example, electron
density from charge balance over the charged, non-electron species:

```cpp
// $JAFF REDUCE charged_specie_index_ne, charged_specie_charge_ne IN charged_specie_indices_ne, charged_specie_charges_ne
ne = $(($charged_specie_charge_ne$ * xn[$charged_specie_index_ne$]))$;
// $JAFF END
```

For the hydrogen network the only charged non-electron species is `H+`
(index 1, charge +1):

```cpp
ne = (((1 * xn[1])));
```

---

## REPLACE Directive

Every command accepts a `REPLACE` modifier inside `$[...]$` that applies a Python
regex substitution to the generated output **after** expansion. It is the
template engine's most powerful feature — it lets you target languages or
naming conventions JAFF does not model natively. Multiple `REPLACE` pairs chain
and apply in sequence.

```text
// $JAFF REPEAT idx, specie IN species $[REPLACE \+ _plus REPLACE - _minus]$
const int idx_$specie$ = $idx$;
// $JAFF END
```

Without `REPLACE`, the raw species names leak invalid identifier characters:

```cpp
const int idx_H+ = 1;
const int idx_e- = 2;
```

With the two `REPLACE` pairs above:

```cpp
const int idx_H = 0;
const int idx_H_plus = 1;
const int idx_e_minus = 2;
```

Capture groups work too — to subscript digits in species names:

```text
// $JAFF REPEAT idx, specie IN species $[REPLACE (H)(\d+) \1_\2]$
arr[$idx$] = "$specie$";
// $JAFF END
```

`"H2"` → `"H_2"`, `"H2O"` → `"H_2O"`, and so on.
