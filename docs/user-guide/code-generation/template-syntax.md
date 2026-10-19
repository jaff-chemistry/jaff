---
tags:
    - User-guide
    - Code-generation
    - Template
---

# Template Syntax

A JAFF template is an **ordinary source file** — C, C++, Fortran, Python, Rust,
Julia, R, or any custom format — containing `$$JAFF` directives. The engine
reads the file line by line and does exactly one of two things with each line:

- **A directive line** (or the content inside a directive block) is **expanded**
  against the loaded [`Network`](../working-with-networks/network.md);
- **Every other line** is copied to the output **verbatim**.

So a template is your real source code, with small generated regions created
by directives. The directive lines themselves are written to the output as
comments (see [Echo](#echo)).

Templates are processed by the [`jaffgen`](jaffgen.md) CLI.

---

## Directive syntax

A directive starts on a line whose first non-blank text is `$$JAFF`. Directive
lines carry **no comment prefix** in the template. A directive with a body is
closed by `$$JAFF END`; the lines between the two are the template, repeated or
substituted per item:

```text
$$JAFF COMMAND arguments { options }
content with $token$ placeholders
$$JAFF END
```

It consists of these pieces:

- **Header** — `$$JAFF`, the command name and its arguments. Arguments are
  separated by commas or spaces.
- **Option block** — an optional `{ … }` after the arguments. The opening `{`
  must be on the header line; the block itself may continue over several
  lines. `#` starts a comment that runs to the end of the line (outside quoted
  strings), in headers and option blocks alike. Line breaks are allowed between
  options, but each `"pattern" -> "replacement"` pair, each `KEY -> value`
  setting and each `USE name` must sit on one line; a break inside one is a
  `ParserError`.
- **Body** — the template lines up to `$$JAFF END`. `DEFINE`, `LANG` and
  `COMMENT` have no body and no `END`.
- **`$token$` placeholders** — substitution points in the body, written between
  dollar signs (`$nspec$`, `$specie$`, `$idx$`).

A multi-line header looks like this:

```cpp
$$JAFF REPEAT idx, rate IN rates {
    REPLACE { "tgas" -> "T" }      # rename the temperature symbol
    REPLACE { "nden" -> "n" }
}
state.rates[$idx$] = $rate$;
$$JAFF END
```

Comments are handy for switching options off without deleting them:

```cpp
$$JAFF REPEAT idx, rhs IN rhses {
    REPLACE {
        "nden\[(\d+)\]" -> "y[\1]"
        # "tgas" -> "T"               # this pair is disabled
    }
    # THERMAL { MODE -> dtdt }        # so is this whole group
}
f[$idx$] = $rhs$;
$$JAFF END
```

Header comments only exist inside the header and its option block. In the body a
`#` is ordinary template text and is copied to the output. Commented-out header
lines are still part of the directive, so with [echo](#echo) on they appear in
the echoed header.

Directive blocks **do not nest**: a `$$JAFF` line inside a body is an error.
Indentation of the directive and of the body lines (spaces and tabs) is
preserved in the output, so directives can sit inside a function body. Output
lines always end in `\n`, also for templates written with `\r\n`.

---

## Commands

| Command   | Purpose                                                       |
| --------- | ------------------------------------------------------------- |
| `SUB`     | Substitute scalar network values (counts, label, …)           |
| `REPEAT`  | Iterate a collection — one line per item, or an inline array  |
| `GET`     | Retrieve properties of a named entity                         |
| `HAS`     | Test whether a species, reaction, or element exists           |
| `REDUCE`  | Build a summed expression over collections                   |
| `END`     | Close the current directive block                             |
| `DEFINE`  | Name a set of options for reuse with `USE` (no body)          |
| `USE`     | Pull a `DEFINE`d option set into an option block              |
| `LANG`    | Set the target language of the file (no body)                 |
| `COMMENT` | Set the comment token used to echo directives (no body)       |

Every command with a body accepts the global `REPLACE` and `USE` options (see
[Options](#options)).

Each example below was rendered against a small test network with the species
`H`, `e-`, `H+`, `He++`, `He+` (in that order) unless stated otherwise.

---

## SUB

Substitute one or more scalar network values.

```text
$$JAFF SUB token1, token2, ...
content with $token1$, $token2$ tokens
$$JAFF END
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
| `dedt`     | `str`  | Language-specific total `dE/dt` expression |
| `dtdt`     | `str`  | Language-specific total `dT/dt` expression |

Integer tokens accept inline arithmetic with `+`, `-`, `*`, `/` (see
[Placeholders](#placeholders)):

```cpp
$$JAFF SUB nspec
#define NSPEC    $nspec$
#define NSPEC_P1 $nspec+1$
$$JAFF END
```

Output:

```cpp
// $$JAFF SUB nspec
#define NSPEC    5
#define NSPEC_P1 6
// $$JAFF END
```

Several tokens can share one directive:

```cpp
$$JAFF SUB nspec, nreact
#define NSPEC    $nspec$
#define NREACT   $nreact$
$$JAFF END
```

---

## REPEAT

```text
$$JAFF REPEAT vars IN property { options }
```

Iterate a collection. With `idx` in the variable list, the body is emitted once
per item (**vertical**); without it, the body's bracketed list is filled inline
(**horizontal**).

### Vertical mode

Every body line holding a placeholder of the REPEAT (`$idx…$`, `$IDX$`, the item
variable or a declared `$cse$`) is written once per item; other lines are copied
once. `$idx$` is the zero-based index; `$idx+N$`, `$idx-N$`, `$idx*N$` and
`$idx/N$` transform it. If a line has `$idx$` tokens, their number must match the
collection's dimensionality (1 for lists, 2 for matrices).

```cpp
$$JAFF REPEAT idx, specie IN species
species_names[$idx$] = "$specie$";
$$JAFF END
```

Output:

```cpp
// $$JAFF REPEAT idx, specie IN species
species_names[0] = "H";
species_names[1] = "e-";
species_names[2] = "H+";
species_names[3] = "He++";
species_names[4] = "He+";
// $$JAFF END
```

For a 2D collection, use two `$idx$` tokens (row, then column):

```cpp
$$JAFF REPEAT idx, element IN element_density_matrix
dens[$idx$][$idx$] = $element$;
$$JAFF END
```

Output:

```cpp
// $$JAFF REPEAT idx, element IN element_density_matrix
dens[0][0] = 1;
dens[0][1] = 0;
dens[0][2] = 1;
dens[0][3] = 0;
dens[0][4] = 0;
dens[1][0] = 0;
dens[1][1] = 0;
dens[1][2] = 0;
dens[1][3] = 1;
dens[1][4] = 1;
// $$JAFF END
```

In Fortran, `$idx$` is still zero-based; write `$idx+1$` where the target array
is one-based.

### Horizontal mode

With no `idx`, the body must contain the variable inside brackets. The engine
fills the bracket with every item, reusing the surrounding bracket, optional
quotes, and the separator it finds:

```cpp
$$JAFF REPEAT specie_charge IN specie_charges
double charges[] = {$specie_charge$, };
$$JAFF END
```

Output:

```cpp
// $$JAFF REPEAT specie_charge IN specie_charges
double charges[] = {0, -1, 1, 2, 1};
// $$JAFF END
```

The recognised brackets are `(`, `{`, `<`, `[` and Fortran's `(/` with the
matching closers `)`, `}`, `>`, `]` and `/)`. A line without a bracketed
variable is copied (with `REPLACE` applied). Between the opening bracket and the
variable only whitespace and an optional quote (`'` or `"`, matching the one after the variable) are
allowed, and the quote is preserved. The separator is what follows the variable
before the closing bracket (`{$x$, }` joins with `", "`; it may consist of `,`,
`;`, `:` and whitespace). When none is written, `", "` is used:

```cpp
$$JAFF REPEAT specie IN species
std::vector<std::string> species = {"$specie$", };
$$JAFF END
```

Output:

```cpp
// $$JAFF REPEAT specie IN species
std::vector<std::string> species = {"H", "e-", "H+", "He++", "He+"};
// $$JAFF END
```

### CSE (common subexpression elimination)

`rates`, `odes`, `rhses`, `radodes`, and `jacobian` support a third variable, `cse`, that
holds repeated sub-expressions as named temporaries. The temporaries are
emitted only where a sub-expression is **actually shared** across items — a tiny
network with no repetition produces none.

```cpp
$$JAFF REPEAT idx, ode, cse IN odes
const double x$idx$ = $cse$;
dydt[$idx$] = $ode$;
$$JAFF END
```

Output (here two sub-expressions are shared, so the `$cse$` line expands to two
temporaries and the ODE lines reference them; without shared sub-expressions
the `$cse$` line produces nothing):

```cpp
// $$JAFF REPEAT idx, ode, cse IN odes
const double x0 = 4.1207223003988966e-10*std::pow(tgas, -0.80000000000000004)*nden[1]*nden[2] - 5.8500000000000005e-11*std::sqrt(tgas)*std::exp(-157809.10000000001/tgas)*nden[0]*nden[1];
const double x1 = nden[1]*std::pow(nden[3], 2);
dydt[0] = x0;
dydt[1] = -x0 - x1;
dydt[2] = -x0;
dydt[3] = -x1;
dydt[4] = x1;
// $$JAFF END
```

The temporaries are named after the identifier that contains `$idx$` on the
`$cse$` line: `x$idx$` gives `x0, x1, …` and `tmp$idx$_value` gives
`tmp0_value, tmp1_value, …`, with or without spaces around `=`. That line must
be the first line of the block that expands, since the temporaries have to be
declared before the other lines use them. The identifier must start with a
letter or underscore before `$idx$`, and the index takes no arithmetic, so
`$idx$x`, `cse[$idx$]` and `cse$idx+1$` are rejected with a `ParserError`. A
`$cse$` line is only special when `cse` is in the variable list.

### Properties

`REPEAT vars IN property` — the variables are written in the order shown:
`idx`, the item variable, then `cse` where the property supports it.

#### Expression-generating properties

These produce code expressions (rendered in the target language). Only those
with a `cse` variable support CSE.

| Property           | Variables              | Description                      | Option groups                               |
| ------------------ | ---------------------- | -------------------------------- | ------------------------------------------- |
| `rates`            | `idx, rate, cse`       | Rate-coefficient expressions     | —                                           |
| `odes`             | `idx, ode, cse`        | Full ODE expressions (dy/dt)     | —                                           |
| `rhses`            | `idx, rhs, cse`        | RHS including the thermal row    | `THERMAL`, `RADIATION`              |
| `jacobian`         | `idx, expr, cse`       | Jacobian matrix elements (2D)    | `THERMAL`, `RADIATION`              |
| `radodes`          | `idx, radode, cse`     | Radiation moment ODEs            | `RADIATION` (`ORDER` only)          |
| `flux_expressions` | `idx, flux_expression` | Flux = rate × reactant densities | —                                           |
| `ode_expressions`  | `idx, ode_expression`  | ODE terms without assignment     | —                                           |

#### Species

| Property                       | Variables                                  | Option groups  |
| ------------------------------ | ------------------------------------------ | -------------- |
| `species`                      | `idx, specie`                              | —              |
| `species_with_normalized_sign` | `idx, specie_with_normalized_sign`         | `SIGN` |
| `specie_masses`                | `idx, specie_mass`                         | —              |
| `specie_charges`               | `idx, specie_charge`                       | —              |
| `specie_masses_ne`             | `idx, specie_mass_ne`                      | —              |
| `specie_charges_ne`            | `idx, specie_charge_ne`                    | —              |
| `neutral_species`              | `idx, neutral_specie`                      | —              |
| `charged_species`              | `idx, charged_specie`                      | —              |
| `neutral_specie_indices`       | `idx, neutral_specie_index`                | —              |
| `neutral_indices`              | `idx, neutral_index`                       | —              |
| `charged_specie_indices`       | `idx, charged_specie_index`                | —              |
| `charged_indices`              | `idx, charge_index`                        | —              |
| `charged_specie_indices_ne`    | `idx, charged_specie_index_ne`             | —              |
| `neutral_specie_masses`        | `idx, neutral_specie_mass`                 | —              |
| `charged_specie_masses`        | `idx, charged_specie_mass`                 | —              |
| `charged_specie_masses_ne`     | `idx, charged_specie_mass_ne`              | —              |
| `charge_truths`                | `idx, charge_truth`                        | —              |
| `charge_truths_ne`             | `idx, charge_truth_ne`                     | —              |

The `_ne` variants exclude the electron. `neutral_indices` and `charged_indices`
list the same values as `neutral_specie_indices` and `charged_specie_indices`
under shorter variable names.

#### Reactions

| Property                 | Variables                      | Description                                  | Option groups |
| ------------------------ | ------------------------------ | -------------------------------------------- | ------------- |
| `reactions`              | `idx, reaction`                | Verbatim reaction strings                    | —             |
| `reactants`              | `idx, reactant`                | Reactant names, one row per reaction (2D)    | —             |
| `products`               | `idx, product`                 | Product names, one row per reaction (2D)     | —             |
| `photo_reactions`        | `idx, photo_reaction`          | Photo-reaction strings                       | —             |
| `photo_reaction_indices` | `idx, photo_reaction_index`    | Indices of the photo reactions               | —             |
| `photo_reaction_truths`  | `idx, photo_reaction_truth`    | Per reaction: is it a photo reaction         | —             |
| `tmins`                  | `idx, tmin`                    | Minimum temperature of each reaction         | —             |
| `tmaxes`                 | `idx, tmax`                    | Maximum temperature of each reaction         | —             |

#### Elements and radiation

| Property                 | Variables                      | Description                                        | Option groups |
| ------------------------ | ------------------------------ | -------------------------------------------------- | ------------- |
| `elements`               | `idx, element`                 | Element symbols                                    | —             |
| `element_density_matrix` | `idx, element`                 | Element content of each species (2D)               | —             |
| `element_truth_matrix`   | `idx, element`                 | Whether each species contains each element (2D)    | —             |
| `radiation_boundaries`   | `idx, radiation_boundary`      | Band edges, length `nbands + 1`; empty without radiation | —       |

### REPEAT option groups

`REPEAT` takes option groups per property (last column above). Anything not
listed for a property is an error. `REPLACE` and `USE` are accepted
everywhere; see [Options](#options) for their general rules.

| Group       | Key       | Values                                       | Description                                                                                       |
| ----------- | --------- | -------------------------------------------- | ------------------------------------------------------------------------------------------------- |
| `THERMAL`   | `MODE`    | `none`, `dedt`, `dtdt`                       | Thermal ODE row: none, `dE/dt` or `dT/dt`. Default `dedt` for `rhses`, `none` for `jacobian`.     |
| `THERMAL`   | `ENERGY`  | `volumetric`, `specific`, `per_particle`, `molar` | Internal-energy form used when `MODE` is `dedt` (default `volumetric`).                      |
| `RADIATION` | `ENABLED` | `True`/`False`                               | Include radiation ODE / Jacobian terms (`rhses`, `jacobian`).                                     |
| `RADIATION` | `ORDER`   | integer                                      | Order of the radiation moments (`rhses`, `jacobian`, `radodes`).                                  |
| `SIGN`      | `POS`     | text                                         | String for the `+` charge marker (`species_with_normalized_sign`, default `j`).                   |
| `SIGN`      | `NEG`     | text                                         | String for the `-` charge marker (`species_with_normalized_sign`, default `k`).                   |

Booleans accept `True` and `False` in any letter case.

A realistic Jacobian directive remaps array accessors and picks the thermal
(`dE/dt`) row and energy form:

```cpp
$$JAFF REPEAT idx, expr, cse IN jacobian {
    REPLACE {
        "nden\[\s*(\d+)\s*\]" -> "state.xn[\1]"
        "tgas" -> "state.T"
    }
    THERMAL { MODE -> dedt  ENERGY -> specific }
    RADIATION { ENABLED -> False }
}
cse$idx$ = $cse$;
J[$idx$][$idx$] = $expr$;
$$JAFF END
```

Output (abridged — the full block has 30 temporaries; the extra column `[5]` is
the thermal row):

```cpp
// $$JAFF REPEAT idx, expr, cse IN jacobian {
    // REPLACE {
        // "nden\[\s*(\d+)\s*\]" -> "state.xn[\1]"
        // "tgas" -> "state.T"
    // }
    // THERMAL { MODE -> dedt  ENERGY -> specific }
    // RADIATION { ENABLED -> False }
// }
cse0 = std::exp(-157809.10000000001/state.T);
cse1 = cse0*state.xn[1];
cse2 = std::sqrt(state.T);
...
J[0][0] = -cse13*cse9 - cse4;
J[0][1] = -cse13*cse17 + cse15 - cse16;
...
J[4][3] = cse27;
// $$JAFF END
```

`POS` and `NEG` set the strings substituted for the `+` and `-` charge markers when
expanding `species_with_normalized_sign`. They default to `j` / `k` (so `H+` → `hj`,
`e-` → `ek`):

```cpp
$$JAFF REPEAT idx, specie_with_normalized_sign IN species_with_normalized_sign
default_$idx$ = "$specie_with_normalized_sign$";
$$JAFF END

$$JAFF REPEAT idx, specie_with_normalized_sign IN species_with_normalized_sign { SIGN { POS -> _plus  NEG -> _minus } }
custom_$idx$ = "$specie_with_normalized_sign$";
$$JAFF END
```

Output:

```cpp
// $$JAFF REPEAT idx, specie_with_normalized_sign IN species_with_normalized_sign
default_0 = "h";
default_1 = "ek";
default_2 = "hj";
default_3 = "hejj";
default_4 = "hej";
// $$JAFF END

// $$JAFF REPEAT idx, specie_with_normalized_sign IN species_with_normalized_sign { SIGN { POS -> _plus  NEG -> _minus } }
custom_0 = "h";
custom_1 = "e_minus";
custom_2 = "h_plus";
custom_3 = "he_plus_plus";
custom_4 = "he_plus";
// $$JAFF END
```

---

## GET

```text
$$JAFF GET property1, property2 IN entity { options }
content with $property1$, $property2$ tokens
$$JAFF END
```

Retrieve properties of a named entity (species, reaction, element, or
auxiliary function). Each property is substituted as `$property$` in the body.
The entity must exist in the network and be of the kind the property needs;
this is checked before anything renders. A reaction is named by its verbatim
text, which contains spaces, so it must be a quoted string, e.g.
`GET reaction_idx IN "H+ + e- -> H"`. A reaction name shared by several
mechanisms is rejected by `GET`.

**Available properties:**

| Property            | Entity             | Returns                             |
| ------------------- | ------------------ | ----------------------------------- |
| `specie_idx`        | species name       | `int` — index in the species array  |
| `specie_mass`       | species name       | `float` — mass in **grams (CGS)**   |
| `specie_charge`     | species name       | `int` — charge                      |
| `specie_latex`      | species name       | `str` — LaTeX name                  |
| `element_idx`       | element symbol     | `int` — index in the element array  |
| `reaction_idx`      | verbatim reaction  | `int` — dense catalogue position in the reaction array (`Reaction.catalogue_index`, not the gapped source `index`) |
| `reaction_tmin`     | verbatim reaction  | `float` — minimum temperature       |
| `reaction_tmax`     | verbatim reaction  | `float` — maximum temperature       |
| `reaction_verbatim` | verbatim reaction  | `str` — verbatim string             |
| `aux_func`          | auxiliary function | `str` — the function definition as target-language code |

```cpp
$$JAFF GET specie_idx IN H+
const int hplus_idx = $specie_idx$;
$$JAFF END

$$JAFF GET specie_mass IN H
const double h_mass = $specie_mass$;   // grams (CGS)
$$JAFF END
```

Output:

```cpp
// $$JAFF GET specie_idx IN H+
const int hplus_idx = 2;
// $$JAFF END

// $$JAFF GET specie_mass IN H
const double h_mass = 1.673773e-24;   // grams (CGS)
// $$JAFF END
```

A reaction is given as a quoted string:

```cpp
$$JAFF GET reaction_idx, reaction_tmin IN "H+ + e- -> H"
idx = $reaction_idx$; tmin = $reaction_tmin$;
$$JAFF END

$$JAFF HAS reaction "H+ + e- -> H"
has = $reaction$;
$$JAFF END
```

Output (this reaction has no temperature range, so `tmin` is `None`):

```cpp
// $$JAFF GET reaction_idx, reaction_tmin IN "H+ + e- -> H"
idx = 1; tmin = None;
// $$JAFF END

// $$JAFF HAS reaction "H+ + e- -> H"
has = 1;
// $$JAFF END
```

<!-- prettier-ignore -->
!!! note "`specie_mass` is in grams, not amu"
    Like [`Specie.mass`](../working-with-networks/species.md), this value is the
    physical mass in CGS grams. For the atomic weight in amu, read it from the
    element data instead.

---

## HAS

```text
$$JAFF HAS kind entity { options }
content with $kind$ token (expands to 0 or 1)
$$JAFF END
```

Substitute `1` if an entity exists in the network, `0` otherwise. The
placeholder is the kind itself. A missing entity is not an error (that is the
point): it gives `0`. Reactions are named by quoted verbatim text, as for `GET`.

**Kinds:** `specie`, `reaction`, `element`.

```cpp
$$JAFF HAS specie e-
const int HAS_ELECTRON = $specie$;
$$JAFF END

$$JAFF HAS element C
const int HAS_CARBON = $element$;
$$JAFF END
```

Output (the test network has the electron but no carbon):

```cpp
// $$JAFF HAS specie e-
const int HAS_ELECTRON = 1;
// $$JAFF END

// $$JAFF HAS element C
const int HAS_CARBON = 0;
// $$JAFF END
```

---

## REDUCE

```text
$$JAFF REDUCE vars IN properties { options }
result = $( ... $var$ ... )$;
$$JAFF END
```

Build a single summed expression over one or more collections by expanding the
`$( ... )$` region once per item and joining the pieces with `+`. Each piece and
the whole sum are wrapped in parentheses, so the result is safe to use as an
operand (`2 * $(...)$`, `1 / $(...)$`). A reduction over an empty collection
expands to `0`.

```cpp
$$JAFF REDUCE specie_charge IN specie_charges
double total_charge = $($specie_charge$)$;
$$JAFF END
```

Output:

```cpp
// $$JAFF REDUCE specie_charge IN specie_charges
double total_charge = ((0) + (-1) + (1) + (2) + (1));
// $$JAFF END
```

Multiple collections can be reduced together (the i-th term takes the i-th
value of each; the first collection sets the number of terms). For example,
electron density from charge balance over the charged, non-electron species:

```cpp
$$JAFF REDUCE charged_specie_index_ne, charged_specie_charge_ne IN charged_specie_indices_ne, charged_specie_charges_ne
ne = $(($charged_specie_charge_ne$ * xn[$charged_specie_index_ne$]))$;
$$JAFF END
```

Output:

```cpp
// $$JAFF REDUCE charged_specie_index_ne, charged_specie_charge_ne IN charged_specie_indices_ne, charged_specie_charges_ne
ne = (((1 * xn[2])) + ((2 * xn[3])) + ((1 * xn[4])));
// $$JAFF END
```

Every `$( … )$` that uses a listed variable is expanded, however many share a
line; one that uses none of them is left as written. Lines are otherwise
copied (with `REPLACE` applied).

**Properties** (variable in brackets):

| Property                      | Variable                    |
| ----------------------------- | --------------------------- |
| `specie_charges`              | `specie_charge`             |
| `specie_masses`               | `specie_mass`               |
| `charge_truths`               | `charge_truth`              |
| `specie_masses_ne`            | `specie_mass_ne`            |
| `specie_charges_ne`           | `specie_charge_ne`          |
| `charge_truths_ne`            | `charge_truth_ne`           |
| `charged_specie_charges`      | `charged_specie_charge`     |
| `charged_specie_charges_ne`   | `charged_specie_charge_ne`  |
| `charged_specie_indices`      | `charged_specie_index`      |
| `charged_specie_indices_ne`   | `charged_specie_index_ne`   |
| `charged_specie_masses`       | `charged_specie_mass`       |
| `charged_specie_masses_ne`    | `charged_specie_mass_ne`    |
| `neutral_specie_indices`      | `neutral_specie_index`      |
| `neutral_specie_masses`       | `neutral_specie_mass`       |
| `photo_reaction_truths`       | `photo_reaction_truth`      |
| `photo_reaction_indices`      | `photo_reaction_index`      |
| `tmins`                       | `tmin`                      |
| `tmaxes`                      | `tmax`                      |

---

## Placeholders

Body lines are substituted as follows:

- `$name$` is replaced by the value of `name` (a SUB token, a GET property, a
  REPEAT variable, …). Unknown names are left alone.
- `$name+N$`, `$name-N$`, `$name*N$`, `$name/N$` apply arithmetic with an
  integer literal `N` to an integer value (`$nspec+1$`, `$specie_idx*2$`).
  Spaces around the operator are allowed. This works in `SUB`, `GET` and `HAS`;
  `/` is integer (floor) division, so `$nspec/2$` with five species gives `2`.
- In `REPEAT`, `$idx$` is the item index and takes the same arithmetic
  (`$idx+1$`, `$idx * 2$`, …). **Every** `$idx…$` token on a line counts
  towards the dimensionality check. `$IDX$` is the first index, through the
  arithmetic of the line's first `$idx…$` token if it has one.
- A `cse` variable's line `cse$idx$ = $cse$` must be the first line that
  expands (see [CSE](#cse-common-subexpression-elimination)).
- In `REDUCE`, `$( … )$` marks the part repeated per item.

A division by zero in a placeholder (`$x/0$`) is reported as a `ParserError` at
the directive.

---

## Options

Options go inside `{ … }` after the header. Each option is either

- `REPLACE { "regex" -> "replacement" … }` — one or more pairs;
- `USE name` — pull in a `DEFINE`d set; or
- `NAME { KEY -> value … }` — a group of settings (`THERMAL`, `SIGN`, …).

Groups do not nest, and keys are written inside their group — there is no
top-level `KEY -> value`. Values are bare words (no spaces, braces, commas or
quotes); anything else must be a double-quoted string.

### REPLACE

`REPLACE` applies Python regex substitutions to the generated output **after**
expansion. It is the template engine's most powerful feature — it lets you
target languages or naming conventions JAFF does not model natively.

```text
REPLACE {
    "pattern" -> "replacement"
    "pattern" -> "replacement"
}
```

- Patterns and replacements are **double-quoted** strings. Only `\"` is
  unescaped (to `"`); everything else reaches the regex untouched, so write
  `\d`, `\+`, `\[` as in any regex.
- Replacements may refer to capture groups with `\1`, `\2`, …. A reference to a
  group the pattern does not have is an error.
- Pairs apply **in order** to every output line of the block — expanded or
  copied — on the line's content after substitution (the indentation is never
  touched, but trailing whitespace is part of the content, so anchor with
  `\s*$` rather than `$` if the template may have trailing spaces). A pattern
  cannot match across lines.
- Patterns are compiled and checked before any rendering.

```cpp
$$JAFF REPEAT idx, specie IN species {
    REPLACE {
        "\+" -> "_plus"
        "-" -> "_minus"
    }
}
const int idx_$specie$ = $idx$;
$$JAFF END
```

Without `REPLACE`, the raw species names leak invalid identifier characters
(`idx_H+`, `idx_e-`). With the two pairs above:

```cpp
// $$JAFF REPEAT idx, specie IN species {
    // REPLACE {
        // "\+" -> "_plus"
        // "-" -> "_minus"
    // }
// }
const int idx_H = 0;
const int idx_e_minus = 1;
const int idx_H_plus = 2;
const int idx_He_plus_plus = 3;
const int idx_He_plus = 4;
// $$JAFF END
```

Capture groups work too:

```cpp
$$JAFF REPEAT idx, specie IN species { REPLACE { "(H)(\d+)" -> "\1_\2" } }
arr[$idx$] = "$specie$";
$$JAFF END
```

would turn `"H2"` into `"H_2"` and `"H2O"` into `"H_2O"` (this network has no such
species, so it is unchanged here).

REPLACE is handy for mapping JAFF's standard symbols (`tgas`, `nden[…]`,
`radeden[…]`, `photden[…]`, `rflux[…]`) onto your code's own names:

```cpp
$$JAFF REPEAT idx, rate IN rates { REPLACE { "tgas" -> "T" } }
state.rates[$idx$] = $rate$;
$$JAFF END
```

Output:

```cpp
// $$JAFF REPEAT idx, rate IN rates { REPLACE { "tgas" -> "T" } }
state.rates[0] = 5.8500000000000005e-11*std::sqrt(T)*std::exp(-157809.10000000001/T);
state.rates[1] = 2.6e-13*std::pow(T/10000.0, -0.80000000000000004);
state.rates[2] = nden[3];
// $$JAFF END
```

### DEFINE and USE

When several directives share the same options, name them once with `DEFINE`
and pull them in with `USE`:

```cpp
$$JAFF DEFINE micro {
    REPLACE {
        "nden\[\s*(\d+)\s*\]" -> "state.xn[\1]"
        "tgas" -> "state.T"      # comments run to the end of the line
    }
}

$$JAFF REPEAT idx, rate IN rates { USE micro }
rate[$idx$] = $rate$;
$$JAFF END

$$JAFF REPEAT idx, ode IN odes {
    USE micro
    REPLACE { "state.T" -> "T_gas" }
}
dydt[$idx$] = $ode$;
$$JAFF END
```

The `DEFINE` line produces no output other than its echo. Output (the second
block's local `REPLACE` pair runs after the one from `micro`, so `state.T` is
renamed to `T_gas`; long lines abridged):

```cpp
// $$JAFF REPEAT idx, rate IN rates { USE micro }
rate[0] = 5.8500000000000005e-11*std::sqrt(state.T)*std::exp(-157809.10000000001/state.T);
rate[1] = 2.6e-13*std::pow(state.T/10000.0, -0.80000000000000004);
rate[2] = state.xn[3];
// $$JAFF END

// $$JAFF REPEAT idx, ode IN odes {
    // USE micro
    // REPLACE { "state.T" -> "T_gas" }
// }
dydt[0] = 4.1207223003988966e-10*std::pow(T_gas, -0.80000000000000004)*state.xn[1]*state.xn[2] - ...;
...
dydt[3] = -state.xn[1]*std::pow(state.xn[3], 2);
dydt[4] = state.xn[1]*std::pow(state.xn[3], 2);
// $$JAFF END
```

Rules:

- `DEFINE name { … }` has no body and no `END`. The option block is required;
  a `DEFINE` without one is an error. The name must not be defined
  twice in a file, and a `DEFINE` cannot contain a `USE`.
- A `DEFINE` must appear **earlier in the file** than any `USE` of it.
- A block may have several `USE` lines, each naming a different set.
- The `USE` sets apply first, in the order written, then the block's own
  options. `REPLACE` lists **concatenate** (the `USE` pairs first); other groups
  **merge key by key**, with the local key winning.
- A `DEFINE` accepts any group some command takes; the group and keys are
  checked against the command that uses them, wherever they came from. Options
  that arrive through a `USE` are reported at the `USE` line in errors.

`END` takes no options.

---

## File settings

Two directives declare facts about the whole file. Each takes no options and no
`END`, may appear **at most once**, and must come **before any other
directive** (plain text such as a shebang line may come first).

### LANG

```text
$$JAFF LANG name
```

Names the target language, overriding the file extension. Accepted names
include the language's aliases (`c++`, `cpp`, `f90`, `py`, …) and are case
insensitive. Without `LANG` the extension decides. A file **with directives**
needs either `LANG` or a known extension; a file with **no directives** and no
known extension is copied unchanged.

```text
$$JAFF LANG f90
$$JAFF REPEAT idx, specie IN species
names($idx+1$) = "$specie$"
$$JAFF END
```

Output (a `.txt` file, so `LANG` is needed):

```text
! $$JAFF LANG f90
! $$JAFF REPEAT idx, specie IN species
names(1) = "H"
names(2) = "e-"
names(3) = "H+"
names(4) = "He++"
names(5) = "He+"
! $$JAFF END
```

### COMMENT

```text
$$JAFF COMMENT "token"
```

Sets the token used to [echo](#echo) directives. It is quoted, so `//`, `!!`
and `#` all work. The default is the comment token of the language. Use it for
non-code files — configuration or parameter files — where the language's own
token is wrong:

```text
#!/usr/bin/env python3
$$JAFF LANG python
$$JAFF COMMENT "#"
$$JAFF REPEAT idx, specie IN species
names[$idx$] = "$specie$"
$$JAFF END
```

Output:

```text
#!/usr/bin/env python3
# $$JAFF LANG python
# $$JAFF COMMENT "#"
# $$JAFF REPEAT idx, specie IN species
names[0] = "H"
names[1] = "e-"
names[2] = "H+"
names[3] = "He++"
names[4] = "He+"
# $$JAFF END
```

`LANG` and `COMMENT` are echoed like any other directive, using the resolved
token.

---

## Echo

By default each directive line is written to the output as a comment: the
directive's indentation, then the comment token, a space, and the rest of the
line. `$$JAFF END` is echoed too; `DEFINE` and `USE` lines appear as well, as
the examples above show. Every line of a multi-line header is echoed, with its
own indentation. Trailing whitespace is dropped, and a missing final newline is
added.

The echo can be switched off with the `echo_directives` argument of
`TemplateParser`. `jaffgen` sets it from `echo` in
[`jaffgen.toml`](jaffgen-toml.md) or the `--echo/--no-echo` flag (default on).

In Fortran, wrapping of long lines leaves comment lines alone.

---

## Errors

Every template problem is raised as a `ParserError` carrying the file name, the
line number and the line text, for example:

```text
ParserError: File: err1.cpp
Line number: 1
Line: $$JAFF REPEAT idx, rate IN rates { THERMAL { MODE -> dedt } }
REPEAT IN rates takes no option group THERMAL (accepted: none)
```

The whole file is **validated before anything renders**, so a typo in the last
directive fails immediately rather than after expensive code generation.
Detected problems include: syntax errors, a missing `END` or unbalanced
braces, an unknown command, property, variable, option group or key, a bad
option value, an invalid regex or `\N` reference, a `USE` of an unknown or
later `DEFINE`, a repeated `DEFINE` or `USE`, a `GET` entity that is not in the
network, of the wrong kind, or a reaction name shared by several mechanisms
(`HAS` only checks its kind; a missing entity simply gives `0`), options on `END`, a `$$JAFF` line inside a
body, a misplaced or repeated `LANG`/`COMMENT`, an unknown language, and no
language at all.

Errors in a multi-line header report the line they occur on; errors about a
body line (the CSE line, the number of `$idx$` tokens) report that body line.
Any other failure while a command renders (for example no radiation bands, or a
division by zero in `$x/0$`) is also a `ParserError`, at the directive, with the
original exception as its cause.
