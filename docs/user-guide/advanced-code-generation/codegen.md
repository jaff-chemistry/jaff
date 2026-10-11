---
tags:
    - User-guide
    - Code-generation
---

# Codegen

Most of the time you generate code through the [normal workflow](../code-generation/index.md),
and when you want a whole project in one call you reach for
[`Builder`](builder.md). `Codegen` sits one level below both: it is the engine
they are built on, and it hands you the **pieces** — the actual lines of code
for rates, fluxes, ODEs, the Jacobian — as strings you drop into your own files.

Use it when writing `Builder`'s plugin or assembling a file your own way.
You give `Codegen` two things up front — a parsed
[`Network`](../working-with-networks/network.md) and a **target language** — and
from then on every `get_*_str()` call returns a ready-to-paste block in that
language.

```python
from jaff import Network
from jaff.codegen import Codegen

net = Network("networks/GOW/GOW.jet")
cg = Codegen(net, lang="python")

print(cg.get_rates_str())     # k[0] = ... , k[1] = ... , ...
print(cg.get_ode_str())       # f[0] = ... , f[1] = ... , ...
```

`Codegen` is exported from the top level as well, so
`#!python from jaff import Codegen` works too.

Every quantity comes in **two forms**, and the rest of this page leans on that
split:

- a **formatted string** — `get_rates_str()`, `get_ode_str()`, … — the
  paste-ready text, with brackets, indices, and line terminators already in the
  target language;
- a **structured form** — `get_indexed_rates()`, `get_indexed_odes()`, … — the
  same expressions as `(index, expression)` pairs, for when you'd rather format
  them yourself.

The string methods are what you want genrally; the indexed ones are covered at
the end.

---

## Choosing a target language

The language you pass to the constructor decides four things at once: the array
brackets, whether indices start at `0` or `1`, the comment marker, and the
assignment operator. Pick it with a name or any common alias.

```python
Codegen(network, lang="c++")
```

| Parameter       | Type      | Default | Description                                 |
| --------------- | --------- | ------- | ------------------------------------------- |
| `network`       | `Network` | —       | Parsed chemical reaction network            |
| `lang`          | `str`     | `"c++"` | Target language — name or alias (see table) |

An unsupported `lang` raises `InvalidLanguageError` listing what _is_ supported,
so a typo fails loudly at construction rather than producing wrong code.

Bracket and token formatting are **not** fixed at construction — each
`get_*_str` method accepts per-call overrides (`brac_format`, `matrix_format`,
`assignment_op`, `line_end`) that apply only to that call and reset afterwards.
See [the Jacobian](#the-jacobian-get_jacobian_str) for a `matrix_format` example.

| Alias(es)           | Canonical | Brackets | Index base | Comment | Assignment |
| ------------------- | --------- | -------- | ---------- | ------- | ---------- |
| `c++`, `cpp`, `cxx` | `cxx`     | `[]`     | 0          | `//`    | `=`        |
| `c`                 | `c`       | `[]`     | 0          | `//`    | `=`        |
| `fortran`, `f90`    | `fortran` | `()`     | 1          | `!`     | `=`        |
| `python`, `py`      | `python`  | `[]`     | 0          | `#`     | `=`        |
| `rust`, `rs`        | `rust`    | `[]`     | 0          | `//`    | `=`        |
| `julia`, `jl`       | `julia`   | `[]`     | 1          | `#`     | `=`        |
| `r`                 | `r`       | `[]`     | 1          | `#`     | `<-`       |

The index base is the one to watch: Fortran, Julia, and R are 1-based, so the
same network emits `k(1)` there and `k[0]` in C++. You don't need to set it by hand —
the language already knows.

---

## The shared knobs

Almost every `get_*_str()` method takes the same handful of formatting
arguments, so they're explained once here. Each defaults to "do what the
language says", so you only pass them to deviate.

| Argument        | What it does                                                                       | Empty / `-1` means                  |
| --------------- | ---------------------------------------------------------------------------------- | ----------------------------------- |
| `idx_offset`    | Number added to every array subscript                                              | `-1` → the language's base          |
| `brac_format`   | Override the 1-D brackets for this one call                                        | `""` → the language's brackets      |
| `assignment_op` | Override the assignment operator                                                   | `""` → the language's `=`/`<-`      |
| `line_end`      | Override the statement terminator                                                  | `""` → the language's `;` / nothing |
| `use_cse`       | Factor out repeated sub-expressions (see [CSE](#common-subexpression-elimination)) | default `True`                      |
| `cse_var`       | Name prefix for the CSE temporaries                                                | per method (`x`, `cse`, `rcse`)     |

<!-- prettier-ignore -->
!!! warning "`idx_offset` defaults are not all the same"
    The lower-level emitters (`get_commons`, `get_rates_str`,
    `get_flux_expressions_str`, `get_ode_expressions_str`) default
    `idx_offset=-1`, i.e. *use the language base* — so Fortran gets 1-based
    subscripts automatically. The assembled-RHS emitters (`get_ode_str`,
    `get_rhs_str`, `get_radode_str`, `get_jacobian_str`) default `idx_offset=0`,
    a literal zero. If you emit those for a 1-based language, pass
    `idx_offset=-1` explicitly to pick the language base back up.

The methods that emit CSE temporaries also take a **type-declaration prefix** so
the temporaries are valid declarations in typed languages — it's called
`var_prefix` on `get_rates_str`/`get_jacobian_str` and `def_prefix` on
`get_ode_str`/`get_rhs_str`/`get_radode_str`. Left empty it becomes the
language's `const double `/`f64 `/… ; for an untyped language (Python, Fortran,
R) it's empty.

---

## The generation pipeline

The methods mirror the math, and they build on each other in the same order a
solver evaluates them: **indices → rates → fluxes → ODEs**, with the Jacobian
off to the side. You can call any one in isolation; they don't share state.

The running examples below are excerpts, shown in the language that reads
clearest for each step.

### Species indices — `#!python get_commons()`

The starting point: one named constant per species mapping its index name
(`fidx`) to its slot, then the network sizes `nspecs` and `nreactions`. Every
other block indexes against these.

```python
cg = Codegen(net, lang="c++")
print(cg.get_commons(definition_prefix="const int "))
```

```cpp
const int idx_h  = 0;
const int idx_h2 = 1;
const int nspecs = 2;
const int nreactions = 5;
```

`definition_prefix` is prepended to every line (here `"const int "`); the index
names themselves come from each specie's `fidx`, so they're already lowercase
and identifier-safe (`H2` → `idx_h2`, `H+` → `idx_hj`).

### Reaction rates — `#!python get_rates_str()`

One assignment per reaction, `k[i] = <rate>`, with CSE temporaries emitted first
so later lines can reference them.

```python
print(cg.get_rates_str(rate_variable="k"))
```

```cpp
const double x0 = exp(-1.0/tgas);
k[0] = 1.8e-11 * x0;
k[1] = photorates(1, G0, av, chi);
k[2] = 3.2e-17;
```

Photochemical rates carry a `photorates($IDX$, ...)` placeholder; `get_rates_str`
substitutes the concrete reaction index (`$IDX$` → `1`) so the line compiles as
written. The CSE prefix here defaults to `x` (giving `x0`, `x1`, …) rather than
the `cse` used elsewhere.

### Fluxes — `#!python get_flux_expressions_str()`

Each reaction's flux: its rate times the densities of its reactants.

```python
print(cg.get_flux_expressions_str(rate_var="k", species_var="y"))
```

```python
flux[0] = k[0] * y[idx_h] * y[idx_h2]
flux[1] = k[1] * y[idx_h]
```

Pass `idx_prefix="idx_"` only if you want an _extra_ prefix on top of the
already-prefixed `fidx`; usually you leave it empty.

### Per-species ODEs — `#!python get_ode_expressions_str()`

The density derivatives written as **sums over a pre-computed `flux` array** —
reactants subtract, products add. This is the cheap form: it assumes `flux[...]`
already exists in the generated code.

```python
print(cg.get_ode_expressions_str())
```

```python
dy[idx_h]  = - flux[0] + flux[2]
dy[idx_h2] = + flux[0] - flux[1]
```

The derivative array is named by `derivative_prefix` + `species_var` (default
`"d"` + `"y"` → `dy`), or set `derivative_var` outright. The subscript is the
species `fidx`, so it lines up with the constants from `get_commons()`.

### The full right-hand side — `#!python get_ode_str()` and `#!python get_rhs_str()`

`get_ode_str()` is the self-contained version of the species ODEs: it inlines
the rate expressions instead of leaning on a `flux` array, and runs CSE across
every species at once.

```python
print(cg.get_ode_str(ode_var="f"))
```

```cpp
const double cse0 = exp(-1.0/tgas);
const double cse1 = 1.8e-11 * cse0;
f[0] = -cse1 * y[0] * y[1] + 3.2e-17 * y[2];
f[1] =  cse1 * y[0] * y[1] - 2.7e-10 * y[1];
```

`get_rhs_str()` is the same idea taken to the whole state vector — species ODEs,
then the energy derivative, then (optionally) radiation — with CSE shared across
_all_ of them so sub-expressions are factored out jointly. The output comes in a
fixed order:

1. CSE temporaries
2. `f[0]` … `f[N-1]` — species density ODEs
3. `f[N]` — energy time-derivative `dE/dt`
4. `f[N+1]` … — radiation ODEs (only when `radiation=True`)

```python
cg.get_rhs_str(
    energy="volumetric",   # or "specific", "per_particle", "molar"
    radiation=True,        # append the radiation moment ODEs
    rad_order=0,           # radiation moment closure order (0–3)
)
```

### Energy derivative — `#!python get_dedt()`

Just the energy time-derivative `dE/dt`, as a single expression with no
assignment or line terminator — handy when you splice it into your own line.

```python
cg.get_dedt(energy="specific")   # erg/g/s, normalised by mass density
cg.get_dtdt()                    # dT/dt in K/s, see get_dtdt
```

### Radiation ODEs — `#!python get_radode_str()`

The radiation moment ODEs on their own, formatted exactly like `get_ode_str()`.
`order` selects the moment closure (`0`–`3`); the CSE prefix defaults to `rcse`.

### The Jacobian — `#!python get_jacobian_str()`

The analytical Jacobian `J[i][j] = ∂f_i/∂y_j`, computed symbolically and emitted
**sparse** — only the non-zero elements appear, ready for a sparse solver.

```python
print(cg.get_jacobian_str(jac_var="J"))
```

```cpp
const double cse0 = exp(-1.0/tgas);
J[0][0] = -1.8e-11 * cse0 * y[1];
J[0][1] = -1.8e-11 * cse0 * y[0];
J[1][1] = -2.7e-10;
```

A few specifics:

- **`thermal="dedt"`** (or **`"dtdt"`**) adds the thermal equation's row and column, coupling the
  chemistry to temperature through the network's EOS. The default is `"none"`.
- Rate functions SymPy can't differentiate symbolically (e.g.
  `#!python photorates(...)`) become named partial calls —
  `#!python photorates_partial_0(...)`, where the suffix is the argument
  differentiated against — for you to supply.
- **`matrix_format`** picks the 2-D indexing style independently of the 1-D
  brackets. Use `","`-variants for a single subscript and the plain/doubled
  forms for nested ones:

    | `matrix_format` | Emits     |
    | --------------- | --------- |
    | `"[]"`          | `J[i][j]` |
    | `"[,]"`         | `J[i, j]` |
    | `"(,)"`         | `J(i, j)` |
    | `"{,}"`         | `J{i, j}` |

    (`"()"`, `"[][]"`, `"<>"`, and their `,`-variants follow the same pattern.) An
    unsupported value raises `InvalidLanguageError`.

---

## Common subexpression elimination

`exp(-1.0/tgas)` showing up in three rates is wasted work. With `use_cse=True`
(the default on every method that supports it) JAFF runs SymPy's
[`cse`](https://docs.sympy.org/latest/modules/simplify/simplify.html#sympy.simplify.cse_main.cse)
over the expressions, lifts each repeated piece into a temporary (`cse0`, `cse1`,
…), and emits those before the assignments that use them. Unused temporaries are
pruned, so nothing dead reaches the output.

```python
cg.get_rates_str(use_cse=False)               # every rate fully expanded inline
cg.get_rates_str(use_cse=True, cse_var="tmp")  # shared pieces lifted to tmp0, tmp1, …
```

Two kinds of rate are left out of CSE because they can't be simplified
symbolically: rates stored as raw code strings, and `photorates(...)` calls
(the `$IDX$` placeholder can't be folded into a shared sub-expression).

---

## The structured form: `#!python get_indexed_*()`

When you don't want a finished string — you're feeding the expressions into a
template engine, renumbering them, or post-processing — every emitter has an
`get_indexed_*` sibling that returns the expressions as data instead of text.

| String method                | Structured sibling               |
| ---------------------------- | -------------------------------- |
| `get_rates_str()`            | `get_indexed_rates()`            |
| `get_flux_expressions_str()` | `get_indexed_flux_expressions()` |
| `get_ode_expressions_str()`  | `get_indexed_ode_expressions()`  |
| `get_ode_str()`              | `get_indexed_odes()`             |
| `get_rhs_str()`              | `get_indexed_rhs()`              |
| `get_radode_str()`           | `get_indexed_radodes()`          |
| `get_jacobian_str()`         | `get_indexed_jacobian()`         |

The rate/ODE/RHS/Jacobian variants return a dict with two keys —
`"extras"["cse"]` (the CSE temporaries) and `"expressions"` (the main
expressions) — each an `IndexedList` of `(index, expression_string)` pairs. The
flux and ode-expression variants return that `IndexedList` directly. The string
methods are just thin formatters over these.

```python
ir = cg.get_indexed_rates()
for idx, expr in ir["expressions"]:
    print(idx, expr)        # ([0], '1.8e-11*x0'), ([1], 'photorates($IDX$, ...)'), ...
```

---

## Inspecting the language tokens

Each language's syntax lives on a `Language` singleton. Resolve one by alias and
read its attributes — useful if you're matching JAFF's output to a hand-written
file.

```python
from jaff.codegen import Language

py = Language("python")
py.comment                       # '#'
Language("fortran").idx_offset   # 1

Language.registered()            # every registered language singleton
Language.comments()              # {'#', '//', '!'}
```

A `Language` carries `lb`/`rb`, `mlb`/`mrb`, `sep`, `assignment_op`, `line_end`,
`code_gen`, `idx_offset`, `comment`, `types`, and `extras`. Add a new target
language by subclassing `Language` with those class attributes — no central
edits required.

---

## Full example

Stitching the pieces into a minimal C++ source file — `Codegen` gives the
bodies, you own the scaffolding:

```python
from jaff import Network
from jaff.codegen import Codegen

net = Network("networks/GOW/GOW.jet")
cg = Codegen(net, lang="c++")

with open("output/network.cpp", "w") as f:
    f.write("#include <cmath>\n\n")
    f.write("// Species indices\n")
    f.write(cg.get_commons(definition_prefix="const int "))
    f.write("\nvoid rates(double* k, double tgas) {\n")
    f.write(cg.get_rates_str())
    f.write("}\n\n")
    f.write("void rhs(double* f, double* y, double tgas) {\n")
    f.write(cg.get_rhs_str())
    f.write("}\n")
```

If you find yourself writing the whole scaffold too, that's the cue to step up
to [`Builder`](builder.md) — it does exactly this from a template.
