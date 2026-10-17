---
tags:
    - User-guide
    - Code-generation
---

# Builder

Most of the time you generate code through the [normal workflow](../code-generation/index.md)
and never touch this page. `Builder` is one of the two **advanced**,
programmatic interfaces — for when you want to drive code generation from Python
yourself rather than from the CLI.

Of the two, `Builder` is the higher-level one. You hand it a loaded
[`Network`](../working-with-networks/network.md), name a template, and call
one method — `build()` — and a complete, runnable solver project is generated.
Everything it does happens through two choices you make:

- **which template** to use — `"python_solve_ivp"`, `"fortran_dlsodes"`, and so
  on. The template decides the language, the solver, and which files come out;
- **where to write** the result — any directory you like, or the current one if
  you say nothing.

If you instead want the individual rate, flux, or Jacobian strings to drop into
your own files by hand, reach for [`Codegen`](codegen.md) — the low-level engine
`Builder` is built on top of.

```python
from jaff import Network, Builder

net = Network("networks/demos/demo2.jet")

b = Builder(net)
b.build(template="python_solve_ivp")
```

`Builder` is exported from the top level, so `#!python from jaff import Builder`
and `#!python from jaff.codegen import Builder` are equivalent — use whichever
reads better alongside the rest of your imports.

**Output**

```text
Building network with template: python_solve_ivp
INFO     Preprocessing commons.py   -> ./commons.py
INFO     Preprocessing rates.py     -> ./rates.py
INFO     Preprocessing fluxes.py    -> ./fluxes.py
INFO     Preprocessing ode.py       -> ./ode.py
INFO     Copying main.py            -> ./main.py
Network built successfully using template 'python_solve_ivp'.
Output files are located in: /path/to/cwd
```

The four templated files (`commons.py`, `rates.py`, `fluxes.py`, `ode.py`) have
your network's species, rates, fluxes, and ODEs filled in; `main.py` is copied
through unchanged as a ready-to-edit driver.

---

## Building a project

You only ever touch two things: the constructor, which takes the network, and
`build()`, which takes the template and the destination.

```python
b = Builder(net)                       # remembers the network
path = b.build(template="python_solve_ivp", output_dir="generated/")
```

| Parameter    | Type            | Default              | Description                                          |
| ------------ | --------------- | -------------------- | ---------------------------------------------------- |
| `template`   | `str`           | `"python_solve_ivp"` | Which template to generate (see table below)         |
| `output_dir` | `str` or `None` | `None`               | Where to write the files; `None` → current directory |

`build()` returns the path it wrote to (the string you passed as `output_dir`,
or the current working directory when you passed nothing), so you can chain it
straight into the next step:

```python
out = b.build(template="python_solve_ivp", output_dir="generated/")
print("Project written to", out)
```

One `Builder` can build as many times as you like — different templates,
different directories — from the same network:

```python
b = Builder(net)
b.build(template="python_solve_ivp", output_dir="py/")
b.build(template="fortran_dlsodes",  output_dir="f90/")
```

<!-- prettier-ignore -->
!!! warning "A wrong template name stops the program"
    If the template doesn't exist, `Builder` prints the list of templates it
    *does* know about and exits the process (`SystemExit`, code `1`) — it does
    not raise a catchable exception:
    ```text
    Error: Template 'pyhton_solve_ivp' not found. Available templates are:
    fortran_dlsodes
    kokkos_ode
    microphysics
    python_solve_ivp
    ```
    So check the spelling against the table below before scripting a build into
    a larger pipeline.

---

## Available templates

Each template is a self-contained recipe: a target language, a solver, and the
set of files it emits.

| Template           | Language     | Solver                      | What you get                                               |
| ------------------ | ------------ | --------------------------- | ---------------------------------------------------------- |
| `python_solve_ivp` | Python       | `scipy.integrate.solve_ivp` | `commons.py`, `rates.py`, `fluxes.py`, `ode.py`, `main.py` |
| `fortran_dlsodes`  | Fortran 90   | DLSODES (ODEPACK)           | Fortran source + driver                                    |
| `kokkos_ode`       | C++17        | VODE (header-only-integrators) | `chemistry_ode.hpp` (RHS + Jacobian), driver, `CMakeLists.txt` |
| `microphysics`     | C++          | AMReX Microphysics          | AMReX-compatible RHS + Jacobian                            |

`python_solve_ivp` is the default and the quickest way to _see_ a network
integrate — no compiler needed, just `python main.py`.

Despite its name, `kokkos_ode` does not use Kokkos: it integrates with the VODE
(BDF) solver from the
[header-only-integrators](https://github.com/BenWibking/header-only-integrators)
library, which its `CMakeLists.txt` downloads on first configure. Build it with
clang (`cmake -S . -B build -DCMAKE_CXX_COMPILER=clang++`); the library
currently uses a C++ construct that recent GCC versions reject.

---

## The Python template, end to end

Pick a real network, point `build()` at a fresh directory, and run the driver:

```python
from jaff import Network, Builder

net = Network("networks/GOW/GOW.jet")
b = Builder(net)
b.build(template="python_solve_ivp", output_dir="generated/")
```

The generated `ode.py` holds the right-hand side, with your species and reactions
already wired in (excerpt):

```python
def ode(t, y, tgas, av):
    k = rates(tgas, av)
    f = fluxes(k, y)

    dydt = [0.0] * nspecs
    dydt[idx_H]  = -flux[0] + flux[2]
    dydt[idx_H2] =  flux[0] - flux[1]
    # ...
    return dydt
```

and `main.py` is the driver you actually run and tweak (excerpt):

```python
from scipy.integrate import solve_ivp
from ode import ode, nspecs

y0 = [1e-4] * nspecs
sol = solve_ivp(lambda t, y: ode(t, y, tgas=10.0, av=1.0),
                t_span=(0, 1e13), y0=y0, method="LSODA")
```

```bash
cd generated/
python main.py
```

---

## Adding a new template

`Builder` only runs templates it can find, but the set is open — adding your own
language or solver is the main reason to reach past the built-ins. A template is
two pieces that share one name:

- **the template files** in `src/jaff/templates/preprocessor/<name>/` — ordinary
  source files with `PREPROCESS_` markers where generated code should land;
- **a plugin** at `src/jaff/plugins/<name>/plugin.py` exposing a `main()` that
  fills those markers in.

The `<name>` is exactly the string you later pass to `build(template=...)`;
`Builder` matches the directory and the plugin module by that name.

### 1. Write the template files

Each file is the real thing you want to emit, with a marked-off hole. A line of
the form `<comment> PREPROCESS_<KEY>` opens a block and `<comment> PREPROCESS_END`
closes it; everything between is replaced by the generated code for `<KEY>`. The
marker lines survive, so the output stays re-processable.

```python hl_lines="9 11"
# src/jaff/templates/preprocessor/my_python/rates.py
import numpy as np
from commons import nreactions


def get_rates(tgas, crate, av):
    k = np.zeros(nreactions)

    # PREPROCESS_RATES

    # PREPROCESS_END
    return k
```

Files with no markers — a `main.py` driver, a `CMakeLists.txt` — need no special
treatment; they are copied through unchanged.

### 2. Write the plugin

The plugin builds a [`Codegen`](codegen.md) for the network, asks it for the
strings each marker needs, and hands the whole lot to the `Preprocessor`. The
list of files and the list of `{KEY: string}` dictionaries line up positionally:

```python
# src/jaff/plugins/my_python/plugin.py
from jaff import Codegen, Preprocessor


def main(network, path_template, path_build=None):
    p = Preprocessor()
    cg = Codegen(network=network, lang="python")

    p.preprocess(
        path_template,
        ["commons.py", "rates.py"],
        [{"COMMONS": cg.get_commons()}, {"RATES": cg.get_rates_str()}],
        comment="#",
        path_build=path_build,
    )
```

`Builder` calls `main(network, path_template=..., path_build=...)`, so keep that
signature. Match `comment` to the language's comment marker (`#`, `//`, `!!`) so
the `PREPROCESS_` lines are recognised.

### 3. Build it

Once both pieces exist under the same name, it's a first-class template:

```python
b.build(template="my_python", output_dir="generated/")
```

<!-- prettier-ignore -->
!!! tip "Going deeper"
    For the marker grammar, the per-language comment styles, and the full set of
    `Codegen` strings you can inject, see
    [Adding template properties](../../development/adding-template-properties.md).
