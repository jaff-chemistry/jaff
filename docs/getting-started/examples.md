---
tags:
    - Introduction
    - Examples
icon: phosphor/book-open
---

# Example Notebooks

The `examples/` directory holds a series of [marimo](https://marimo.io)
notebooks that walk through JAFF from loading a network to running generated
solvers. Each one explains every step as it goes, and later notebooks build on
earlier ones, so they read best in order.

| Notebook                     | Network(s)                 | What it covers                                                                                                                                         |
| ---------------------------- | -------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `01_network_basics.py`       | `demo1`                    | Loading a network, the mixed input formats in one file, the load pipeline, species, catalogue queries, stoichiometry and a hand-made charge check      |
| `02_reactions_and_rates.py`  | `demo1`, `h_photo`         | Reaction identity and types, symbolic rates, printing a rate in seven languages, flux expressions, rate plots, multi-range rates, photo cross-sections |
| `03_comparing_networks.py`   | `demo1`, `demo2`           | Matching species and reactions by their serialized form, and why shared reactions can still have different rates                                      |
| `04_codegen_jaffgen.py`      | `demo2`                    | Template collections and directives; generating, building and running a Python and a Fortran solver with `jaffgen`                                     |
| `05_kokkos_cpp.py`           | `COthin`                   | Generating the `kokkos_ode` C++ solver, building it with CMake, running it, and its element/charge conservation report                                 |
| `06_rate_tables.py`          | `uclchem_small_gas`        | Tabulating rate coefficients with `write_table`, the HDF5 layout, interpolation accuracy, and the accuracy/size trade-off                               |

## Running them

marimo notebooks are plain Python files. Open one in the interactive editor,
where cells re-run automatically when their inputs change:

```bash
uv run marimo edit examples/01_network_basics.py
```

or run one straight through without a browser. A cell that raises, including
a shell command whose exit code is not 0, makes the run exit with an error:

```bash
uv run python examples/01_network_basics.py
```

Generated files are written to a temporary directory, so running a notebook
leaves the repository untouched.

## Requirements

Notebooks `01`–`03` and `06` need only JAFF itself. The code-generation
notebooks shell out to external tools:

- `04_codegen_jaffgen.py` builds the Fortran solver with `gfortran` and `make`;
- `05_kokkos_cpp.py` needs `cmake`, `clang++` and network access, because CMake
  downloads the [header-only-integrators](https://github.com/BenWibking/header-only-integrators)
  library on first configure. Use clang: the library currently uses a C++
  construct that recent GCC versions reject.

## Tested in CI

`tests/test_example_notebooks.py` runs every notebook headless as part of the
regular test suite and requires it to exit 0. Notebooks `04` and `05` run only
on Linux, and are skipped when a tool they need is missing.
