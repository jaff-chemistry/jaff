---
tags:
    - Development
icon: phosphor/stack
---

# Codebase Structure

This page maps the `src/jaff` source tree, explains what each package owns, and shows how data flows through the library from a raw network file to generated code.

## Package Map

```
src/jaff/
├── core/                       # Domain model
│   ├── network/                # network.py (Network — main entry point)
│   │   ├── _spec.py            # NetworkSpec — normalized Network params
│   │   └── _args.py            # NetworkArgs — raw CLI arg accumulator
│   ├── reaction/               # reaction.py (Reaction) · reactions.py (Reactions)
│   ├── species/                # specie.py (Specie) · species.py (Species)
│   ├── elements/               # element.py (Element) · elements.py (Elements)
│   ├── parsers/                # File parsers (network + auxiliary)
│   │   ├── network/            # Multi-format network file parser
│   │   │   ├── _engine.py      # NetworkParser — drives format plugins
│   │   │   ├── _typing/        # parsedListProps, krome/prizmoFormatProps
│   │   │   └── _formats/       # One subpackage per format
│   │   │       ├── _parser.py  # Parser ABC (the only ABC) + register / all_parsers
│   │   │       ├── _record.py  # Record / ParsedRecord / ParseResult
│   │   │       ├── __init__.py # re-exports Parser, register, all_parsers, ...
│   │   │       ├── krome/      # parser.py + header.py · var.py · reaction.py handlers
│   │   │       ├── prizmo/     # parser.py + vars.py · reaction.py handlers
│   │   │       ├── udfa/       # parser.py + reaction.py handler
│   │   │       ├── uclchem/    # parser.py + reaction.py handler
│   │   │       └── kida/       # parser.py + reaction.py handler
│   │   └── auxiliary_func/     # .jfunc auxiliary function parser
│   │       ├── _engine.py      # AuxiliaryFunctionParser
│   │       └── _typing/        # AuxiliaryFunctionsDict
│   └── _typing/                # Shared core TypedDicts (Network/Element/Reaction)
│
├── physics/                    # Symbolic ODE/flux generation + physics helpers
│   ├── _equations.py           # get_sfluxes, get_sodes, get_sradodes
│   ├── photo_reactions/        # Photochemistry: cross sections, radiation, shielding
│   │   ├── _photochemistry.py  # get_xsec / get_verner_xsec / shielding — lookups
│   │   ├── _radiation.py       # Radiation moment equations
│   │   ├── _typing/            # TypedDicts (XsecsProps, ...)
│   │   └── shielding/          # Shielding-function registry (@_register, by reaction metadata)
│   │       ├── _base.py        # ShieldingFunction ABC (name, reaction attrs)
│   │       ├── global_/        # Global models, reaction=None (e.g. leiden.py)
│   │       └── H2__PHOTON__H_H/ # Local H2 self-shielding (db1996, hg2015) + shared _utils
│   ├── _typing/                # TypedDicts (Numeric, ...)
│   └── constants.py            # Physical constants (astropy Quantities)
│
├── plotting/                   # Publication-style seaborn plotting
│   ├── _api.py                 # plot_rates / plot_xsecs — free functions (reactions, exprs, arrays)
│   ├── plotter.py              # Plotter.render_series — seaborn-objects renderer
│   ├── _theme.py               # seaborn theme, palettes, scoped/global application
│   ├── _frames.py              # tidy DataFrame builders
│   ├── _units.py               # energy/xsec unit conversion + axis labels
│   └── _xsec.py                # trim / dynamic-scale helpers
│
├── codegen/                    # Code generation pipeline
│   ├── codegen.py              # SymPy → C/C++/Fortran/Python/Rust/Julia/R
│   ├── preprocessor.py         # Template marker substitution
│   ├── builder.py              # Plugin-based orchestration
│   └── _template_engine.py     # JAFF directive rendering
│
├── io/                         # Serialization and logging
│   ├── _io.py                  # .jaff gzip-JSON read/write; data table export
│   └── _logger.py              # JaffLogger + progress bars
│
├── config/                     # Package-wide path constants
│   └── _config.py              # SRC_DIR, DATA_DIR, XSECS/SHIELDING dirs, ...
│
├── drivers/                    # Config / data format adapters
│   ├── toml.py                 # TOML config reader
│   ├── csv.py                  # CSV I/O
│   ├── hdf5.py                 # HDF5 I/O
│   ├── sqlite.py               # SQLite I/O
│   └── pooch.py                # Download/cache remote cross-section data files
│
├── cli/                        # Command-line entry points (Typer)
│   ├── _helper.py              # Shared argument helpers (funcfile_arg)
│   ├── jaffgen/                # jaffgen — template-driven code generation
│   │   ├── _engine.py          # JaffGen pipeline + Typer `generate` command
│   │   ├── _structs.py         # State / ResolvedPath
│   │   └── _config_table.py    # [[table]] config → HDF5/CSV output
│   └── jaffx/                  # jaffx — network inspection / export
│       └── _engine.py          # JaffX handlers + nested Typer commands
│
├── plugins/                    # Named solver plugins
│   ├── python_solve_ivp/       # SciPy solve_ivp wrapper
│   ├── fortran_dlsodes/        # Fortran DLSODES solver
│   ├── kokkos_ode/             # Kokkos GPU ODE solver
│   └── microphysics/           # AMReX microphysics driver
│
├── templates/                  # Source templates consumed by plugins
│   ├── generator/<name>/       # JAFF directive template files
│   └── preprocessor/<name>/    # Marker substitution templates
│
├── types/                      # Base data structures
│   ├── _catalogue.py           # Catalogue[T] — O(1) list + dict lookup
│   ├── _vector.py              # Typed numeric container
│   ├── _indexed.py             # IndexedList / IndexedValue
│   └── _hdf5.py                # HDF5 type helpers
│
├── common/                     # Shared utilities
│   ├── _helper.py              # Element/mass table loading
│   ├── _integrators.py         # Dependency resolution (DFS)
│   ├── _sympy_json.py          # Versioned SymPy ↔ JSON encoding
│   ├── _fastlog.py             # Fast structured logging
│   └── _welcome.py             # MOTD / version banner
│
├── errors/
│   └── _parser.py              # ParserError hierarchy
│
├── data/                       # Raw data assets
│   ├── .downloads/             # Compressed, hash-verified pooch downloads + registry.txt (not bundled)
│   ├── atom_mass.csv           # Element mass table (bundled)
│   ├── xsecs/                  # Photo cross-section data (downloaded via drivers/pooch.py, not bundled)
│   │   ├── leiden.hdf5         # Leiden PDR cross sections (one group per reaction)
│   │   ├── norad.hdf5          # NORAD/OP ground-state photoionisation
│   │   └── verner_1996.csv     # Verner (1996) analytic-fit parameters
│   └── shielding/              # Line-shielding tables (downloaded via drivers/pooch.py, not bundled)
│       └── leiden.hdf5         # Leiden line shielding (one group per reaction)
│
├── db/                         # Prebuilt SQLite database
│   └── jaff.db                 # Mass + photo cross-section (Leiden/NORAD + Verner) tables, built from data/
│
└── _utils/                     # Standalone maintenance scripts
    ├── generate_mass_table.py          # Build mass tables in jaff.db from data/atom_mass.csv
    ├── download_nahar_xsecs.py         # Download NORAD/OP ground-state photoionisation .dat files
    ├── collapse_xsecs_hdf5.py          # Merge per-reaction files into leiden.hdf5 / norad.hdf5
    ├── split_xsecs_photodecay.py       # Split source diss/ion datasets into the photodecay channel
    ├── generate_photo_xsecs_table.py   # Build photo_reaction_cross_sections table in jaff.db
    ├── generate_ion_xsecs_table.py     # Build verner_cross_sections table in jaff.db
    ├── compress_hdf5.py                # Gzip local data files for the mirror + registry lines
    └── build_shielding_hdf5.py         # Collapse Leiden shielding tables into shielding/leiden.hdf5
```

Downloaded data files are cached, compressed and hash-checked, under
`src/jaff/data/.downloads/`. `drivers/pooch.py` then installs a copy at the path
shown above: HDF5 files are rewritten uncompressed and contiguous so they load
quickly (the Leiden cross sections take ≈365 MB), other files are copied. A file
you place at an install path yourself is never overwritten; JAFF logs a warning
instead. Set `JAFF_OFFLINE=1` to skip all downloads.

To publish regenerated data, edit or rebuild the (uncompressed) files under
`src/jaff/data/` locally, then write compressed copies plus their registry
lines with
`python -m jaff._utils.compress_hdf5 src/jaff/data/xsecs/leiden.hdf5 ... --outdir upload/ --registry upload/registry.txt`
and upload the contents of `upload/` to the mirror.

## Architecture Diagram

```mermaid
%%{init: {"flowchart": {"useMaxWidth": false}}}%%
flowchart TD
    subgraph input_sg ["Input"]
        NF["Network file\nKROME · PRIZMO · UDFA\nKIDA · UCLChem · .jaff"]
        JF[".jfunc\nauxiliary functions"]
        CFG["jaffgen.toml / CLI"]
    end

    subgraph parse_sg ["Parsing  —  core.parsers"]
        NE["NetworkParser\nauto-detect format\nformat plugins → dicts"]
        AE["AuxiliaryParser\n@var / @function\nSymPy expressions"]
    end

    subgraph model_sg ["Domain Model  —  core"]
        NET["Network\nassemble · validate\nSpecies · Reactions · Elements"]
    end

    subgraph codegen_sg ["Code Generation"]
        EQ["physics\nsfluxes · sodes · sradodes"]
        CG["Codegen\nSymPy → C · C++ · F90\nPy · Rust · Julia · R"]
        TP["TemplateParser  —  jaffgen path\ntemplates/generator/\nSUB · REPEAT · REDUCE directives"]
        OUT_G["Generated output files"]
        PP["Preprocessor  —  builder path\ntemplates/preprocessor/\n!! KEY marker substitution"]
        BL["Builder\nplugin dispatch"]
        OUT_B["Plugin output files"]
    end

    NF --> NE
    JF --> AE
    CFG --> NE

    NE --> NET
    AE --> NET

    NET --> EQ --> CG

    CG --> TP --> OUT_G
    CG --> PP --> BL --> OUT_B
```

## Data Flow — End to End

The table below traces a single `jaffgen` invocation from command line to output files.

| Step | Component                                | What happens                                                                                       |
| ---- | ---------------------------------------- | -------------------------------------------------------------------------------------------------- |
| 1    | `cli/jaffgen/_engine.py`                 | Parse CLI args (Typer), read `jaffgen.toml`, resolve config: CLI > jaffgen.toml > Network defaults |
| 2    | `core/parsers/network/_engine.py`        | Auto-detect format via registered plugins; convert each reaction line to a `parsedListProps` dict  |
| 3    | `core/parsers/auxiliary_func/_engine.py` | Parse `.jfunc` file (if present); resolve `@var`/`@function` blocks into SymPy expressions         |
| 4    | `core/network/network.py`                | Build `Species`, `Reactions`, `Elements` catalogues; validate duplicates, sinks, isomers           |
| 5    | `physics/_equations.py`                  | Compute symbolic fluxes (`sfluxes`) and ODE RHS (`sodes`) using SymPy                              |
| 6    | `codegen/codegen.py`                     | Translate SymPy expressions into assignment strings for the chosen language                        |
| 7    | `codegen/preprocessor.py`                | Walk template files; replace `!! PREPROCESS_KEY … !! PREPROCESS_END` blocks with generated strings |
| 8    | `codegen/builder.py`                     | Invoke the named plugin's `#!python main()` to write final output files to the build directory     |

## Key Design Decisions

**One `Parser` subclass per format, plain handler objects.**
Each network format is a single `Parser` subclass (the only ABC, in `_formats/_parser.py`) living in its own subpackage's `parser.py` under `core/parsers/network/_formats/`. A class registers itself with the `@register` decorator; `NetworkParser` discovers all parsers via `all_parsers()`, ordered by each parser's `priority` (not file or import order). A `Parser` owns a list of plain (no base class) **handlers** — one per line-type — each exposing a static `global_re`/`local_re` pair for detection/extraction and a `parse` (reaction) or `apply` (directive) method; the engine buckets raw `Record`s by owning parser and calls each parser's `process(records)`, which returns a `ParseResult` of `ParsedRecord`s + globals. Adding a new format means adding one subpackage — no edits to the engine or shared code. See [Adding a Parser](adding-parsers.md).

**SymPy as the intermediate representation.**
All rate expressions, fluxes, and ODEs live as SymPy objects inside `Network`. Code generation (`Codegen`) calls SymPy's language-specific printers (`ccode`, `cxxcode`, `fcode`, etc.), so adding a new target language is a single `Language` subclass in `jaff/codegen/_languages.py`.

**Plugin-based code generation.**
`Builder` discovers plugins at `jaff.plugins.<name>.plugin` and calls their `#!python main()`. Each plugin owns its template files and knows nothing about the parser. This keeps solver-specific logic out of the core library.

**`Catalogue[T]` for all domain collections.**
`Species`, `Reactions`, and `Elements` all inherit from `Catalogue`, giving O(1) lookup by integer index, slice, string name, _and_ serialized canonical name. The serialized form (e.g. `"+/H/H/O"` for H₂O⁺) enables duplicate detection that is independent of input name formatting.

**`.jaff` binary format.**
Networks can be saved as gzip-compressed JSON (`.jaff` files) via `io/_io.py`. On load, SymPy expressions are reconstructed from the versioned compact encoding in `common/_sympy_json.py`. This avoids re-parsing large networks on repeated runs.

## Utility Scripts

`src/jaff/_utils/` holds standalone, easy-to-run scripts for maintaining the bundled data. They are **not** part of the runtime data flow — they are run by hand (or during maintenance) to regenerate the assets in `data/` and `db/jaff.db`.

The cross-section scripts are ordered as a pipeline: download raw NORAD data,
collapse the per-reaction files into combined HDF5 files, then build the
SQLite lookup tables that JAFF queries at runtime.

| Script                          | Purpose                                                                                                                                                                |
| ------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `generate_mass_table.py`        | Read `data/atom_mass.csv` and (re)build the element mass tables inside `db/jaff.db`.                                                                                   |
| `download_nahar_xsecs.py`       | Download NORAD/OP (Nahar, OSU) ground-state photoionisation cross sections (Z = 1..26) into `data/xsecs/op/` using serialized reaction names.                          |
| `collapse_xsecs_hdf5.py`        | Merge the per-reaction Leiden and NORAD files into combined `leiden.hdf5` / `norad.hdf5` (one group per reaction, photon energy in eV, σ in cm²).                      |
| `split_xsecs_photodecay.py`     | Split the source dissociation/ionisation datasets into the single `photodecay` channel used by the collapsed HDF5 files.                                               |
| `generate_photo_xsecs_table.py` | Build the `photo_reaction_cross_sections` table in `db/jaff.db` from the collapsed HDF5 files (`photo_absorption` flag, `decay_type` + `file.hdf5::<group>` pointers). |
| `generate_ion_xsecs_table.py`   | Build the `verner_cross_sections` table in `db/jaff.db` from the Verner (1996) analytic-fit parameters in `data/xsecs/verner_1996.csv`.                                |
| `build_shielding_hdf5.py`       | Collapse the per-species Leiden line-shielding tables into `data/shielding/leiden.hdf5` (one group per reaction).                                                      |
| `compress_hdf5.py`              | Write gzip-compressed copies of local data files (HDF5 datasets ≥ `--min-size`; others byte-copied) under `--outdir` and print/merge their `registry.txt` lines for the mirror. |

Run a script as a module from the project root, e.g.:

=== "python"

    ```bash
    python -m jaff._utils.generate_mass_table
    ```

=== "uv"

    ```bash
    uv run python -m jaff._utils.generate_mass_table
    ```
