---
tags:
    - User-guide
    - Code-generation
---

# jaffgen

`jaffgen` is JAFF's command-line entry point for the **full code-generation
pipeline**. It gathers a set of [template files](template-syntax.md), loads a
chemical reaction network, runs the
[template engine](template-syntax.md) over each file, and writes the generated
source into an output directory.

```bash
jaffgen --network <file> [options]
```

---

## How it works

A `jaffgen` run is five phases:

1. **Resolve configuration** — merge a
   [`jaffgen.toml`](../configuration/index.md) (or the one `--template`
   bundles) and the command line (see [below](#configuration-priority)), and
   check the `[jaffgen]` section.
2. **Load the network** — exactly one `Network`: it merges the run's
   `[network]` settings over its own `jaff.toml` and validates them.
3. **Gather template files** — from `--indir`, `--files`, and `--template`.
   These **combine** — a single run can pull from all three.
4. **Generate** — run [`TemplateParser`](template-syntax.md) on every gathered
   file and write the result to a new [run folder](#run-folders) in
   `--outdir`, preserving each file's path relative to its source directory,
   together with any data tables
   ([`[[jaffgen.table]]`](../configuration/data-tables.md), shielding tables).
5. **Record** — write the run's whole config to `config.toml` in the run folder.

The mental model: _jaffgen is the template engine driven over a file set with
one loaded network behind it._

---

## Input sources

The three input flags are additive, not exclusive:

| Flag         | Adds                                                          |
| ------------ | ------------------------------------------------------------- |
| `--indir`    | Every file in a directory (recursive), keeping its structure  |
| `--files`    | Specific individual files                                     |
| `--template` | A built-in collection from `jaff/templates/generator/<name>/` |

For `--template`, every file in the template's folder is rendered, keeping its
path relative to the folder; files without `$$JAFF` directives (a solver's
driver, its Makefile, bundled library sources) come out unchanged.

---

## Arguments

### Network and network options

| Argument       | Description                                                                        |
| -------------- | ---------------------------------------------------------------------------------- |
| `--network`    | Network file path **or** a built-in network name; required (CLI or `jaffgen.toml`) |
| `--network-config` | Path to another [`jaff.toml`](../configuration/index.md), merged over `<network_dir>/jaff.toml` |
| `--label`      | Override the network label (defaults to the file stem)                             |
| `--funcfile`   | Path to a `.jfunc` auxiliary file; `true` scans the network dir; `false` skips     |
| `--expand-nuclei` | `--expand-nuclei` / `--no-expand-nuclei` — expand `n_<element>_nuc` element-nucleus sums (default: on) |
| `--strict`     | `--strict` / `--no-strict` — exit instead of warning when network validation finds a problem: unconserved mass or charge, sink/source species, missing recombinations, or isomers |
| `--duplicate-policy` | Resolve duplicate rate coefficients over the same temperature range: `preserve-first` (default), `preserve-last`, or `error` |

When `--network` is a built-in network name (a sub-directory of `networks/`),
it resolves to the single `.jet` file inside that directory — the directory
must contain exactly one. A built-in name takes precedence over a same-named
path on disk. jaffgen forwards the resolved absolute path to `Network`, which
performs this resolution.

### Output

| Argument   | Description                                                           |
| ---------- | --------------------------------------------------------------------- |
| `--outdir` | Output directory; each run writes into a new [run folder](#run-folders) in it (created if absent; defaults to `<repo>/generated/`) |
| `--echo`   | `--echo` / `--no-echo` — copy `$$JAFF` directives into the output as comments (default: on; overrides `echo` in `jaffgen.toml`) |

### Input sources

| Argument     | Description                                             |
| ------------ | ------------------------------------------------------- |
| `--indir`    | Directory of template files                             |
| `--files`    | One or more individual template files (comma-separated: `a.txt,b.txt`) |
| `--template` | Built-in template collection name                       |

### Configuration file

| Argument   | Description                                                                                                  |
| ---------- | ------------------------------------------------------------------------------------------------------------ |
| `--config` | Path to a [`jaffgen.toml`](../configuration/index.md); without it, a `jaffgen.toml` bundled with `--template` is used |
| `--set`    | `--set key.path=value` (repeatable) — set any [config key](../configuration/reference.md), e.g. `--set network.reactions.t_cutoff=clip`; wins over the other options |

---

## Configuration priority

Each setting is resolved highest-wins:

1. `--set key.path=value` items
2. Other command-line options (e.g. `--network`, `--strict`)
3. [`jaffgen.toml`](../configuration/index.md) value (from `--config`, or bundled with the template)
4. `--network-config`, if given
5. The network's `jaff.toml`
6. The config defaults

Paths taken from a config file are resolved **relative to that file's
directory**; paths given on the CLI are resolved relative to the current
directory. Each option sets one config key; see the
[option table](../configuration/index.md#from-the-command-line).

The run folder also receives `config.toml`: the run's whole config after
merging, every key and default included — an exact record of what the run used,
loadable again with `--config`. It records the output directory actually used,
the default one included, not the run folder, so a rerun from it writes a new
run folder next to this one. (The `jaffgen.toml` itself is not copied.)

### Run folders

Each run writes into a new folder in the output directory, named
`<network>_<template>_<YYYYmmdd-HHMMSS>`, so earlier runs are never
overwritten:

```text
output/
├── GOW_microphysics_20261031-101500/
└── GOW_microphysics_20261031-103012/
```

- `<network>` is the network's [`label`](../configuration/index.md) (by
  default the network file's name without its extension).
- `<template>` is the `--template` name; without one, the `--indir` folder's
  name; with only `--files`, `custom`.
- Characters other than letters, digits, `_`, `.` and `-` become `_`.
- A second run in the same second gets `-2` appended (then `-3`, ...).

The last log line names the run folder ("Generated files can be found at ...").

---

## Examples

### Generate from a template directory

```bash
jaffgen \
    --network networks/GOW/GOW.jet \
    --indir   templates/           \
    --outdir  output/
```

For predefined networks, the folder name can also be passed as the network

```bash
jaffgen \
    --network GOW                  \
    --indir   templates/           \
    --outdir  output/
```

### Process specific files

The language of each file comes from the file itself: a `$$JAFF LANG <language>`
line, or else its extension (see
[Template syntax](template-syntax.md)). Files named here with an extension JAFF
knows, such as `.cpp`, need nothing extra:

```bash
jaffgen \
    --network networks/GOW/GOW.jet \
    --files   rates.cpp,odes.cpp   \
    --outdir  output/
```

### Use a built-in template collection

```bash
jaffgen \
    --network  networks/GOW/GOW.jet \
    --template microphysics         \
    --outdir   output/
```

<!-- prettier-ignore -->
!!! note "Templates that bundle build files"
    The `microphysics` collection ships non-source files (e.g. `Make.package`
    and an extension-less `_parameters`) alongside its `.cpp`/`.H` templates.
    Files with no `$$JAFF` directives and no known extension (a `Makefile`,
    `CMakeLists.txt`, a `.toml`) are copied unchanged. A file that does contain
    directives names its language with a `$$JAFF LANG` line (and, if needed, a
    `$$JAFF COMMENT` token for echoed directives), as `_parameters` does.

### Combine a template collection with custom files

```bash
jaffgen \
    --network  networks/GOW/GOW.jet \
    --template microphysics         \
    --files    custom_rhs.cpp       \
    --outdir   output/
```

### Drive everything from a config file

A `jaffgen.toml` can carry the network, inputs, output, and radiation settings, so
the command collapses to:

```bash
jaffgen --config jaffgen.toml
```

A bundled template can even supply its own config. The `microphysics` collection
ships a `jaffgen.toml` with a `[network.radiation]` block, so `--template microphysics`
uses it and turns on photochemistry without any extra flags. See the
[configuration reference](../configuration/reference.md) for every key.

---

## Built-in templates

Collections live under `src/jaff/templates/generator/`. List them with:

```bash
ls src/jaff/templates/generator/
```
