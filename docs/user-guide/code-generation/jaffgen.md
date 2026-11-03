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

A `jaffgen` run is four phases:

1. **Resolve configuration** — merge three sources by priority (see
   [below](#configuration-priority)): explicit CLI flags, a
   [`jaffgen.toml`](jaffgen-toml.md) file, then `Network` constructor defaults.
2. **Gather template files** — from `--indir`, `--files`, and `--template`.
   These **combine** — a single run can pull from all three.
3. **Load the network** — exactly one `Network`, built from the resolved
   network/label/funcfile/radiation settings.
4. **Generate** — run [`TemplateParser`](template-syntax.md) on every gathered
   file and write the result to `--outdir`, preserving each file's path
   relative to its source directory.

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

For `--template`, generator files are collected first; any **preprocessor** file
(`jaff/templates/preprocessor/<name>/`) whose relative path does not clash is
appended — so the generator always wins on a path collision.

---

## Arguments

### Network and network options

| Argument       | Description                                                                        |
| -------------- | ---------------------------------------------------------------------------------- |
| `--network`    | Network file path **or** a built-in network name; required (CLI or `jaffgen.toml`) |
| `--network-config` | Path to a network [`jaff.toml`](../working-with-networks/jaff-toml.md) (temperature cutoffs); defaults to `<network_dir>/jaff.toml` |
| `--label`      | Override the network label (defaults to the file stem)                             |
| `--funcfile`   | Path to a `.jfunc` auxiliary file; `true` scans the network dir; `false` skips     |
| `--expand-nuclei` | `--expand-nuclei` / `--no-expand-nuclei` — expand `n_<element>_nuc` element-nucleus sums (default: on) |
| `--errors`     | `--errors` / `--no-errors` — exit on conservation violations instead of warning    |
| `--duplicate-policy` | Resolve duplicate rate coefficients over the same temperature range: `preserve-first` (default), `preserve-last`, or `error`. Overrides `jaffgen.toml` and the network `jaff.toml` |

When `--network` is a built-in network name (a sub-directory of `networks/`),
it resolves to the single `.jet` file inside that directory — the directory
must contain exactly one. A built-in name takes precedence over a same-named
path on disk. jaffgen forwards the resolved absolute path to `Network`, which
performs this resolution.

### Output

| Argument   | Description                                                           |
| ---------- | --------------------------------------------------------------------- |
| `--outdir` | Output directory (created if absent; defaults to `<repo>/generated/`) |
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
| `--config` | Path to a [`jaffgen.toml`](jaffgen-toml.md); also auto-detected if a `jaffgen.toml` is among the gathered files |

---

## Configuration priority

Each setting is resolved highest-wins:

1. Explicit CLI flag (e.g. `--network`)
2. [`jaffgen.toml`](jaffgen-toml.md) value (from `--config`, or auto-detected)
3. `Network` constructor default

Paths taken from a `jaffgen.toml` are resolved **relative to the config file's
directory**; paths given on the CLI are resolved relative to the current
directory.

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
auto-detects it and turns on photochemistry without any extra flags. See
[`jaffgen.toml`](jaffgen-toml.md) for the full schema.

---

## Built-in templates

Collections live under `src/jaff/templates/generator/`. List them with:

```bash
ls src/jaff/templates/generator/
```

### `honghu`: a network for Honghu.jl

The `honghu` collection writes one Julia file, `network.jl`, that
[Honghu.jl](https://github.com/chongchonghe/Honghu.jl), a radiation-hydrodynamics code
for prototyping Quokka's algorithms, includes as a chemical network. The file is a
module named `JAFF_<label>` holding one type, `<label> <: Honghu.AbstractNetwork`, so
the network's label must be a valid Julia name. The type gets the species, their
masses, the reactions and the band edges, and `rhses` and `jacobian` with the thermal
row in volumetric `dE/dt` and the radiation terms on. Honghu supplies the temperature,
`tgas = E / (3/2 k_B sum_s n_s)`, which is the ideal gas JAFF writes the Jacobian with.

The template needs radiation bands, `mode = "nph"` and `rsl = "c_hat"`: Honghu stores
each band's photon density divided by `c_hat / c`, so the photo rates read `c_hat`, as
in Quokka. The shipped `jaffgen.toml` sets them; a `--config` of your own replaces it,
so repeat them there:

```toml
[network]
label = "HydrogenIliev2"

[network.radiation]
bands = [13.6, "inf"]
profile_index = 0
mode = "nph"
rsl = "c_hat"
```

```bash
jaffgen --config jaffgen.toml --network hydrogen_iliev2.jet --template honghu --outdir out/
```
