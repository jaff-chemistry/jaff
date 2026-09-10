---
tags:
    - User-guide
    - Code-generation
---

# Configuration File (`jaffgen.toml`)

A `jaffgen.toml` declares a [`jaffgen`](jaffgen.md) run once, so you don't repeat
CLI flags every time. It is loaded when you pass `--config <file>`, **or**
automatically when a file named `jaffgen.toml` turns up among the gathered template
files — which is how a bundled template (like `microphysics`) can ship its own
settings.

---

## Priority order

Every setting is resolved highest-wins, so the config file fills gaps the CLI
leaves and overrides constructor defaults:

1. Explicit CLI argument (e.g. `--network`)
2. `jaffgen.toml` value
3. `Network` constructor default

<!-- prettier-ignore -->
!!! note "Relative paths are resolved from the config file's directory"
    Any path that comes from the `jaffgen.toml` (network, funcfile, input/output
    dirs, table files) is resolved relative to **where the `jaffgen.toml` lives**,
    not the current working directory. Paths passed on the CLI are resolved
    relative to the CWD.

The smallest useful config is a single section — the bundled `microphysics`
template, for instance, ships only a `[network.radiation]` block:

```toml
[network.radiation]
bands = [13.6, "inf"]
power_law_index = 0
energy_density = false
rsl = 2.99792458e10
```

---

## `[jaffgen]` section

Controls the pipeline itself — mirrors the `jaffgen` CLI flags.

```toml
[jaffgen]
output_dir   = "../generated"          # where generated files are written
input_dir    = "."                     # directory of template files
input_files  = ["extra.cpp"]           # individual files (combined with input_dir)
template     = "microphysics"          # built-in template collection name
network      = "networks/GOW/GOW.jet"  # network file or built-in network name
default_lang = "cxx"                   # fallback language for unknown extensions
```

| Key            | Type        | Description                                                         |
| -------------- | ----------- | ------------------------------------------------------------------- |
| `output_dir`   | `str`       | Output directory (created if absent)                                |
| `input_dir`    | `str`       | Directory of template files to process                              |
| `input_files`  | `list[str]` | Individual template files; combined with `input_dir` and `template` |
| `template`     | `str`       | Built-in collection under `jaff/templates/generator/`               |
| `network`      | `str`       | Network file path, or a built-in network name                       |
| `default_lang` | `str`       | Fallback language for unrecognised extensions                       |

---

## `[network]` section

Sets `Network` constructor options.

```toml
[network]
label            = "GOW-2017"
funcfile         = "networks/GOW/GOW.jfunc"
replace_nH       = true
errors           = false
duplicate_policy = "preserve-first"
```

| Key                | Type            | Default            | Description                                                                                                                                                                                                                                       |
| ------------------ | --------------- | ------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `label`            | `str`           | file stem          | Human-readable network name                                                                                                                                                                                                                       |
| `funcfile`         | `str` or `bool` | `true`             | Path to a `.jfunc` auxiliary file; `true` scans the network dir, `false` skips loading                                                                                                                                                            |
| `replace_nH`       | `bool`          | `true`             | Expand `nh` / `nhe` shorthands in rate expressions                                                                                                                                                                                                |
| `errors`           | `bool`          | `false`            | Treat conservation violations as fatal                                                                                                                                                                                                            |
| `duplicate_policy` | `str`           | `"preserve-first"` | Resolve duplicate rate coefficients over the same temperature range: `preserve-first`, `preserve-last`, or `error`. The `--duplicate-policy` CLI flag overrides this; this in turn overrides the network `jaff.toml` `[network].duplicate_policy` |

Besides these scalar keys, `[network]` also holds the subtables below
(`radiation`, `rates`, `reactions`).

---

## `[network.rates]` and `[network.reactions."<serialized>"]` (temperature cutoffs)

The per-reaction and global **temperature-cutoff** behaviour uses the same
`[network.rates]` / `[network.reactions."<srxn>"]` schema documented for
[`jaff.toml`](../working-with-networks/jaff-toml.md). Placing it here lets a
particular `jaffgen` run override the network's own `jaff.toml` defaults:

```toml
[network.rates]
T_cutoff = "extrapolate"                   # global default for this run

[network.reactions."CO._PHOTON__C.O"]
T_cutoff = "clip"                          # per-reaction override
```

A per-reaction value wins over any global; within a scope, `jaffgen.toml`
overrides `jaff.toml`. See
[Network Configuration → Resolution order](../working-with-networks/jaff-toml.md#resolution-order).

---

## `[network.radiation]` section

Configures the photochemistry radiation field. Present this block to enable
photochemistry radiation ode and jacobian radiation generation terms; omit it
(or give an empty `bands`) to leave it off.

```toml
[network.radiation]
bands            = [13.6, "inf"]    # band edges in eV; "inf" for an open upper bound
power_law_index  = 0                # spectral power-law index
energy_density   = false            # true = energy density; false = photon density
rsl              = 2.99792458e10    # speed of light (cm/s). Used to configure reduced speed of light for solvers
background_field = "draine"         # reference field used to scale chi_pe
```

| Key                | Type             | Default                 | Description                                                                                       |
| ------------------ | ---------------- | ----------------------- | ------------------------------------------------------------------------------------------------- |
| `bands`            | `list`           | `[]`                    | Band boundaries in eV; omit to disable photochemistry                                             |
| `power_law_index`  | `int or float`   | `0`                     | Spectral index for band integration                                                               |
| `energy_density`   | `bool`           | `false`                 | Radiation density variable type. `radeden` when `true` else `photden`                             |
| `rsl`              | `float` or `str` | `constants.c.cgs.value` | Speed of light override (maps to the `c` constructor arg). Becomes a symbol if passed as a string |
| `background_field` | `str`            | `"draine"`              | Reference radiation field (HDF5 group name) used to scale the photoelectric-band `chi_pe` symbol  |

`power_law_index` is used to configure the weight factor of the photo-reaction cross-sections (Refer to the [Photochemistry](../designing-networks/photochemistry.md) section for more information).

`background_field` only matters when the [dust module](#networkdust-section) is
enabled; it names the reference field that `chi_pe` is scaled against.

---

## `[network.dust]` section

Present this (possibly empty) table to enable the **dust module**. It maps to the
`dust=True` constructor argument and activates dust-driven physics — currently
photoelectric emission, which supplies the [`chi_pe`](../designing-networks/photochemistry.md#self-consistent-photoelectric-field-chi_pe)
symbol (the local field scaled to the photoelectric band).

```toml
[network.dust]
# presence alone enables the dust module; no keys are required yet
```

The dust module needs radiation enabled — a `[network.radiation]` block with
non-empty `bands` — because `chi_pe` is built from the radiation bands and the
`background_field` reference; generation aborts otherwise.

---

## `[network.reactions.<serialized>.shielding]` section

Attaches a shielding factor to one photo-reaction, keyed by the reaction's
[serialized form](../working-with-networks/reactions.md). The factor multiplies
that reaction's rate coefficient at runtime; the conceptual model, the formulae
and the original papers are in the
[Shielding](../designing-networks/photochemistry.md#shielding) section. The
reaction **must** be a photo-reaction or generation aborts.

The serialized key contains `.` separators, so it **must be quoted** in the
table header — otherwise TOML reads the dots as nested tables. Photo-reactions
also carry the `_PHOTON` agent in their serialized form.

```toml
# Leiden tabulated line shielding
[network.reactions."CO._PHOTON__C.O".shielding]
type        = "leiden"           # default if omitted
radiation   = "ISRF"
shielded_by = ["self", "H2"]

# H2 self-shielding (Hartwig et al. 2015)
[network.reactions."H2._PHOTON__H.H".shielding]
type      = "hg2015"
min_ncol  = 1.0e-35
min_vdisp = 1.0e-20
```

Common key:

| Key    | Type  | Default    | Description                                                                 |
| ------ | ----- | ---------- | --------------------------------------------------------------------------- |
| `type` | `str` | `"leiden"` | Shielding function: `"leiden"`, `"db1996"`, or `"hg2015"`. Case-insensitive |

`type = "leiden"` keys:

| Key           | Type   | Default  | Description                                                                                                    |
| ------------- | ------ | -------- | -------------------------------------------------------------------------------------------------------------- |
| `shielded_by` | `list` | required | Shielding species; allowed: `"self"`, `"H2"`, `"H"`, `"C"`, `"N2"`, `"CO"`. Per-species factors are multiplied |
| `radiation`   | `str`  | `"ISRF"` | Radiation-field subgroup in the Leiden table                                                                   |

`type = "db1996"` / `"hg2015"` keys (only on the `H2._PHOTON__H.H` reaction):

| Key         | Type    | Default | Description                          |
| ----------- | ------- | ------- | ------------------------------------ |
| `min_ncol`  | `float` | `1e-50` | Lower floor used in the fit (cm⁻²)   |
| `min_vdisp` | `float` | `1e-50` | Lower floor used in the fit (cm s⁻¹) |

## `[[table]]` section

A `[[table]]` array entry converts a data table from one format to another as
part of the generation run — typically to ship the lookup table that the
generated [interpolation functions](table-interpolation.md) read at runtime. One
block describes one conversion, with a `[table.source]` and a `[table.target]`.
Supported directions: **HDF5 → HDF5**, **CSV → HDF5**, and **CSV → CSV**.

### How the conversion works

The engine loads the source into a flat tree, builds a target tree from your
`[table.target]` headings, then writes it out:

1. **Source tree.** The source is flattened to a lookup keyed by absolute path.
   An HDF5 source becomes `{ "/co/TCO": <dataset>, "/co/L0CO": <dataset>, … }`; a
   CSV source becomes `{ "T0": <column>, "NeffCO": <column>, … }`.
   `path = "default"` is shorthand for the network's own rate table,
   `<network_dir>/<network_stem>.hdf5`.

2. **Target headings are output paths.** Every `[table.target]` key beginning
   with `/` is a path in the **output** HDF5 file. What you place under that
   heading says where its data comes from:
    - **`h5path = "/old/path"`** (HDF5 → HDF5) — move the source dataset/group at
      `/old/path` to this heading's path. The whole source tree is copied first,
      so datasets you don't remap pass through unchanged; a remapped source path
      is removed from its old location. Omit `h5path` to leave a dataset where it
      already is.
    - **`h5path = ["/a", "/b", …]`** (HDF5 → HDF5, *composite*) — fold several
      source datasets into a **single** dataset at this heading's path. The
      folded-in sources are consumed (removed from their old locations). A `type`
      key picks the layout:
        - **`type = "compound"`** (default) — a record/table dataset, one named
          column per source. Column names come from a `names = [...]` list, or
          from each source path's final component when `names` is omitted. All
          sources must share a length.
        - **`type = "ndarray"`** — a single array of shape `(Ncols, *xlens)`:
          `Ncols` source columns stacked along a new leading axis, each reshaped
          to the grid whose axis lengths are the lengths of the datasets listed
          under this heading's `regrid.x`. Each source's flat size must equal the
          product of those axis lengths.
    - **`col = "ColName"`** (CSV → HDF5) — write the named CSV column as the
      dataset at this heading's path. Only columns named by a `col` are written;
      nested paths create the intermediate groups (e.g. `/co/1d/Temp`).

3. **Attributes.** A `.attrs` sub-table on a heading attaches HDF5 attributes to
   that path. Each attribute **value** is one of:
    - a `"/target/path.property"` **reference** — a statistic computed from the
      **target** tree at write time. Supported properties: `max`, `min`, `mean`,
      `median`, `length`. Because they read the target tree, the referenced path
      must exist in the output (e.g. the remapped path, not the original source
      path).
    - a plain **literal** — number, boolean, or string (e.g. `units = "K"`,
      `Ndim = 2`, `spacing = "log"`).
    - a **list** mixing either of the above, resolved element-wise — e.g.
      `Nx = ["/co/TCO.length", "/co/NeffCO.length"]` (a list of references) or
      `spacing = ["log", "log"]` (a list of literals).

<!-- prettier-ignore -->
!!! note "Distinguishing references from literals"
    A string is treated as a computed reference only when it ends in
    `.<property>` for one of the supported property names; every other string
    (including `"log"` and `"K"`) is stored verbatim.

`default_group` sets the output root group (default `/`). For CSV sides,
`delimiter` and `comment` configure parsing/writing.

### HDF5 → HDF5

Copy the source tree into a new file, remapping selected paths and attaching
computed attributes. Here `/co/TCO` is republished as `/temperature`:

```toml
[[table]]
[table.source]
path = "default"             # the network's own HDF5 rate table

[table.target]
path          = "GOW.hdf5"
default_group = "/"

[table.target."/temperature"]
h5path = "/co/TCO"           # move source /co/TCO here (other datasets copy through)

[table.target."/temperature".attrs]
tmax = "/temperature.max"    # computed from the data now at /temperature
npts = "/temperature.length"
```

#### Composite datasets

Fold several source datasets into one. As a compound (named-column) dataset:

```toml
[table.target."/c0/data"]
h5path = ["/co/LLTECO", "/co/alphaCO", "/co/nhalfCO"]
type   = "compound"                       # default; may be omitted
names  = ["LLTE", "alpha", "nhalf"]       # optional; else path stems are used
```

Or as a single `(Ncols, *xlens)` N-D array, gridded over the `regrid.x` axes:

```toml
[table.target."/c0/data"]
h5path = ["/co/LLTECO", "/co/alphaCO"]
type   = "ndarray"

[table.target."/c0/data".regrid]
x = ["/co/TCO", "/co/NeffCO"]             # shape becomes (2, len TCO, len NeffCO)
```

### CSV → HDF5

Write named CSV columns as HDF5 datasets — the usual way to turn a
`co_1d.csv`-style table into the HDF5 file an interpolation routine reads:

```toml
[[table]]
[table.source]
delimiter = " "
comment   = "#"
path      = "networks/GOW/co_1d.csv"

[table.target]
path          = "GOW.hdf5"
default_group = "/"

[table.target."/co/1d/Temp"]
col = "T0"                   # CSV column "T0" → dataset /co/1d/Temp

[table.target."/co/1d/Temp".attrs]
max       = "/co/1d/Temp.max"
min       = "/co/1d/Temp.min"
t0_length = "/co/1d/Temp.length"
```

### CSV → CSV

Select specific columns from one CSV and rewrite them, optionally changing the
delimiter:

```toml
[[table]]
[table.source]
delimiter = " "
comment   = "#"
path      = "networks/GOW/co_1d.csv"
cols      = ["T0", "NeffCO"]     # only these columns are kept

[table.target]
delimiter = ","
path      = "GOW.csv"
```

---

## Dummy example

The reference configuration below exercises every section.

```toml
[jaffgen]
output_dir   = "../generated"
input_dir    = "."
input_files  = ["../new.cpp", "test.cpp"]
template     = "microphysics"
network      = "networks/GOW/GOW.jet"
default_lang = "cxx"

[network]
label      = "GOW-generator"
funcfile   = "networks/GOW/GOW.jfunc"
replace_nH = true
errors     = false

[network.radiation]
bands           = [13.6, "inf"]
power_law_index = 0
energy_density  = false
rsl             = 2.99792458e10

# HDF5 → HDF5
[[table]]
[table.source]
path = "default"

[table.target]
path          = "GOW.hdf5"
default_group = "/"

[table.target."/temperature"]
h5path = "/co/TCO"

[table.target."/temperature".attrs]
tmax = "/temperature.max"
tmin = "/temperature.min"

# CSV → HDF5
[[table]]
[table.source]
delimiter = " "
comment   = "#"
path      = "networks/GOW/co_1d.csv"

[table.target]
path          = "GOW.hdf5"
default_group = "/"

[table.target."/co/1d/Temp"]
col = "T0"

# CSV → CSV
[[table]]
[table.source]
delimiter = " "
comment   = "#"
path      = "networks/GOW/co_1d.csv"
cols      = ["T0", "NeffCO"]

[table.target]
delimiter = " "
comment   = "#"
path      = "GOW.csv"
```
