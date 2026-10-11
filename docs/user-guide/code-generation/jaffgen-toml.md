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
profile_index = 0
mode = "nph"
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
expand_nuclei    = true
errors           = false
duplicate_policy = "preserve-first"
```

| Key                | Type            | Default            | Description                                                                                                                                                                                                                                       |
| ------------------ | --------------- | ------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `label`            | `str`           | file stem          | Human-readable network name                                                                                                                                                                                                                       |
| `funcfile`         | `str` or `bool` | `true`             | Path to a `.jfunc` auxiliary file; `true` scans the network dir, `false` skips loading                                                                                                                                                            |
| `expand_nuclei`    | `bool`          | `true`             | Expand `n_<element>_nuc` element-nucleus shorthands in rate expressions                                                                                                                                                                           |
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

### Per-reaction `pi_database`

Photoionization reactions may override the global database:

```toml
[network.reactions."C._PHOTON__C+.e-"]
pi_database = "verner"
```

The key is only valid on photoionization reactions (any other reaction raises a
`ParserError`); a value from `jaffgen.toml` wins over `jaff.toml`. The override
is keyed by the reaction as written in the network (its `serialized` key), even
when `use_proxy_photoreaction` maps it to a different database key. Per-reaction
overrides from `jaffgen.toml` are not applied when the network is loaded from a
`.jaff` file; the value stored in the `.jaff` is used.
See [`pi_database`](#networkradiation-section) for the global choice.

---

## `[network.radiation]` section

Configures the photochemistry radiation field. Present this block to enable
photochemistry radiation ode and jacobian radiation generation terms; omit it
(or give an empty `bands`) to leave it off.

```toml
[network.radiation]
bands            = [13.6, "inf"]    # band edges in eV; "inf" for an open upper bound
profile_index  = 0                # spectral power-law index; or one per band, e.g. [0, 1]
mode   = "nph"            # "nph" = photon number density; "u" = energy density
rsl              = 2.99792458e10    # speed of light (cm/s). Used to configure reduced speed of light for solvers
background_field = "draine"         # reference field used to scale chi_pe
use_proxy_photoreaction = false     # use proxy photo-reactions when computing cross-sections
pi_database = "norad"               # photoionization xsecs: norad | verner | leiden
```

| Key                       | Type                     | Default                 | Description                                                                                                      |
| ------------------------- | ------------------------ | ----------------------- | ---------------------------------------------------------------------------------------------------------------- |
| `bands`                   | `list`                   | `[]`                    | Band boundaries in eV; omit to disable photochemistry                                                            |
| `profile_index`           | `int`, `float` or `list` | `0`                     | Spectral index for band integration; a list gives one index per band (length `len(bands) - 1`)                   |
| `mode`                    | `str`                    | `"nph"`                 | Radiation density variable type: `"nph"` (photon number density, `photden`) or `"u"` (energy density, `radeden`) |
| `rsl`                     | `float` or `str`         | `constants.c.cgs.value` | Speed of light override (maps to the `c` `RadiationProps` arg). Becomes a symbol if passed as a string           |
| `background_field`        | `str`                    | `"draine"`              | Reference radiation field (HDF5 group name) used to scale the photoelectric-band `chi_pe` symbol                 |
| `use_proxy_photoreaction` | `bool`                   | `false`                 | Whether to use proxy photo-reactions when computing cross-sections instead of bypassing them                     |
| `pi_database`             | `str`                    | `"norad"`               | Photoionization cross-section database: `"norad"`, `"verner"` or `"leiden"` (case-insensitive)                   |

`profile_index` is used to configure the weight factor of the photo-reaction cross-sections (Refer to the [Photochemistry](../designing-networks/photochemistry.md) section for more information). A scalar applies the same index to every band; a list such as `profile_index = [0, 1]` sets one index per band and must have exactly `len(bands) - 1` entries.

`pi_database` selects the photoionization cross-section database:
`norad` (default; NORAD/Nahar R-matrix ground-state cross-sections including
resonances, with theoretical thresholds), `verner` (Verner et al. 1996 analytic
fits, integrated symbolically), or `leiden` (Heays et al. 2017). Values are
case-insensitive. It affects photoionization only; photodissociation and
photoabsorption always use Leiden. When the chosen database lacks a reaction,
jaff falls back to `norad`, then `verner`, then `leiden`; after the network
finishes loading, a single summary warning lists every fallback
(`requested -> used: key, ...`). If no database has the reaction, radiation is
enabled and the reaction has no custom rate, loading fails with a `ParserError`;
otherwise the reaction gets no cross-section and keeps its own rate. A single
reaction can override this choice, see
[per-reaction `pi_database`](#per-reaction-pi_database).

`background_field` only matters when the [dust module](#networkdust-section) is
enabled; it names the reference field that `chi_pe` is scaled against.

---

## `[network.dust]` section

Present this table to enable the **dust module**. Its presence maps to a
`dust_props=DustProps(...)` constructor argument and activates dust-driven
physics — currently photoelectric emission, which supplies the
[`chi_pe`](../designing-networks/photochemistry.md#self-consistent-photoelectric-field-chi_pe)
symbol (the local field scaled to the photoelectric band).

```toml
[network.dust]
rv               = 3.1              # extinction curve Rv: one of 3.1, 4.0, 5.5
u_reduction      = "extinction"     # radiation energy-density reduction kind
f_reduction      = "extinction"     # radiation flux reduction kind
pe_threshold_low = 6                # photoelectric band lower edge (eV)
pe_threshold_high = 13.6            # photoelectric band upper edge (eV)
```

| Key                 | Type    | Default        | Description                                                                                            |
| ------------------- | ------- | -------------- | ------------------------------------------------------------------------------------------------------ |
| `rv`                | `float` | `3.1`          | Extinction-curve total-to-selective ratio; one of `3.1`, `4.0`, `5.5`                                  |
| `u_reduction`       | `str`   | `"absorption"` | Radiation energy-density reduction kind: `extinction`, `absorption`, `scattering`, `transport`, `none` |
| `f_reduction`       | `str`   | `"transport"`  | Radiation flux reduction kind (same set of values as `u_reduction`)                                    |
| `pe_threshold_low`  | `float` | `6.0`          | Photoelectric band lower edge in eV (grain work function)                                              |
| `pe_threshold_high` | `float` | `13.6`         | Photoelectric band upper edge in eV (hydrogen ionisation edge)                                         |

The dust module needs radiation enabled — a `[network.radiation]` block with
non-empty `bands` — because `chi_pe` is built from the radiation bands and the
`background_field` reference; generation aborts otherwise.

`background_field` only matters when the [dust module](#networkdust-section) is
enabled; it names the reference field that `chi_pe` is scaled against.

---

## `[network.dust]` section

Present this (possibly empty) table to enable the **dust module**. It maps to the
`dust_props=DustProps(...)` constructor argument and activates dust-driven physics — currently
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

## `[network.eos]` section

Selects the equation of state used for the internal-energy equation and the
Jacobian temperature column. The table maps to the
`eos_props=EosProps(**table)` constructor argument; `type` picks the EOS and
the remaining keys are that type's parameters. Without this table an ideal gas
with `gamma = 1.6666666666667` is used.

```toml
[network.eos]
type  = "ideal"
gamma = 1.4
```

| `type`        | Keys                                                                         | Default                   |
| ------------- | ---------------------------------------------------------------------------- | ------------------------- |
| `ideal`       | `gamma` (float > 1)                                                          | `gamma = 1.6666666666667` |
| `multi_gamma` | `default_gamma` (float > 1), `gamma_map` (table of species name → float > 1) | —                         |

See [eos](../../api/core/network/thermodynamics.md#eosprops) for the full list of types.

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
    - **`h5path = ["/a", "/b", …]`** (HDF5 → HDF5, _composite_) — fold several
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
expand_nuclei = true
errors     = false

[network.radiation]
bands           = [13.6, "inf"]
profile_index = 0
mode  = "nph"
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
