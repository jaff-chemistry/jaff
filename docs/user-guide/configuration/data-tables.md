---
tags:
    - User-guide
    - Configuration
---

# Data tables (`[[jaffgen.table]]`)

A `[[jaffgen.table]]` array entry converts a data table from one format to another as
part of the generation run — typically to ship the lookup table that the
generated [interpolation functions](../code-generation/table-interpolation.md) read at runtime. One
block describes one conversion, with a `[jaffgen.table.source]` and a `[jaffgen.table.target]`.
Supported directions: **HDF5 → HDF5**, **CSV → HDF5**, and **CSV → CSV**.

## How the conversion works

The engine loads the source into a flat tree, builds a target tree from your
`[jaffgen.table.target]` headings, then writes it out:

1. **Source tree.** The source is flattened to a lookup keyed by absolute path.
   An HDF5 source becomes `{ "/co/TCO": <dataset>, "/co/L0CO": <dataset>, … }`; a
   CSV source becomes `{ "T0": <column>, "NeffCO": <column>, … }`.
   `path = "default"` is shorthand for the network's own rate table,
   `<network_dir>/<network_stem>.hdf5`.

2. **Target headings are output paths.** Every `[jaffgen.table.target]` key beginning
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

## HDF5 → HDF5

Copy the source tree into a new file, remapping selected paths and attaching
computed attributes. Here `/co/TCO` is republished as `/temperature`:

```toml
[[jaffgen.table]]
[jaffgen.table.source]
path = "default"             # the network's own HDF5 rate table

[jaffgen.table.target]
path          = "GOW.hdf5"
default_group = "/"

[jaffgen.table.target."/temperature"]
h5path = "/co/TCO"           # move source /co/TCO here (other datasets copy through)

[jaffgen.table.target."/temperature".attrs]
tmax = "/temperature.max"    # computed from the data now at /temperature
npts = "/temperature.length"
```

### Composite datasets

Fold several source datasets into one. As a compound (named-column) dataset:

```toml
[jaffgen.table.target."/c0/data"]
h5path = ["/co/LLTECO", "/co/alphaCO", "/co/nhalfCO"]
type   = "compound"                       # default; may be omitted
names  = ["LLTE", "alpha", "nhalf"]       # optional; else path stems are used
```

Or as a single `(Ncols, *xlens)` N-D array, gridded over the `regrid.x` axes:

```toml
[jaffgen.table.target."/c0/data"]
h5path = ["/co/LLTECO", "/co/alphaCO"]
type   = "ndarray"

[jaffgen.table.target."/c0/data".regrid]
x = ["/co/TCO", "/co/NeffCO"]             # shape becomes (2, len TCO, len NeffCO)
```

## CSV → HDF5

Write named CSV columns as HDF5 datasets — the usual way to turn a
`co_1d.csv`-style table into the HDF5 file an interpolation routine reads:

```toml
[[jaffgen.table]]
[jaffgen.table.source]
delimiter = " "
comment   = "#"
path      = "networks/GOW/co_1d.csv"

[jaffgen.table.target]
path          = "GOW.hdf5"
default_group = "/"

[jaffgen.table.target."/co/1d/Temp"]
col = "T0"                   # CSV column "T0" → dataset /co/1d/Temp

[jaffgen.table.target."/co/1d/Temp".attrs]
max       = "/co/1d/Temp.max"
min       = "/co/1d/Temp.min"
t0_length = "/co/1d/Temp.length"
```

## CSV → CSV

Select specific columns from one CSV and rewrite them, optionally changing the
delimiter:

```toml
[[jaffgen.table]]
[jaffgen.table.source]
delimiter = " "
comment   = "#"
path      = "networks/GOW/co_1d.csv"
cols      = ["T0", "NeffCO"]     # only these columns are kept

[jaffgen.table.target]
delimiter = ","
path      = "GOW.csv"
```

<!-- prettier-ignore -->
!!! note "Relative paths"
    A table's `source.path` and `target.path` are resolved against the
    directory of the `jaffgen.toml` that declares it; the target lands in the
    run's output directory.  Table entries are passed to the converter as
    written: the [configuration schema](reference.md) does not check their keys.
