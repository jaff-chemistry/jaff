---
tags:
    - User-guide
    - Network
---

# jaffx

`jaffx` is JAFF's command-line tool for quick network inspection and
rate-coefficient export. Think of it as a thin shell around
[`Network`](network.md): every invocation **loads a network and calls one of its
methods**, then exits. It deliberately does **not** run the code-generation
pipeline. Use [`jaffgen`](../code-generation/jaffgen.md) when you want
generated source.

```bash
jaffx <command> <subcommand> --network <file> [options]
```

---

## How it maps to `Network`

There is nothing in `jaffx` that you cannot do from Python — each subcommand is a
one-line wrapper. Holding that mapping in mind is the whole mental model:

| CLI invocation                | Equivalent Python                     |
| ----------------------------- | ------------------------------------- |
| `jaffx get num-species`       | `Network(...).species.count`          |
| `jaffx get num-reactions`     | `Network(...).reactions.count`        |
| `jaffx export hdf5 -f f.hdf5` | `Network(...).to_hdf5("f.hdf5", ...)` |
| `jaffx export txt -f f.txt`   | `Network(...).to_txt("f.txt", ...)`   |
| `jaffx export jaff -f f.jaff` | `Network(...).to_jaff("f.jaff")`      |

As with `to_hdf5`, an HDF5 output path that does not end in `.hdf5` or `.hdf`
has its suffix replaced: `-f f.h5` writes `f.hdf5`.

Every run first loads the network (printing the JAFF banner and the usual load
log), so the load-time validation warnings you'd see from `Network(...)` show up
here too.

---

## Shared network arguments

Every subcommand loads a network, so all of them accept the same loading
options. These mirror the [`Network` constructor](network.md#constructor):

| Argument       | Description                                                                      |
| -------------- | -------------------------------------------------------------------------------- |
| `--network`    | Path to the network file, or a predefined network name (**required**)           |
| `--label`      | Override the network label (defaults to the file stem)                           |
| `--funcfile`   | Path to a `.jfunc` auxiliary file; `true` scans the network dir (default); `false` skips loading |
| `--expand-nuclei` | `--expand-nuclei` / `--no-expand-nuclei` — expand `n_<element>_nuc` element-nucleus density sums |
| `--duplicate-policy` | Resolve duplicate rate coefficients over the same temperature range: `preserve-first` (default), `preserve-last`, or `error`. Overrides the network `jaff.toml` `[network.reactions] duplicate_policy` |
| `--set`        | Set any network config key, `--set key.path=value` (repeatable; the value is TOML, a bare word is a string), e.g. `--set network.reactions.t_cutoff=clip` |

There is no `--strict` option: strict validation is set through `--set`. With
it, a validation problem ends the run with exit status 1:

```bash
jaffx get num-species --network networks/GOW/GOW.jet --set network.strict=true
```

```text
WARNING  Electron recombination not found for H2+
ERROR    Recombination errors found
```

---

## `jaffx get`

Query a scalar property of the network. Output is written through the logger, so
it carries the usual `INFO` prefix.

### `get num-species`

```bash
jaffx get num-species --network networks/GOW/GOW.jet
```

```text
INFO     Total number of species: 18
```

### `get num-reactions`

```bash
jaffx get num-reactions --network networks/GOW/GOW.jet
```

```text
INFO     Total number of reactions: 50
```

---

## `jaffx export`

Export the network in one of three formats. The first two tabulate rate
coefficients against temperature; the third serializes the whole network.

### `export hdf5` and `export txt`

`export hdf5` writes an HDF5 rate table; `export txt` writes the same data as a
whitespace-separated text table. Both take the identical option set and forward
it straight to [`Network.to_hdf5`](network.md#export-and-caching) /
`Network.to_txt`:

| Argument          | Description                                                                  |
| ----------------- | ---------------------------------------------------------------------------- |
| `--file`, `-f`    | Output file path (**required**)                                              |
| `--tmin`          | Minimum tabulation temperature (default: minimum over reactions)             |
| `--tmax`          | Maximum tabulation temperature (default: maximum over reactions)             |
| `--nT`            | Initial number of temperature points (default 64; adaptive refinement may add more) |
| `--err-tol`       | Relative interpolation error tolerance (default 0.01). Adaptive sampling is always applied: an unset or zero value uses 0.01 |
| `--rate-min`      | Small floor added to the relative-error denominator so near-zero rates do not divide by zero (default 1e-30); rates are not clamped to it |
| `--rate-max`      | Rates above this ceiling are clipped to prevent overflow (default 1e100)     |
| `--fast-log` / `--no-fast-log` | Sample equally in `#!python fast_log2(T)` space instead of `#!python log(T)` |
| `--include-all` / `--no-include-all` | Include every reaction, marking non-tabulatable ones `NaN`   |
| `--verbose` / `--no-verbose`, `-v` | Print progress during adaptive refinement                      |

```bash
# HDF5 table, adaptively refined to 0.1% interpolation error
jaffx export hdf5 \
    --network networks/GOW/GOW.jet \
    --file    rates.hdf5           \
    --tmin    10                   \
    --tmax    1e4                  \
    --nT      200                  \
    --err-tol 1e-3

# Plain-text table starting from a 100-point grid (refined to the default 1%)
jaffx export txt \
    --network networks/GOW/GOW.jet \
    --file    rates.txt            \
    --tmin 10 --tmax 1e4 --nT 100
```

Only rates that depend solely on temperature are tabulated; reactions involving
`av`, `crate`, or undefined functions are dropped, and so are rates that are
constant in temperature (under `--include-all` every reaction is kept, the
non-tabulatable ones as `NaN`).

Unlike `Network.to_hdf5(..., err_tol=None)`, the CLI cannot switch adaptive
sampling off: `--err-tol` always falls back to 0.01.

### `export jaff`

Serialize the entire parsed network to a gzip-compressed JSON `.jaff` file —
exactly [`Network.to_jaff`](network.md#jaff-serialized-format). Re-loading it with
`Network("...jaff")` skips parsing and the expensive SymPy assembly, so this is
the recommended cache for large networks.

```bash
jaffx export jaff \
    --network uclchem_small_gas \
    --file    uclchem.jaff
```

<!-- prettier-ignore -->
!!! warning "Same serialization limit as `Network.to_jaff`"
    Networks whose rates contain an undefined function — most often an
    unresolved `photorates(...)` (a photo-reaction loaded without `radiation`)
    or a custom `interp(...)` — cannot be serialized. GOW is one: `jaffx export
    jaff --network networks/GOW/GOW.jet ...` currently stops with a Python
    traceback ending in
    `ValueError: Cannot serialize: expression contains undefined SymPy function(s)`
    rather than a one-line `Error:` message.

---

## Examples

```bash
# Quick species count
jaffx get num-species --network networks/h_photoionization/h_photo.jet

# 500-point HDF5 table with tight adaptive sampling
jaffx export hdf5 \
    --network networks/GOW/GOW.jet \
    --file    GOW_rates.hdf5       \
    --tmin 10 --tmax 1e5 --nT 500  \
    --err-tol 1e-4

# Text table including all reactions (NaN for non-tabulatable)
jaffx export txt \
    --network    networks/GOW/GOW.jet \
    --file       GOW_rates.txt        \
    --include-all

# Cache the network for fast re-loading
jaffx export jaff \
    --network uclchem_small_gas \
    --file    uclchem.jaff
```
