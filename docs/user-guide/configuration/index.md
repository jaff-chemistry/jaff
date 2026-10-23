---
tags:
    - User-guide
    - Configuration
---

# Configuration

JAFF reads its settings from two TOML files that share **one schema**:

| File | Holds | Read by |
| --- | --- | --- |
| `jaff.toml` | `[network]`: how a network is built (reactions, radiation, dust, thermodynamics, …) | [`Network`](../working-with-networks/network.md), [`jaffx`](../working-with-networks/jaffx.md), [`jaffgen`](../code-generation/jaffgen.md) |
| `jaffgen.toml` | `[jaffgen]`: a code-generation run, plus any `[network]` settings for it | [`jaffgen`](../code-generation/jaffgen.md) |

A `jaff.toml` lives next to its network file and describes the network itself;
a `jaffgen.toml` describes one code-generation run. Both accept every
`[network]` key, so a run can adjust the network it generates code for. Every
key is listed in the [reference](reference.md).

```toml
# networks/GOW/jaff.toml
[network.reactions]
t_cutoff = "extrapolate"
```

```toml
# run/jaffgen.toml
[jaffgen]
template     = "microphysics"
network_file = "../networks/GOW/GOW.jet"
output_dir   = "out"

[network.radiation]
bands = [13.6, "inf"]
c     = "c_hat"
```

---

## Layers and precedence

Every source of settings is a **layer** shaped like the file; JAFF merges the
layers, a later one winning, then validates the result. Lowest precedence
first:

| `jaffgen` | Python `Network(...)` |
| --- | --- |
| 1. the network's `jaff.toml` | 1. the network's `jaff.toml` (unless `use_jaff_toml=False`) |
| 2. the `jaffgen.toml` (`--config`, or the one a template bundles) | 2. the `config=` layers, in order |
| 3. command-line options | 3. keyword arguments (`label=`, `strict=`, `radiation_props=`, …) |
| 4. `--set key.path=value` items | |

When two layers set the same table, its keys merge one by one; a list or a
plain value is replaced whole. Tables whose keys depend on their `type` are
replaced whole too: a `shielding` table and the
`[network.thermodynamics.eos]` table, so switching the EOS `type` in one layer
never mixes in keys of another type.

`[jaffgen] network_config` (or `--network-config`) names the `jaff.toml` to use
instead of the one next to the network file.

---

## Global and per-reaction settings

`[network.reactions]` applies to every reaction;
`[network.reactions.reaction."<serialized>"]` applies to one reaction and wins
over the global value of the same key, whichever layer set either:

```toml
[network.reactions]
t_cutoff = "extrapolate"

[network.reactions.photo]
pi_database = "verner"

[network.reactions.reaction."CO._PHOTON__C.O"]
t_cutoff = "clip"

[network.reactions.reaction."CO._PHOTON__C.O".photo]
shielding = { type = "leiden", shielded_by = ["self", "H2"] }
```

<!-- prettier-ignore -->
!!! note "Reaction keys are serialized names"
    The per-reaction key is the reaction's **serialized form**:
    `<reactants>__<products>`, each side a `.`-joined, alphabetically sorted
    list of species names; pseudo-species keep their underscore (`_PHOTON`,
    `_CR`). `CO + PHOTON -> C + O` is `CO._PHOTON__C.O`. The `.` separators
    mean the key **must be quoted**. A key naming a reaction that is not in the
    network is an error.

`photo` settings are only valid on photo reactions.

---

## Paths

A relative path is resolved against the directory of the file it is written
in (`jaff.toml` or `jaffgen.toml`); a path given on the command line or from
Python is resolved against the working directory.

---

## From Python

```python
from jaff import Network

net = Network("GOW")                                    # its jaff.toml, if any
net = Network("GOW", config="tuned.toml")               # a file layer on top
net = Network("GOW", config={"network": {"reactions": {"t_cutoff": "clip"}}})
net = Network("GOW", config=["base.toml", {"network": {"label": "run-2"}}])
net = Network("GOW", strict=True, use_jaff_toml=False)  # keyword arguments win

net.config                                              # the validated config
```

A `config=` layer is a file path, a `dict` shaped like a `jaff.toml`, or a
`JaffTomlConfig`; of a `JaffTomlConfig`, only the keys explicitly set count. A
keyword argument given (not `None`) wins over every layer.
`radiation_props=`, `dust_props=` and `eos_props=` take the matching config
model or a `dict`.

---

## From the command line

`jaffgen` and `jaffx` turn their options into a layer:

| Option | Key |
| --- | --- |
| `--network` | `jaffgen.network_file` |
| `--template` | `jaffgen.template` |
| `--indir` / `--files` | `jaffgen.input_dir` / `jaffgen.input_files` |
| `--outdir` | `jaffgen.output_dir` |
| `--network-config` | `jaffgen.network_config` |
| `--echo` / `--no-echo` | `jaffgen.echo` |
| `--label` | `network.label` |
| `--funcfile` | `network.funcfile` |
| `--expand-nuclei` / `--no-expand-nuclei` | `network.expand_nuclei` |
| `--strict` / `--no-strict` | `network.strict` |
| `--duplicate-policy` | `network.reactions.duplicate_policy` |

`--set key.path=value` (repeatable) sets any key and wins over the other
options. The key path and value are written as in TOML; a value that is not a
TOML value is taken as a string:

```bash
jaffgen --network GOW --template microphysics \
    --set network.reactions.t_cutoff=clip \
    --set 'network.radiation.bands=[6, 13.6, "inf"]' \
    --set 'network.reactions.reaction."CO._PHOTON__C.O".t_cutoff=extrapolate'
```

`jaffx` accepts the `network.*` keys only.

---

## What a run used: `config.toml`

`jaffgen` writes the run's whole config, every key and default after merging
all layers, to `config.toml` in the output directory. It records exactly what
the run used and loads again with `--config`.

---

## Errors

The config is validated after merging. Every problem is reported at once, by
key path:

```text
ParserError: Invalid config:
  network.strict: Input should be a valid boolean
  network.radiation.rsl: unknown key
  network.reactions.t_cutoff: Input should be 'clip' or 'extrapolate'
  network.thermodynamics.eos: unknown EOS type 'polytrope'; choose one of ideal, multi_gamma, fermi_degenerate, relativistic_fermi_degenerate
```

Types, allowed values and unknown keys are also checked in your editor as you
type. Checks across keys (a `profile_index`
list needs one entry per band; a shielding model bound to one reaction) run
when the config loads, and checks against the network (a reaction key that
names no reaction; `photo` settings on a non-photo reaction) when the network
is built.

---

## Migrating from the old keys

Old keys are errors (`unknown key`). They map to:

| Old | New |
| --- | --- |
| `[network] errors` / `--errors` / `Network(errors=)` | `[network] strict` / `--strict` / `Network(strict=)` |
| `[network] config` | `[jaffgen] network_config` |
| `[network] duplicate_policy` | `[network.reactions] duplicate_policy` |
| `[network.rates] T_cutoff` | `[network.reactions] t_cutoff` |
| `[network.reactions."X"] T_cutoff` | `[network.reactions.reaction."X"] t_cutoff` |
| `[network.reactions."X"] pi_database`, `override_custom_rate`, `override_custom_cross_section`, `shielding` | `[network.reactions.reaction."X".photo] …` |
| `[network.radiation] pi_database`, `override_custom_rate`, `override_custom_cross_section`, `use_proxy_photoreaction` | `[network.reactions.photo] …` |
| `[network.radiation] rsl` | `[network.radiation] c` |
| `[network.dust] pe_threshold_low` / `pe_threshold_high` | `[network.dust.pe] threshold_low` / `threshold_high` |
| `[network.eos]` | `[network.thermodynamics.eos]` |
| `[[table]]` | `[[jaffgen.table]]` |
| `RadiationProps` / `DustProps` / `EosProps` | `RadiationConfig` / `DustConfig` / the EOS config models (or a `dict`) |
| a copy of `jaffgen.toml` in the output | `config.toml` |
