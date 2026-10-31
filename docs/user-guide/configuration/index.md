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
key is listed in the [reference](reference.md); [editor setup](editor-setup.md)
gives completion, descriptions and error checking while you write the files.

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

Every source of settings is a **layer**; JAFF merges the layers, a later one
winning, then validates the result once. The `jaff.toml` next to a network
holds that network's own defaults and always applies. Lowest precedence first:

| `jaffgen` | Python `Network(...)` |
| --- | --- |
| 1. the network's `jaff.toml` | 1. the network's `jaff.toml` |
| 2. the `jaffgen.toml` (`--config`, or the one a template bundles) | 2. `config=`: one more source (a TOML file, a `dict` of the `[network]` keys, or a `NetworkConfig`) |
| 3. command-line options | 3. keyword arguments: the `[network]` keys (`label=`, `reactions=`, `radiation=`, …) |
| 4. `--set key.path=value` items | |

For `jaffgen`, the `[network]` settings of the `jaffgen.toml` and the command
line (and of a `--network-config` file) are merged over the network's own
`jaff.toml` by the network, which validates them once; `jaffgen` itself checks
the `[jaffgen]` section. Without `--config`, the `jaffgen.toml` a template
bundles is the file layer, merged before anything is checked.

When two layers set the same table, its keys merge one by one; a list or a
plain value is replaced whole. Tables whose keys depend on their `type` are
replaced whole too: a `shielding` table and the
`[network.thermodynamics.eos]` table, so switching the EOS `type` in one layer
never mixes in keys of another type.

`[jaffgen] network_config` (or `--network-config`) names another `jaff.toml`,
merged over the one next to the network file.

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
net = Network("GOW", config="tuned.toml")               # tuned.toml merged over it
net = Network("GOW", config={"reactions": {"t_cutoff": "clip"}})
net = Network("GOW", config=other_net.config, label="run-2")
net = Network("GOW", strict=True)                       # keyword arguments win

net.config                                              # the validated config
```

`config=` is one source merged over the network's `jaff.toml`: a TOML file (its
`[network]` table), a `dict` of the `[network]` table's keys, or a
`NetworkConfig` (only the keys explicitly set on it count). A keyword argument
given (not `None`) is merged over `config=`. The keyword
arguments are the keys of the `[network]` table, named and nested as in the
file: `strict`, `label`, `funcfile`, `expand_nuclei`, and the tables
`reactions`, `radiation`, `dust` and `thermodynamics`, each a `dict` or its
config model (`ReactionsConfig`, `RadiationConfig`, `DustConfig`,
`ThermodynamicsConfig`):

```python
Network(
    "networks/GOW/GOW.jet",
    reactions={"duplicate_policy": "error", "photo": {"use_proxy_photoreaction": True}},
    radiation={"bands": [6.0, 13.6, "inf"]},
    thermodynamics={"eos": {"type": "ideal", "gamma": 1.4}},
)
```

---

## From the command line

`jaffgen` and `jaffx` turn their options into a layer. `jaffx` has only the
options marked in its column:

| Option | Key | `jaffx` |
| --- | --- | --- |
| `--network` | `jaffgen.network_file` | the network to load |
| `--template` | `jaffgen.template` | |
| `--indir` / `--files` | `jaffgen.input_dir` / `jaffgen.input_files` | |
| `--outdir` | `jaffgen.output_dir` | |
| `--network-config` | `jaffgen.network_config` | |
| `--echo` / `--no-echo` | `jaffgen.echo` | |
| `--label` | `network.label` | yes |
| `--funcfile` | `network.funcfile` | yes |
| `--expand-nuclei` / `--no-expand-nuclei` | `network.expand_nuclei` | yes |
| `--strict` / `--no-strict` | `network.strict` | |
| `--duplicate-policy` | `network.reactions.duplicate_policy` | yes |

`--set key.path=value` (repeatable) sets any key and wins over the other
options. The key path and value are written as in TOML; a value that is not a
TOML value is taken as a string:

```bash
jaffgen --network GOW --template microphysics \
    --set network.reactions.t_cutoff=clip \
    --set 'network.radiation.bands=[6, 13.6, "inf"]' \
    --set 'network.reactions.reaction."CO._PHOTON__C.O".t_cutoff=extrapolate'
```

`jaffx` accepts the `network.*` keys only; it has no `--strict`, so use
`--set network.strict=true`.

---

## What a run used: `config.toml`

`jaffgen` writes the run's whole config, every key and default after merging
all layers, to `config.toml` in the run's folder in the output directory. It
records what the run used and loads again with `--config`; its `output_dir` is
the output directory, not the run folder, so a rerun writes a new folder beside
it. As TOML has no null, keys left unset
(`None`) and empty tables are omitted; the file starts with a `#:schema` line
for [editors](editor-setup.md) and puts `[jaffgen]` first.

<!-- prettier-ignore -->
!!! note "Data-table paths in `config.toml`"
    `[[jaffgen.table]]` entries are written as they were given, so a relative
    `source.path` is resolved against the directory of the loaded file: when
    `config.toml` is reloaded from the run folder, that is the run folder
    rather than the original `jaffgen.toml`'s directory. Use absolute source
    paths (or `"default"`) if you plan to rerun from `config.toml`.

---

## Errors

The config is validated after merging. Every problem is reported at once, by
key path. From the command line (`jaffgen`, `jaffx`) it is printed after
`Error:` and the command exits with status 1; from Python, `Network(...)`
raises a `ParserError` with the same message:

```text
Error: Invalid config:
  network.strict: Input should be a valid boolean
  network.reactions.t_cutoff: Input should be 'clip' or 'extrapolate'
  network.radiation.rsl: unknown key
  network.thermodynamics.eos: unknown EOS type 'polytrope'; choose one of ideal, multi_gamma, fermi_degenerate, relativistic_fermi_degenerate
```

Types, allowed values and unknown keys are also checked in your editor as you
type ([editor setup](editor-setup.md)). Checks across keys (band edges must
increase, at least two of them, with `"inf"` only last; a `profile_index` list
needs one entry per band; `[network.dust.pe]` `threshold_low` below
`threshold_high`; a shielding model bound to one reaction) run when the config
loads, and checks against the network (a reaction key that
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
| `Network(radiation_props=…)` / `dust_props=…` | `Network(radiation=…)` / `dust=…` |
| `Network(eos_props=…)` | `Network(thermodynamics={"eos": …})` |
| `Network(duplicate_policy=…)` | `Network(reactions={"duplicate_policy": …})` |
| `Network(use_proxy_photoreaction=…)` | `Network(reactions={"photo": {"use_proxy_photoreaction": …}})` |
| a copy of `jaffgen.toml` in the output | `config.toml` |
| `Network(config="x.toml")` replaced the network's `jaff.toml` | `x.toml` is merged over the `jaff.toml`, which always applies |
| `Network(config={"network": {...}})` | `Network(config={...})`: a `dict` holds the `[network]` keys |
| `Network(config=[a, b])`, `Network(config=JaffTomlConfig(...))` | one source: `Network(config=a, **b)`, or a `NetworkConfig` |
| `Network(use_jaff_toml=False)` | removed: the `jaff.toml` next to a network always applies |
