---
tags:
    - User-guide
    - Network
---

# Network Configuration (`jaff.toml`)

A `jaff.toml` carries **network-scoped** settings that belong to a network
rather than to a code-generation run — the per-reaction and global
**temperature-cutoff** behaviour and the **duplicate policy**. It lives alongside the network file and is
independent of [`jaffgen.toml`](../code-generation/jaffgen-toml.md), which
configures the `jaffgen` pipeline.

JAFF loads it from, in order:

1. the `config=` argument to [`Network`](network.md) (or `--network-config` on
   [`jaffgen`](../code-generation/jaffgen.md));
2. otherwise a `jaff.toml` auto-detected in the **network file's directory**.

All keys live under the `[network]` table.

---

## Temperature cutoffs

Every reaction carries a temperature validity range `[Tmin, Tmax]`. What happens
when the gas temperature `tgas` falls **outside** that range is controlled by the
reaction's _temperature cutoff_ behaviour:

| Cutoff        | Behaviour outside `[Tmin, Tmax]`                                                    |
| ------------- | ---------------------------------------------------------------------------------- |
| `clip`        | `tgas` is clamped to the nearest bound, so the rate is frozen at its boundary value |
| `extrapolate` | `tgas` is left untouched, so the rate expression is evaluated (extrapolated) as-is  |

`clip` is the default. Under `clip`, `tgas` in each rate expression is replaced by
`max(min(tgas, Tmax), Tmin)` (only the bounds that are defined are applied). Under
`extrapolate`, no clamp is inserted and the raw `tgas` symbol survives into the
generated code.

```toml
# Global default applied to every reaction
[network.rates]
T_cutoff = "clip"          # "clip" (default) or "extrapolate"

# Per-reaction override, keyed by the reaction's serialized form
[network.reactions."CO._PHOTON__C.O"]
T_cutoff = "extrapolate"
```

Values are case-insensitive; anything other than `clip` or `extrapolate` raises a
`ParserError`.

<!-- prettier-ignore -->
!!! note "Reaction keys are serialized names"
    The per-reaction key is the reaction's **serialized form**:
    `<reactants>__<products>`, where each side is a `.`-joined list of species
    names sorted alphabetically. Special pseudo-species keep their underscore
    prefix (e.g. `_PHOTON`, `_CR`). For `CO + PHOTON -> C + O` this is
    `CO._PHOTON__C.O`. The `.` separators mean the key **must be quoted** in the
    TOML table header. The cutoff is baked into the rate expression at load
    time, so a `.jaff` file saved afterwards already carries the resolved
    behaviour and does not need the config on reload.

For a reaction defined over a **single** temperature range, `clip` replaces
`tgas` with `max(min(tgas, Tmax), Tmin)` as above. When a reaction is defined
over **multiple disjoint ranges** (merged into a piecewise rate — see
[Temperature ranges and piecewise rates](reactions.md#temperature-ranges-and-piecewise-rates)),
`T_cutoff` governs only the **outer** edges: below the first range and above the
last, `clip` holds the boundary rate while `extrapolate` lets the end segments
run on. Gaps *between* adjacent ranges are always linearly interpolated,
independent of `T_cutoff`.

---

## Photoionization database

Photoionization reactions may override the global database:

```toml
[network.reactions."C._PHOTON__C+.e-"]
pi_database = "verner"
```

The key is only valid on photoionization reactions (any other reaction raises a
`ParserError`); a value from `jaffgen.toml` wins over `jaff.toml`.

`jaff.toml` has no `[network.radiation]` table. The global choice lives in
[`jaffgen.toml`](../code-generation/jaffgen-toml.md#networkradiation-section)
(`[network.radiation] pi_database`) or, in Python, in
`RadiationProps(bands=[...], pi_database="verner")`. The default is `norad`.

---

## Resolution order

The same settings can also appear in [`jaffgen.toml`](../code-generation/jaffgen-toml.md)
under an identical `[network.rates]` / `[network.reactions."<srxn>"]` schema.
When both files supply a value, they are resolved per reaction as (highest wins):

1. `jaffgen.toml` per-reaction `T_cutoff`
2. `jaff.toml` per-reaction `T_cutoff`
3. `jaffgen.toml` global (`[network.rates]`) `T_cutoff`
4. `jaff.toml` global (`[network.rates]`) `T_cutoff`
5. built-in default `clip`

A **per-reaction** setting always wins over any **global** one; within the same
scope, `jaffgen.toml` overrides `jaff.toml`. In short: `jaffgen.toml` lets a
particular code-generation run override the network's own defaults, without
editing `jaff.toml`.

---

## Duplicate policy

When two rate coefficients share the same reaction, mechanism (`type`), **and**
temperature range `[Tmin, Tmax]`, `duplicate_policy` decides how the clash is
resolved during network construction:

| Policy           | Behaviour                                                     |
| ---------------- | ------------------------------------------------------------ |
| `preserve-first` | Keep the first coefficient seen; drop later duplicates (default) |
| `preserve-last`  | Keep the last coefficient seen; discard earlier ones         |
| `error`          | Raise a `ParserError` — the clash must be fixed in the source network |

```toml
[network]
duplicate_policy = "preserve-last"
```

Resolution order (highest wins):

1. `--duplicate-policy` CLI flag (`jaffgen` / `jaffx`)
2. `jaffgen.toml` `[network].duplicate_policy` (`jaffgen` only)
3. `jaff.toml` `[network].duplicate_policy`
4. built-in default `preserve-first`

Unlike `T_cutoff`, this is a single network-wide key — there is no per-reaction
form.
