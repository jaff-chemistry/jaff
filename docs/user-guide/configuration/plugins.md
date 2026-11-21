---
tags:
    - User-guide
    - Configuration
---

# Plugins

A **plugin** is an opt-in physics model kept in the `plugins/` folder at the root of
the JAFF repository: a shielding factor, a reaction's energy change, a heating or
cooling mechanism, an emission or absorption spectrum. A plugin does nothing until a
`jaff.toml` or `jaffgen.toml` selects it by name, in a table whose `type` is the
plugin's name and whose other keys are its settings:

```toml
[network.thermodynamics]
heating_cooling_rate = [{ type = "my_cooling", species = "CO" }]
```

Plugins are used from a clone of the repository (`plugins/`, `templates/` and
`networks/` are not part of an installed package yet). To write one, see
[Writing plugins](../../development/writing-plugins.md).

## Where each kind is selected

| Kind | Config key | Value | Per | Units |
| --- | --- | --- | --- | --- |
| Shielding | `[network.reactions.photo] shielding` or `[network.reactions.reaction."<s>".photo] shielding` | one table | photo reaction | dimensionless factor |
| Energy change | `[network.reactions.reaction."<s>".thermodynamics] energy_change` | one table | reaction (only per reaction) | erg per reaction event |
| Heating / cooling rate | `[network.thermodynamics] heating_cooling_rate` | list of tables | network | erg cm⁻³ s⁻¹ |
| Emissivity | `[network.radiation] emissivity` | list of tables | network | erg cm⁻³ s⁻¹ eV⁻¹, a spectrum in `ph_E` |
| Absorption coefficient | `[network.radiation] absorption_coefficient` | list of tables | network | cm⁻¹, a spectrum in `ph_E` |

`<s>` is the serialized reaction, e.g. `"H2._PHOTON__H.H"`. Positive energy changes
and rates heat the gas, negative ones cool it.

```toml
[network.reactions.reaction."H.H__H2".thermodynamics]
energy_change = { type = "formation_heating", fraction = 0.33 }

[network.thermodynamics]
heating_cooling_rate = [
    { type = "line_cooling", species = "CO" },
    { type = "line_cooling", species = "H2O" },
]

[network.radiation]
bands = [6.0, 13.6]
emissivity = [{ type = "free_free" }]
absorption_coefficient = [{ type = "dust_extinction", d2g = 0.01 }]
```

(The plugin names above are illustrative: JAFF ships no plugins of these kinds yet.)

## How plugins combine with the `.jfunc` file

An [auxiliary function](../designing-networks/auxiliary-functions.md) can describe the
same physics. The rule depends on whether the quantity belongs to one reaction or to
the whole network:

| Quantity | `.jfunc` function | Rule |
| --- | --- | --- |
| Energy change of reaction *N* | `deltaE<N>` | **The `.jfunc` wins.** Otherwise the plugin, otherwise the value in the network file, otherwise 0. |
| Heating / cooling rate | `heatingCoolingRate` | **Summed**: the `.jfunc` rate plus every listed plugin. |
| Emissivity | `globalEmissivity` | **Summed**, then integrated over each radiation band. |
| Absorption coefficient | `globalAbsorptionCoefficient` | **Summed**, then averaged over each radiation band. |

A per-reaction value is one number, so one source must win: the `.jfunc`, being
written for this network, takes preference. Network-wide mechanisms are independent
physical processes, so they add up: listing a plugin next to a `.jfunc` function
counts both.

## Lists

`heating_cooling_rate`, `emissivity` and `absorption_coefficient` are lists, and every
entry is one term of the sum:

- the same plugin may be listed more than once **with different settings** (e.g. one
  `line_cooling` entry per species): each counts;
- an entry that **repeats another exactly** (same `type`, same settings) is an error,
  as it would count one process twice;
- a list given in a later [config layer](index.md) (e.g. `config=` or a keyword
  argument over the `jaff.toml`) **replaces** the earlier list whole; it is not
  appended to.

## Errors

The config is checked when it is read, with the key path in the message:

- an unknown `type`, or a missing `type`;
- an unknown or ill-typed setting, or a missing required one;
- a reaction-bound plugin selected for another reaction;
- `energy_change` in `[network.reactions]` (it is per reaction only).

## Editor completion

The editor schemas (see [Editor setup](editor-setup.md)) list each plugin's `type` and
settings. They are generated from the plugins present when the schemas were built: after
adding a plugin, regenerate them so that your editor knows it:

```bash
uv run --extra dev python -m jaff.config.json_schema
```

Without a plugin of a kind, the schema accepts any table with a `type`, and JAFF checks
it on load.
