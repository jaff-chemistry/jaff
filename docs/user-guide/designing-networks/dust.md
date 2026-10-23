---
tags:
    - User-guide
    - Network
---

# Dust

JAFF ships an optional **dust module** that couples interstellar grains to the
network's radiation field. Once enabled it does two things: it drives grain
**photoelectric heating** (via the self-consistent `chi_pe` field) and it
**attenuates the radiation field** as photons are absorbed and scattered by
dust. This page describes the grain model, the tabulated cross sections, and how
to turn the module on.

<!-- prettier-ignore -->
!!! note "Dust needs radiation"
    The dust module is built on top of the radiation subsystem — both
    `chi_pe` and the extinction terms are derived from the network's radiation
    bands and its reference `background_field`. A network with dust but no
    radiation configured aborts with a `ParserError`. See
    [Photochemistry](photochemistry.md) for how the radiation field is set up.

---

## What dust does

The module contributes two distinct pieces of physics.

**Photoelectric heating (`chi_pe`).** Far-UV photons eject electrons from grain
surfaces, heating the gas. The heating rate is parametrised by `chi_pe`: the
local radiation field in the photoelectric band, scaled to a reference
(Draine/ISRF) field. JAFF computes `chi_pe` from the network's own radiation
bands, so photoelectric heating stays consistent with the radiation transport
actually being solved. This is the same `chi_pe` symbol documented in
[Photochemistry → Self-consistent photoelectric field](photochemistry.md#self-consistent-photoelectric-field-chi_pe);
the dust module is what supplies it.

**Radiation extinction.** Dust removes energy from the radiation field. When the
radiation-transport ODEs are generated, each band's density (0th, "u") and flux
(1st, "f") moments lose a dust-attenuation term of the form

$$
\left.\frac{\partial u_i}{\partial t}\right|_\text{dust}
  = -\, Z_d\, c\, u_i\, n_\mathrm{H}\, \langle\sigma\rangle_i^{\,(u)}
$$

and analogously for the flux moment with $\langle\sigma\rangle_i^{\,(f)}$. Here
$n_\mathrm{H}$ is the total hydrogen-nucleus number density, $c$ the (reduced)
speed of light, $Z_d$ a runtime dust-abundance scaling symbol supplied by the
host code, and $\langle\sigma\rangle_i$ the band-averaged dust cross section per
H nucleus. Which cross-section kind attenuates the density and the flux moment
is chosen independently by `u_reduction` and `f_reduction` (set either to
`"none"`/`None` to drop that term).

---

## Grain model and cross sections

Dust optical properties come from a bundled table
(`src/jaff/data/dust/dust.hdf5`). The `rv` parameter selects a Milky-Way grain
model by its total-to-selective extinction ratio $R_V$, reading the matching
`mw_rv{rv}` group (`mw_rv3.1`, `mw_rv4.0`, or `mw_rv5.5`).

The table stores a single **extinction** cross section per H nucleus,
$C_\text{ext}/H$ (cm²/H), along with the scattering albedo $\omega$ and the
scattering asymmetry $\langle\cos\theta\rangle$ over a photon-energy grid (eV).
From these, `jaff.physics.dust.tabular.Tabular` derives the four cross-section
kinds:

| Kind         | Definition                                              | Meaning                                          |
| ------------ | ------------------------------------------------------- | ------------------------------------------------ |
| `extinction` | $C_\text{ext}$                                          | total removal from the beam (absorption + scattering) |
| `absorption` | $(1-\omega)\,C_\text{ext}$                              | energy absorbed by grains                        |
| `scattering` | $\omega\,C_\text{ext}$                                  | energy scattered out of the beam                 |
| `transport`  | $(1-\omega\langle\cos\theta\rangle)\,C_\text{ext}$      | momentum-transfer (transport) cross section      |

For each radiation band, JAFF forms the **band-averaged** cross section per H by
weighting with the radiation photon-density profile over the band edges:

$$
\langle\sigma\rangle_i =
  \frac{\displaystyle\int_{E_\text{lo}}^{E_\text{hi}} \sigma(E)\,\dfrac{\partial n_\gamma}{\partial E}\, dE}
       {\displaystyle\int_{E_\text{lo}}^{E_\text{hi}} \dfrac{\partial n_\gamma}{\partial E}\, dE}
$$

This is the quantity that enters the extinction terms above, so the same band
grid used for photochemistry also controls dust extinction.

---

## Parameters

The dust module is configured by the `[network.dust]` table (in Python, a
`DustConfig` or a `dict` of the same keys); its photoelectric band by the
`[network.dust.pe]` sub-table:

```toml
[network.dust]
rv          = 3.1
u_reduction = "absorption"
f_reduction = "transport"

[network.dust.pe]
threshold_low  = 6.0    # eV
threshold_high = 13.6   # eV
```

| Key                 | Type    | Default        | Description                                                                                                                               |
| ------------------- | ------- | -------------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| `rv`                | `float` | `3.1`          | Milky-Way grain model, selected by $R_V$. One of `3.1`, `4.0`, `5.5`; picks the `mw_rv{rv}` table in `dust.hdf5`.                          |
| `u_reduction`       | `str`   | `"absorption"` | Cross-section kind used to attenuate the radiation **density** (u) moment. One of `extinction`, `absorption`, `scattering`, `transport`, or `none` to disable. |
| `f_reduction`       | `str`   | `"transport"`  | Cross-section kind used to attenuate the radiation **flux** (f) moment; same choices as `u_reduction`.                                    |
| `pe.threshold_low`  | `float` | `6.0`          | Lower edge of the photoelectric band (eV) — the grain work function.                                                                      |
| `pe.threshold_high` | `float` | `13.6`         | Upper edge of the photoelectric band (eV) — the hydrogen ionisation edge.                                                                 |

The `[threshold_low, threshold_high]` window is the band over which
`chi_pe` is integrated; the defaults span the classic 6–13.6 eV FUV band.

---

## Enabling dust

The **presence** of a `[network.dust]` table enables the module; all its keys
are optional, so an empty table is enough. Dust attenuates the radiation field,
so configure `[network.radiation]` too:

```toml
[network.radiation]
bands            = [6, 11.2, 13.6]
background_field = "draine"

[network.dust]
u_reduction = "extinction"
f_reduction = "extinction"
```

This works in a network's `jaff.toml` as well as in a `jaffgen.toml`. See the
[`[network.dust]` reference](../configuration/reference.md#networkdust) for the
full key list, and [Photochemistry](photochemistry.md) for the radiation setup
that dust builds on.

From Python, pass the same settings as `dust_props=` (a `DustConfig` or a
`dict`), together with radiation:

```python
from jaff import Network
from jaff.physics.dust.config import DustConfig
from jaff.physics.photo_reactions.config import RadiationConfig

net = Network(
    "GOW++",
    radiation_props=RadiationConfig(
        bands=[6.0, 11.2, 13.6],   # band edges in eV
        background_field="draine",  # reference field for chi_pe
    ),
    dust_props=DustConfig(u_reduction="extinction", f_reduction="extinction"),
)
```

Without a `[network.dust]` table or `dust_props`, the dust module is off.

---

## Python API

```python
from jaff import Network

net = Network(
    "GOW++",
    radiation_props={"bands": [6.0, 11.2, 13.6]},
    dust_props={"rv": 3.1},
)

# The dust container and its sub-objects
net.dust                      # jaff.physics.dust.Dust (or None if disabled)
net.dust.rv                   # selected R_V grain model
net.dust.pe.chi               # symbolic chi_pe expression
net.dust.tabular.photon_energy  # eV grid backing the cross sections

# Band-averaged dust cross section per H (cm^2/H) over an energy window
net.dust.tabular.avg_cross_section_per_hnuc("absorption", (6.0, 13.6))
```
