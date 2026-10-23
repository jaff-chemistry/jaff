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

The dust module is configured with a `DustProps` model:

```python
from jaff.physics import DustProps

DustProps(
    rv=3.1,
    u_reduction="absorption",
    f_reduction="transport",
    pe_threshold_low=6.0,
    pe_threshold_high=13.6,
)
```

| Parameter           | Type          | Default          | Description                                                                                                                               |
| ------------------- | ------------- | ---------------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| `rv`                | `float`       | `3.1`            | Milky-Way grain model, selected by $R_V$. One of `3.1`, `4.0`, `5.5`; picks the `mw_rv{rv}` table in `dust.hdf5`.                          |
| `u_reduction`       | `str \| None` | `"absorption"`   | Cross-section kind used to attenuate the radiation **density** (u) moment. One of `extinction`, `absorption`, `scattering`, `transport`, or `none`/`None` to disable. |
| `f_reduction`       | `str \| None` | `"transport"`    | Cross-section kind used to attenuate the radiation **flux** (f) moment; same choices as `u_reduction`.                                    |
| `pe_threshold_low`  | `float`       | `6.0`            | Lower edge of the photoelectric band (eV) — the grain work function.                                                                      |
| `pe_threshold_high` | `float`       | `13.6`           | Upper edge of the photoelectric band (eV) — the hydrogen ionisation edge.                                                                 |

The `[pe_threshold_low, pe_threshold_high]` window is the band over which
`chi_pe` is integrated; the defaults span the classic 6–13.6 eV FUV band.

---

## Enabling dust

### From Python

Pass a `DustProps` instance as `dust_props=` to the `Network` constructor,
together with a `radiation_props=` so the radiation field exists:

```python
from jaff import Network
from jaff.physics import RadiationProps, DustProps

net = Network(
    "GOW++",
    radiation_props=RadiationProps(
        bands=[6.0, 11.2, 13.6],   # band edges in eV
        background_field="draine",  # reference field for chi_pe
    ),
    dust_props=DustProps(
        rv=3.1,
        u_reduction="extinction",
        f_reduction="extinction",
    ),
)
```

Omitting `dust_props` (leaving it `None`) leaves the dust module off.

### From `jaffgen.toml`

For code generation, the mere **presence** of a `[network.dust]` table enables
the module. Its keys mirror `DustProps`:

```toml
[network.radiation]
bands            = [6, 11.2, 13.6]
background_field = "draine"

[network.dust]
rv                = 3.1
u_reduction       = "extinction"
f_reduction       = "extinction"
pe_threshold_low  = 6      # eV
pe_threshold_high = 13.6   # eV
```

All keys are optional and fall back to the `DustProps` defaults; an empty
`[network.dust]` table is enough to turn the module on. See the
[`[network.dust]` reference](../configuration/reference.md#networkdust)
in the configuration guide for the full key list, and
[Photochemistry](photochemistry.md) for the radiation setup that dust builds on.

---

## Python API

```python
from jaff import Network
from jaff.physics import RadiationProps, DustProps

net = Network(
    "GOW++",
    radiation_props=RadiationProps(bands=[6.0, 11.2, 13.6]),
    dust_props=DustProps(rv=3.1),
)

# The dust container and its sub-objects
net.dust                      # jaff.physics.dust.Dust (or None if disabled)
net.dust.rv                   # selected R_V grain model
net.dust.pe.chi               # symbolic chi_pe expression
net.dust.tabular.photon_energy  # eV grid backing the cross sections

# Band-averaged dust cross section per H (cm^2/H) over an energy window
net.dust.tabular.avg_cross_section_per_hnuc("absorption", (6.0, 13.6))
```
