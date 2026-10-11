---
tags:
    - Api
---

# Radiation

`jaff.physics.RadiationProps`, `jaff.physics.Radiation`, `jaff.physics.RadiationGroup`

The radiation field is split into contiguous photon-energy **bands**.
`RadiationProps` holds the configuration and checks it, `Network` builds a
`Radiation` from it (`net.radiation`), and `Radiation` creates one
`RadiationGroup` per band. Each group stores its band-averaged quantities and the
rate data for every photo-reaction in that band.

See [Photochemistry](../../../user-guide/designing-networks/photochemistry.md) for the
physics, and [`[network.radiation]`](../../../user-guide/code-generation/jaffgen-toml.md)
for the `jaffgen.toml` equivalent.

## Photon spectrum

Within band $i$ the photon-number spectrum is a power law with that band's own
spectral index $\alpha_i$:

$$
n(E) \propto E^{\alpha_i - 2}, \qquad u(E) = E\,n(E) \propto E^{\alpha_i - 1}
$$

A scalar `profile_index` gives every band the same $\alpha$. A list gives one
$\alpha_i$ per band. Each band is normalised on its own, so the spectrum is a
histogram and is discontinuous at band edges by design. Common values are
$\alpha = 1$ (flat energy spectrum) and $\alpha = 2$ (flat photon spectrum).

---

## RadiationProps

```python
RadiationProps(
    bands=[],
    profile_index=0.0,
    mode="nph",
    c=constants.c.cgs.value,
    background_field="draine",
    pi_database="norad",
)
```

| Argument           | Type                      | Default                 | Description                                                                                                                                                                                                     |
| ------------------ | ------------------------- | ----------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `bands`            | `list[float]` (+ `"inf"`) | `[]`                    | Band edges in eV, ascending, at least two. `"inf"` in the last slot means an open upper band (`sympy.oo`)                                                                                                       |
| `profile_index`    | `float` or `list[float]`  | `0.0`                   | Spectral index $\alpha$. A list must have `len(bands) - 1` entries, one per band                                                                                                                                |
| `mode`             | `str`                     | `"nph"`                 | `"nph"` tracks photon number density (`photden`, cm⁻³); `"u"` tracks energy density (`radeden`, erg cm⁻³)                                                                                                       |
| `c`                | `float` or `str`          | `constants.c.cgs.value` | Speed of light in cm/s. A string (e.g. `"c_hat"`) becomes a symbol, for a reduced speed of light                                                                                                                |
| `background_field` | `str`                     | `"draine"`              | Reference field used to scale `chi_pe`: `bb_4000`, `bb_10000`, `bb_20000`, `draine`, `habing`, `mathis`, `solar`, `tw_hydra`                                                                                    |
| `pi_database`      | `str`                     | `"norad"`               | Photoionization cross-section database: `norad`, `verner` or `leiden` (case-insensitive). Falls back `norad` → `verner` → `leiden` with one summary warning after loading; photodissociation always uses Leiden |

The arguments are validated at construction; the table of `ParserError` /
`RuntimeError` conditions is in
[Valid band / index combinations](../../../user-guide/designing-networks/photochemistry.md#valid-band-index-combinations).
The validated values are stored as attributes of the same names, with `"inf"`
replaced by `sympy.oo` and `mode` / `background_field` / `pi_database` lower-cased.

---

## Radiation

Built by `Network` from `radiation_props`; available as `net.radiation`
(`None` when radiation is off).

### Attributes

| Attribute          | Type                      | Description                                                                              |
| ------------------ | ------------------------- | ---------------------------------------------------------------------------------------- |
| `bands`            | `list`                    | Band edges in eV (from `RadiationProps.bands`)                                           |
| `nbands`           | `int`                     | Number of bands, `len(bands) - 1`                                                        |
| `mode`             | `str`                     | `"nph"` or `"u"`                                                                         |
| `c`                | `float` or `sympy.Symbol` | Speed of light used in `k = c · den · <σ>`                                               |
| `den`              | `sympy.IndexedBase`       | Band density variable, shape `(nbands,)`, named `photden` (`"nph"`) or `radeden` (`"u"`) |
| `groups`           | `list[RadiationGroup]`    | One group per band, in ascending energy order                                            |
| `background_field` | `BackgroundField`         | Reference field built from `RadiationProps.background_field`                             |
| `E_sym`            | `sympy.Symbol`            | Photon-energy symbol `E` (eV) used in the profiles                                       |
| `nph_profile`      | `sympy.Piecewise`         | Photon-number profile over all bands, `E**(α_i - 2)` in band _i_                         |
| `energy_profile`   | `sympy.Piecewise`         | Energy-density profile over all bands, `E**(α_i - 1)` in band _i_                        |
| `photden_tot`      | `float`                   | $\int n(E)\,dE$ over the full range, i.e. the sum of the groups' `photden`               |

In the piecewise profiles a point on an interior edge belongs to the band above
it. Energies below `bands[0]` use the first band's profile and energies above
`bands[-1]` use the last band's.

### Methods

| Method                                    | Returns           | Description                                                                                                        |
| ----------------------------------------- | ----------------- | ------------------------------------------------------------------------------------------------------------------ |
| `get_photden_profile(ph_energy)`          | `numpy.ndarray`   | Evaluates `E**(α_i - 2)` on an energy grid (eV), each point using the index of its band (same edge rules as above) |
| `get_eden_profile(ph_energy)`             | `numpy.ndarray`   | `ph_energy * get_photden_profile(ph_energy)`, i.e. `E**(α_i - 1)`                                                  |
| `ordered_index(idx, order)`               | `tuple[int, int]` | Positions `(density, flux)` of band `idx` in the flat radiation ODE array for layout `order` (0–3)                 |
| `set_reaction_rate_coefficient(reaction)` | `None`            | Band-averages the reaction's tabulated cross section, fills `grp.props[reaction]` and sets `reaction.rate`         |
| `set_custom_rate(reaction)`               | `None`            | Splits a user-supplied `reaction.rate` across bands in proportion to the band integral of `dRad`                   |

`Network` calls the last two methods while building the network; you normally
only read their results. The `order` layouts are:

| `order` | Layout                                     |
| ------- | ------------------------------------------ |
| `0`     | `[den_0, flux_0, den_1, flux_1, ...]`      |
| `1`     | `[flux_0, den_0, flux_1, den_1, ...]`      |
| `2`     | `[den_0, den_1, ..., flux_0, flux_1, ...]` |
| `3`     | `[flux_0, flux_1, ..., den_0, den_1, ...]` |

---

## RadiationGroup

One band, `[lower, upper]` in eV.

| Attribute        | Type                   | Description                                                                              |
| ---------------- | ---------------------- | ---------------------------------------------------------------------------------------- |
| `index`          | `int`                  | Band index in `Radiation.groups`                                                         |
| `sym`            | `sympy.Indexed`        | This band's density variable, `den[index]`                                               |
| `lower`, `upper` | `float` / `sympy.oo`   | Band edges in eV                                                                         |
| `band`           | `tuple`                | `(lower, upper)`                                                                         |
| `dE`             | `float` or `None`      | `upper - lower`; `None` for an open band                                                 |
| `profile_idx`    | `float`                | This band's spectral index $\alpha_i$                                                    |
| `nph_profile`    | `sympy.Expr`           | `E**(profile_idx - 2)`                                                                   |
| `energy_profile` | `sympy.Expr`           | `E**(profile_idx - 1)`                                                                   |
| `photden`        | `float`                | $\int_\text{lower}^\text{upper} n(E)\,dE$, the normalisation for band averages           |
| `eavg`           | `float`                | Band-average photon energy $\langle E\rangle_i = \int E\,n\,dE / \int n\,dE$, in **erg** |
| `props`          | `dict[Reaction, dict]` | Per-reaction band data, keys below                                                       |

`props[reaction]` keys:

| Key         | Description                                                                                                                                                                                                   |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `k`         | Symbolic rate coefficient for this band                                                                                                                                                                       |
| `xsec`      | Photon-number-weighted band-average cross section $\langle\sigma\rangle_i$ (cm²); `None` for custom-rate reactions                                                                                            |
| `xsec_frac` | Tabulated reactions: $\langle\sigma\rangle_i$ divided by the full-spectrum average (`reaction.rad_xsecs`), so values need not sum to 1. Custom-rate reactions: this band's share of the total `dRad` integral |
| `delta_rad` | Band integral of `dRad`, the radiation energy added per reaction event (erg)                                                                                                                                  |

The same values per reaction are available as a table through
[`Reaction.band_xsecs`](../../core/reaction/band_xsecs.md) (there `eavg` is in eV).

---

## Example

```python
import numpy as np
from jaff import Network
from jaff.physics import RadiationProps

net = Network(
    "networks/h_photoionization/h_photo.jet",
    radiation_props=RadiationProps(bands=[13.6, 20.0, 100.0], profile_index=[0, 1]),
)
rad = net.radiation

[grp.profile_idx for grp in rad.groups]   # [0, 1]
[grp.nph_profile for grp in rad.groups]   # [E**(-2), 1/E]
rad.nph_profile                           # Piecewise((E**(-2), E < 20.0), (1/E, True))
float(rad.photden_tot)                    # 1.6329...  (= 0.02353 + ln 5)
rad.get_photden_profile(np.array([15.0, 50.0]))   # [0.00444, 0.02]
rad.ordered_index(1, order=2)             # (1, 3)
```
