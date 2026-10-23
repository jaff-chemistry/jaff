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
physics, and [`[network.radiation]`](../../../user-guide/configuration/index.md)
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
    override_custom_rate=False,
    override_custom_cross_section=False,
    absorption_attenuates_flux=True,
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
| `override_custom_rate` | `bool`                | `False`                 | When `True`, a photo reaction with a `.jfunc` `chemRate<N>` has its rate recomputed from its cross sections instead. A reaction's own `override_custom_rate` (in `jaffgen.toml` / `jaff.toml`) takes precedence |
| `override_custom_cross_section` | `bool`       | `False`                 | When `True`, a photo reaction with a `.jfunc` `photoCrossSection<N>` uses the database cross section instead. A reaction's own `override_custom_cross_section` takes precedence |
| `absorption_attenuates_flux` | `bool`          | `True`                  | Whether the `globalAbsorptionCoefficient` also attenuates the flux equations ($-c\,\langle\alpha\rangle_i\,rflux_i$); `False` applies it to the density equations only |

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
| `ph_E`             | `sympy.Symbol`            | Photon-energy symbol `ph_E` (eV) used in the profiles                                    |
| `nph_profile`      | `sympy.Piecewise`         | Photon-number profile over all bands, `E**(α_i - 2)` in band _i_                         |
| `energy_profile`   | `sympy.Piecewise`         | Energy-density profile over all bands, `E**(α_i - 1)` in band _i_                        |
| `photden_tot`      | `float`                   | $\int n(E)\,dE$ over the full range, i.e. the sum of the groups' `photden`               |
| `emissivity`       | `sympy.Expr` or `None`    | The network `.jfunc` `globalEmissivity`: the volume emissivity (erg cm⁻³ s⁻¹ eV⁻¹ in `ph_E`, integrated over solid angle), or `None` |
| `total_emissivity` | `sympy.Expr` or `None`    | Property: the sum of the groups' `emissivity`, i.e. its integral over the whole band range |
| `absorption_coefficient` | `sympy.Expr` or `None` | The network `.jfunc` `globalAbsorptionCoefficient` (cm⁻¹ as a function of `ph_E`), or `None` |

In the piecewise profiles a point on an interior edge belongs to the band above
it. Energies below `bands[0]` use the first band's profile and energies above
`bands[-1]` use the last band's.

### Methods

| Method                                    | Returns           | Description                                                                                                        |
| ----------------------------------------- | ----------------- | ------------------------------------------------------------------------------------------------------------------ |
| `get_photden_profile(ph_energy)`          | `numpy.ndarray`   | Evaluates `E**(α_i - 2)` on an energy grid (eV), each point using the index of its band (same edge rules as above) |
| `get_eden_profile(ph_energy)`             | `numpy.ndarray`   | `ph_energy * get_photden_profile(ph_energy)`, i.e. `E**(α_i - 1)`                                                  |
| `ordered_index(idx, order)`               | `tuple[int, int]` | Positions `(density, flux)` of band `idx` in the flat radiation ODE array for layout `order` (0–3)                 |
| `integrate_over_groups(spectrum)`         | `list`            | Integrates an emitted spectrum in `ph_E` (a reaction's `spectralYield<N>`, or the network `globalEmissivity`) over each group, in the units of the group density: `∫ S dE` in `"u"` mode, `∫ S/E dE` (photons) in `"nph"` mode |
| `average_over_groups(spectrum)`           | `list`            | Averages a coefficient defined at each photon energy (the `globalAbsorptionCoefficient`, cm⁻¹) over each group, weighted by $E\,n(E)$ in `\"u\"` mode and $n(E)$ in `\"nph\"` mode |
| `get_group_average_cross_section(xsec, ph_E, group, ph_profile=None, energy_weighted=False)` | `float` | Mean of a tabulated or symbolic σ over one band, weighted by $n(E)$ (normalised by `photden`) or, with `energy_weighted`, by $E\,n(E)$ (normalised by `eden`) |
| `get_total_average_cross_section(xsec, ph_E, ph_profile=None)` | `float` | Photon-number-weighted mean of σ over the whole band range |
| `get_reaction_rate_props(reaction)`        | `tuple`           | Band-averages the reaction's cross section; returns `(k_tot, xsec_tot, grouped)`: the total rate coefficient, the spectrum-averaged cross section (cm²) and one `GroupedPhotoReactionProps` per band |

`Reaction` calls the last two methods while it is built and stores the results
on `reaction.rate_coefficient`, `reaction.photo.xsec.total` and
`reaction.photo.rgrouped`; you normally only read those. The `order` layouts are:

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
| `eden`           | `float`                | $\int_\text{lower}^\text{upper} E\,n(E)\,dE$ (eV), the normalisation for energy-weighted band averages |
| `eavg`           | `float`                | Band-average photon energy $\langle E\rangle_i = \int E\,n\,dE / \int n\,dE$, in **erg** |
| `emissivity`     | `sympy.Expr` or `None` | The network `globalEmissivity` integrated over this band, in the band density's units (erg cm⁻³ s⁻¹ in `"u"` mode, photons cm⁻³ s⁻¹ in `"nph"`); added to its density equation |
| `absorption_coefficient` | `sympy.Expr` or `None` | The network `globalAbsorptionCoefficient` averaged over this band (cm⁻¹), energy-weighted in `\"u\"` mode and photon-number weighted in `\"nph\"`; removed from its density equation as $c\,\langle\alpha\rangle_i\,den_i$ |

Per-reaction band data is not stored on the group: each photo reaction holds
one `GroupedPhotoReactionProps` per band in `reaction.photo.rgrouped`, in the
same order as `Radiation.groups` (`None` for a reaction with a `.jfunc`
`chemRate<N>`, which is not coupled to the bands), with these attributes:

| Attribute   | Description                                                                                                                                                                                                   |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `k`         | Band loss coefficient: times the reactant densities, the radiation the reaction removes from the band per unit time. It is $c \cdot den_i \cdot \langle\sigma\rangle_i$ (with any shielding), with the absorption cross section photon-number weighted in `"nph"` mode and **energy** weighted in `"u"` mode. The reaction rate itself always uses the photon-number weighting |
| `xsec`      | Photon-number-weighted band-average cross section $\langle\sigma\rangle_i$ (cm²); `None` for custom-rate reactions                                                                                            |
| `xsec_frac` | Tabulated reactions: $\langle\sigma\rangle_i$ divided by the full-spectrum average (`reaction.photo.xsec.total`), so values need not sum to 1. |

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
