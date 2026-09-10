---
tags:
    - User-guide
    - Network
---

# Photochemistry

JAFF provides first-class support for photodissociation and photoionisation reactions. This page explains the physical model, how JAFF reads cross-section data, and how rate coefficients are derived from a discretised radiation field.

---

## Overview

A photochemical reaction has the form

$$
\text{Species} + h\nu \longrightarrow \text{Products}
$$

The reaction rate per unit volume is

$$
\mathcal{R} = \prod_i {n_i}^{\alpha} \int_0^\infty \sigma(E_\nu)\, c\, \dfrac{ \partial n_\gamma}{\partial E_\nu}\, dE_{\nu}
$$

where $\sigma(E_\nu)$ is the reaction cross section ($\text{cm}^2$), $c$ is the speed of light ($\text{cm s}^{-1}$), and $n_\gamma$ is the photon number density ($\text{cm}^{-3}$). The subscript $\gamma$ signifies the local radiation field and a subscript of $\nu$ signifies the frequency of the photon. JAFF uses high-resolution cross-section tables from the [Leiden Observatory Database](https://home.strw.leidenuniv.nl/~ewine/photo/index.html) and the [NORAD Database](https://norad.astronomy.osu.edu/) for photo-dissociaton and photo-ionization reactions respectively to store and calculate band averaged cross-sections.

When the radiation energy density per band is calculated, JAFF replaces the continuous integral with a sum over discrete **radiation frequency bands**:

$$
\mathcal{R} = \prod_i {n_i}^{\alpha} \sum_i k_i \cdot \rho_i
$$

where $\rho_i$ is the radiation energy or photon density in band $i$ and $k_i$ is a pre-computed rate coefficient for that band.

<!-- prettier-ignore -->
!!! tip "Radiation source terms"
    If a reaction adds photons to the local radiation field, it can be specified using a custom `deltaRad<N>` function. Details are mentioned in [Auxiliary functions file](auxiliary-functions.md)

---

## Radiation Field Discretisation

The radiation field is divided into contiguous energy bands. You specify the band boundaries in `jaffgen.toml` (see [Configuration File](../code-generation/jaffgen-toml.md)):

```toml
[network.radiation]
bands             = [13.6, "inf"]   # band edges in eV; "inf" for open upper bound
power_law_index   = 0               # photon-number spectrum index α
energy_density    = false           # use photon density (false) or energy density (true)
rsl               = 2.99792458e10   # speed of light (cm/s)
```

### Photon spectrum

The photon number density spectrum is assumed to be a power law in photon energy:

$$
n(E) \propto E^{\alpha - 2}
$$

where $\alpha$ is `power_law_index`. Setting $\alpha = 0$ gives $n(E) \propto E^{-2}$, i.e. equal energy per logarithmic bin which is the default assumption.

### Band-averaged cross section

For band $i$ spanning $[E_\text{lo}, E_\text{hi}]$:

$$
\langle\sigma\rangle_i = \frac{\displaystyle\int_{E_\text{lo}}^{E_\text{hi}} \sigma(E_\nu)\, \dfrac{ \partial n_\gamma}{\partial E_\nu}\, dE_\nu}{\displaystyle\int_{E_\text{lo}}^{E_\text{hi}} \dfrac{ \partial n_\gamma}{\partial E_\nu}\, dE_\nu}
$$

### Rate coefficient

In **photon-density mode** (`energy_density = false`):

$$
k_i = c \cdot \langle\sigma\rangle_i
$$

In **energy-density mode** (`energy_density = true`):

$$
k_i = \frac{c \cdot \langle\sigma\rangle_i}{\langle E\rangle_i}
$$

where $\langle E\rangle_i$ is the photon-number–weighted average energy in the band.

---

## Generated Rate Expressions

After band integration, each photochemical reaction contributes one term per band to the ODE right-hand side:

$$
\frac{dn_i}{dt} \supset -\prod_i {n_i}^{\alpha} \sum_i k_i \cdot \rho_i
$$

The generated code uses the radiation density array `den[i]` as the runtime variable for $\rho_i$. The ODE for the radiation field itself can be obtained from `#!python Network.sradodes()`.

---

## Marking a Reaction as Photochemical

In the `PRIZMO` format, append `PHOTO` followed by the energy threshold (eV) instead of a rate expression:

```text
H -> H+ + E    []    PHOTO, 13.60
```

JAFF will look up the matching cross section in its bundled database and integrate it over the configured radiation bands.

---

## Cross-Section Data

JAFF bundles cross sections from three sources. At network-load time the
serialized reaction key (`Reactant1.Reactant2__Product1.Product2`, the
[serialized form of the reaction](../../api/core/reaction/index.md); photo
reactions carry the `_PHOTON` agent, e.g. `H._PHOTON__H+.e-`) is looked
up in `jaff.db`, and the cross-section arrays are attached to the reaction's
[`xsecs_dict`](../working-with-networks/reactions.md).

| Source                                                                                         | Processes                                           | Notes                        |
| ---------------------------------------------------------------------------------------------- | --------------------------------------------------- | ---------------------------- |
| **Leiden** PDR database ([van Dishoeck et al.](https://home.strw.leidenuniv.nl/~ewine/photo/)) | photoabsorption, photodissociation, photoionization | Tabulated cross sections     |
| **NORAD** / OP ([Nahar, OSU](https://norad.astronomy.osu.edu/))                                | photoionization (ground state)                      | Tabulated, per ion Z = 1..26 |
| **Verner et al. 1996** ([ADS](https://ui.adsabs.harvard.edu/abs/1996ApJ...465..487V/abstract)) | photoionization                                     | Analytic fits σ(E)           |

### Storage layout

Tabulated cross sections are stored as two collapsed HDF5 files, one group per
serialized reaction, all datasets co-sorted by ascending photon energy. These
files are not bundled in the package: on first network load JAFF downloads them
(via `pooch`) from a remote mirror into the local `data/xsecs` directory, then
reuses the cached copies on subsequent runs:

```text
src/jaff/data/xsecs/leiden.hdf5      # one group per reaction (absorption + decay)
src/jaff/data/xsecs/norad.hdf5       # one group per reaction (ionization decay)
src/jaff/data/xsecs/verner_1996.csv  # analytic-fit parameters
```

`photon_energy` is stored in **eV** and every cross-section dataset
(`photoabsorption` / `photodecay`) in **cm²**, where `photodecay` is the
reaction's single ionization-or-dissociation channel.

These assets feed two SQLite tables in `jaff.db`, which is what JAFF actually
queries at runtime (see [Codebase Structure](../../development/codebase-structure.md)):

- `photo_reaction_cross_sections` — one row per reaction, with a
  `photo_absorption` flag, a `decay_type` (`"ionization"` / `"dissociation"`)
  and a `file.hdf5::<group>` pointer into the Leiden or NORAD HDF5 file.
- `verner_cross_sections` — the Verner analytic σ(E) expression as a
  SymPy-parseable string (symbol `E`, photon energy in erg, σ in cm²).

### What lands on the reaction

For tabulated sources, `reaction.xsecs_dict` is an `XsecsProps` dict carrying
the `photon_energy` grid (eV) plus `photo_absorption` and the single
`photodecay` channel (cm² arrays, or `None`); `_equations["decay_type"]`
records whether that channel is ionization or dissociation. The radiation
integrator reads these arrays directly and integrates them numerically over
each band; for photoionization it falls back to the Verner analytic fit when no
tabulated entry exists.

<!-- prettier-ignore -->
!!! note "Threshold energy"
    The `PHOTO, <eV>` threshold in the network file selects which band edges
    apply to a reaction. Reactions whose threshold lies above the upper edge of
    a band are not assigned a rate coefficient for that band.

---

## Shielding

Photo-rates computed from cross sections assume an _unattenuated_ radiation
field. In a real cloud the photons that drive a reaction are absorbed by
intervening gas — including the dissociating species shielding itself
("self-shielding"). JAFF models this with a dimensionless **shielding factor**
$f_\text{sh} \in [0, 1]$ that multiplies the rate coefficient:

$$
k_i \;\longrightarrow\; k_i \cdot f_\text{sh}
$$

The factor is applied to **every band** of a reaction (it multiplies each
band's $k_i$ before the band sum), so a single $f_\text{sh}$ attenuates the
whole reaction rate. $f_\text{sh}$ is a symbolic expression in runtime
quantities — column densities (`ncol_<species>`) and, for H2, the velocity
dispersion (`vdisp`) — that the host code supplies at integration time.

When a shielding function depends on **several** shielding species, the
per-species factors are **multiplied together**:

$$
f_\text{sh} = \prod_{s} f_s\big(N_s\big)
$$

### Enabling shielding

Shielding is opt-in per reaction, declared in `jaffgen.toml` under
`[network.reactions."<serialized>".shielding]` (see the
[configuration reference](../code-generation/jaffgen-toml.md#networkreactionsserializedshielding-section)).
The reaction **must be a photo-reaction**; the `type` key selects the shielding
function (default `"leiden"`). The serialized key contains `.` separators and so
**must be quoted** in the TOML header.

```toml
# Leiden tabulated line shielding for CO photodissociation
[network.reactions."CO._PHOTON__C.O".shielding]
type        = "leiden"
radiation   = "ISRF"           # radiation-field subgroup (default "ISRF")
shielded_by = ["self", "H2"]   # shielding species; "self" = the reactant (CO)

# H2 self-shielding via the Hartwig et al. (2015) fit
[network.reactions."H2._PHOTON__H.H".shielding]
type      = "hg2015"
min_ncol  = 1.0e-35            # optional floors (see below)
min_vdisp = 1.0e-20
```

### Shielding types

| `type`     | Function                               | Reactions          | Reference                                                                                                                                                                                |
| ---------- | -------------------------------------- | ------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `"leiden"` | Leiden tabulated tables                | any photo-reaction | [Leiden photodissociation database](https://home.strw.leidenuniv.nl/~ewine/photo/); [Heays et al. 2017, A&A 602, A105](https://ui.adsabs.harvard.edu/abs/2017A%26A...602A.105H/abstract) |
| `"db1996"` | H2 self-shielding fit ($\alpha = 2$)   | `H2._PHOTON__H.H`  | [Draine & Bertoldi 1996, ApJ 468, 269](https://ui.adsabs.harvard.edu/abs/1996ApJ...468..269D/abstract) (DOI [10.1086/177689](https://doi.org/10.1086/177689))                            |
| `"hg2015"` | H2 self-shielding fit ($\alpha = 1.1$) | `H2._PHOTON__H.H`  | [Hartwig et al. 2015, MNRAS 452, 1233](https://ui.adsabs.harvard.edu/abs/2015MNRAS.452.1233H/abstract) (DOI [10.1093/mnras/stv1368](https://doi.org/10.1093/mnras/stv1368))              |

#### Leiden tabulated shielding (`type = "leiden"`)

The shielding factor is read from the collapsed Leiden tables
(`data/shielding/leiden.hdf5`, one group per reaction, downloaded on first use
alongside the cross sections). Parameters:

| Key           | Required | Default  | Description                                                                                                     |
| ------------- | -------- | -------- | --------------------------------------------------------------------------------------------------------------- |
| `shielded_by` | **yes**  | —        | List of shielding species. Allowed: `"self"` (the reactant's own column), `"H2"`, `"H"`, `"C"`, `"N2"`, `"CO"`. |
| `radiation`   | no       | `"ISRF"` | Radiation-field subgroup in the table (e.g. `"ISRF"`, `"bb-10000"`, `"Ly-alpha"`).                              |

For each species in `shielded_by`, JAFF emits a per-reaction
`shielding_<reaction>.hdf5` table next to the generated code and one
[interpolation call](../code-generation/table-interpolation.md)
`interp_<index>_shielding_<species>(ncol_<species>)`; the total factor is their
product. `"self"` resolves to the reaction's reactant, so it interpolates over
that species' own column density (e.g. `ncol_CO` for `CO._PHOTON__C.O`).

#### H2 self-shielding (`type = "db1996"` / `"hg2015"`)

Both apply only to the `H2._PHOTON__H.H` (H2 → H + H) reaction and evaluate the
standard three-term analytic fit

$$
f_\text{sh} = \frac{0.965}{(1 + x/b_5)^{\alpha}}
  + \frac{0.035}{(1 + x)^{1/2}}\,
    \exp\!\big(-8.5\times10^{-4}\,(1 + x)^{1/2}\big)
$$

with $x = N_{\text{H2}}/5\times10^{14}\,\text{cm}^{-2}$ and
$b_5 = b/10^5\,\text{cm s}^{-1}$, where the Doppler parameter
$b = \sqrt{2}\,\sigma_v$. The exponent $\alpha$ is the only difference between
the two: $\alpha = 2$ for `db1996`, $\alpha = 1.1$ for `hg2015`. The runtime
inputs are the H2 column density `ncol_H2` and velocity dispersion `vdisp`.

| Key         | Required | Default | Description                                                          |
| ----------- | -------- | ------- | -------------------------------------------------------------------- |
| `min_ncol`  | no       | `1e-50` | Lower floor applied in the fit to avoid a zero denominator (cm⁻²).   |
| `min_vdisp` | no       | `1e-50` | Lower floor applied in the fit to avoid a zero denominator (cm s⁻¹). |

<!-- prettier-ignore -->
!!! note "Shielding is exposed programmatically too"
    `jaff.physics.Photochemistry.shielding(reaction, network)` returns the
    symbolic factor and caches it on
    `reaction._metadata["shielding"]["value"]`; the radiation integrator reuses
    that cached value across bands.

---

## Self-consistent photoelectric field (`chi_pe`)

Grain photoelectric heating is driven by far-UV photons in the **photoelectric
band** — from the grain work function (`6 eV`) up to the hydrogen ionisation
edge (`13.6 eV`). Rate and heating expressions parametrise this by `chi_pe`: the
local radiation field in that band, scaled to a reference (Draine/ISRF) field.

Whereas the plain [`chi`](network-formats.md#rate-expression-variables) symbol is
a runtime input (an assumed external field), `chi_pe` is computed **from the
network's own radiation bands**, so photoelectric heating stays consistent with
the radiation transport actually being solved.

### Definition

`chi_pe` is the ratio of two energy densities in the `[6, 13.6] eV` band:

$$
\chi_\mathrm{pe} = \frac{u_\mathrm{bands}(6\text{–}13.6\,\mathrm{eV})}{u_\mathrm{background}(6\text{–}13.6\,\mathrm{eV})}
$$

- **numerator** — the energy density carried by the network's radiation bands,
  summing over each band the fraction of its energy that overlaps the
  photoelectric band. In photon-density mode each band's `photden` is multiplied
  by the band-average photon energy `eavg` to get an energy density.
- **denominator** — the energy density of the reference `background_field` in the
  same band, obtained from its tabulated photon-flux spectrum.

Both are energy densities in `erg cm⁻³`, so `chi_pe` is dimensionless.

### Enabling it

`chi_pe` needs **radiation transport and the dust module** both on. From Python:

```python
net = Network(
    "networks/GOW/GOW.jet",
    rad_bands=[6.0, 11.2, 13.6, "inf"],  # band edges in eV
    dust=True,                            # enables the dust module
    background_field="draine",            # reference field for the scaling
)
```

From a [`jaffgen.toml`](../code-generation/jaffgen-toml.md#networkdust-section)
add a `[network.dust]` table and a `background_field` key under
`[network.radiation]`. If either radiation or dust is missing when a rate
references `chi_pe`, generation aborts with a `ParserError`.

Reference the symbol in any rate or `.jfunc` expression exactly like `chi`; JAFF
substitutes the computed expression when the network is built.

---

## Python API

```python
from jaff import Network

# Enable photochemistry by declaring radiation bands
net = Network(
    "networks/h_photoionization/h_photo.jet",
    rad_bands=[13.6, float("inf")],   # band edges in eV
    rad_powerlaw_index=0,
    rad_energy_density=False,
)

# Inspect photochemical reactions
photo = net.reactions.photo_reactions()
for rxn in photo:
    print(rxn.verbatim)

# Cross sections are attached at load time
rxn = photo[0]
rxn.xsecs_dict["photon_energy"]     # eV grid
rxn.xsecs_dict["photodecay"]        # cm² array (or None)
rxn.plot_xsecs()                    # visualise σ(E)

# Symbolic radiation ODEs
print(net.sradodes())
```

Cross-section lookup is also exposed directly via the
`jaff.physics.Photochemistry` class. Constructing it downloads the cross-section
and line-shielding data files on first use (cached thereafter), so instantiate
once and reuse:

```python
from jaff.physics import Photochemistry

photo = Photochemistry()
photo.get_xsec(rxn)         # XsecsProps from the tabulated databases
photo.get_verner_xsec(rxn)  # analytic Verner σ(E) (sympy) or None
```
