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

---

## Radiation Field Discretisation

The radiation field is divided into contiguous energy bands. You specify the band boundaries in the `[network.radiation]` table of a `jaff.toml` or `jaffgen.toml` (see the [configuration reference](../configuration/reference.md#networkradiation)); the table enables radiation:

```toml
[network.radiation]
bands             = [13.6, "inf"]   # band edges in eV; "inf" for open upper bound
profile_index   = 0               # photon-number spectrum index α (or one per band, e.g. [0, 1])
mode    = "nph"           # photon number density ("nph") or energy density ("u")
c                 = 2.99792458e10   # speed of light (cm/s), or a symbol name such as "c_hat"
```

### Photon spectrum

The photon number density spectrum is assumed to be a power law in photon energy:

$$
n(E) \propto E^{\alpha - 2}
$$

where $\alpha$ is `profile_index`. Setting $\alpha = 0$ gives $n(E) \propto E^{-2}$, i.e. equal energy per logarithmic bin which is the default assumption.

`profile_index` may also be a list with one $\alpha_i$ per band (e.g. `profile_index = [0, 1]` for two bands). Each band $i$ then uses $n(E) \propto E^{\alpha_i - 2}$ and is normalised on its own, so the spectrum is a histogram and is discontinuous at band edges by design.

### Valid band / index combinations

The configuration is checked when it loads, the spectrum when the radiation field is built:

| Condition | Result |
| --- | --- |
| `profile_index` is not an `int`/`float` or a non-empty list of them (`bool` is rejected) | `ParserError` |
| `profile_index` is a list whose length is not `len(bands) - 1` | `ParserError` |
| `bands` has fewer than two edges | `ParserError` |
| `bands` edges are not strictly increasing | `ParserError` |
| `"inf"` appears anywhere but as the last edge of `bands` | `ParserError` |
| `mode = "u"`, `bands[0] = 0` and the **first** band's $\alpha \le 0$ | `RuntimeError` (band-average energy diverges at $E \to 0$) |
| `mode = "u"`, `bands[-1] = "inf"` and the **last** band's $\alpha \ge 0$ | `RuntimeError` (band-average energy diverges at $E \to \infty$) |
| `bands[0] < 1` eV and the first band's $\alpha \le 1$ | warning: the photon-number normalisation $\int E^{\alpha-2}\,dE$ is ill-conditioned near $E \to 0$ |

Only the first band can reach $E = 0$ and only the last can reach $E = \infty$, so the
divergence checks use those bands' indices; interior bands always converge.

### Band-averaged cross section

For band $i$ spanning $[E_\text{lo}, E_\text{hi}]$:

$$
\langle\sigma\rangle_i = \frac{\displaystyle\int_{E_\text{lo}}^{E_\text{hi}} \sigma(E_\nu)\, \dfrac{ \partial n_\gamma}{\partial E_\nu}\, dE_\nu}{\displaystyle\int_{E_\text{lo}}^{E_\text{hi}} \dfrac{ \partial n_\gamma}{\partial E_\nu}\, dE_\nu}
$$

### Rate coefficient

In **photon-density mode** (`mode = "nph"`):

$$
k_i = c \cdot \langle\sigma\rangle_i
$$

In **energy-density mode** (`mode = "u"`):

$$
k_i = \frac{c \cdot \langle\sigma\rangle_i}{\langle E\rangle_i}
$$

where $\langle E\rangle_i$ is the photon-number–weighted average energy in the band.

The rate counts photons, so it always uses the photon-number–weighted
$\langle\sigma\rangle_i$ above: two photons the reaction is equally likely to absorb
contribute equally, whatever their energies.

### Radiation removed from the band

Each reaction absorbs one photon, so in **photon-density mode** the band loses photons
at the reaction rate. In **energy-density mode** the band loses the *energy* of the
absorbed photons, and each takes its own energy $E$, not the band average. The energy
loss therefore weights the cross section by photon energy:

$$
\frac{du_i}{dt} \supset -\,c\, u_i\, \langle\sigma\rangle_i^{(E)} \prod n^{\alpha},
\qquad
\langle\sigma\rangle_i^{(E)} = \frac{\displaystyle\int_{E_\text{lo}}^{E_\text{hi}} \sigma(E)\, E\, \dfrac{\partial n_\gamma}{\partial E}\, dE}{\displaystyle\int_{E_\text{lo}}^{E_\text{hi}} E\, \dfrac{\partial n_\gamma}{\partial E}\, dE}
$$

The two weightings agree for a flat $\sigma$ and differ when $\sigma$ varies across the
band, e.g. at a threshold inside it, where the absorbed photons are more energetic than
the band average. With a separate photo-absorption cross section, the absorption cross
section is the one weighted. Each band's loss coefficient is stored as
`reaction.photo.rgrouped[i].k` (see [`Radiation`](../../reference/physics/radiation/index.md)).

---

## Generated Rate Expressions

After band integration, each photochemical reaction contributes one term per band to the ODE right-hand side:

$$
\frac{dn_i}{dt} \supset -\prod_i {n_i}^{\alpha} \sum_i k_i \cdot \rho_i
$$

The generated code uses the radiation density array as the runtime variable for $\rho_i$: `photden[i]` in photon-density mode (`mode = "nph"`) and `radeden[i]` in energy-density mode (`mode = "u"`). The ODEs for the radiation field itself can be obtained from `#!python Network.sradodes()`; they evolve both the band densities and the band fluxes `rflux[i]`, which lose radiation to the same reactions.

---

## Marking a Reaction as Photochemical

In the `PRIZMO` format, write `PHOTO` (case-insensitive) instead of a rate expression, optionally followed by the band edges in eV (`PHOTO, lo[, hi]` or `photo(lo, hi)`; see [Network Formats](network-formats.md#photoreaction-syntax)):

```text
H -> H+ + E    []    PHOTO, 13.60
```

JAFF will look up the matching cross section in its bundled database and integrate it over the configured radiation bands.

---

## Cross-Section Data

JAFF bundles cross sections from three sources. At network-load time the
serialized reaction key (`Reactant1.Reactant2__Product1.Product2`, the
[serialized form of the reaction](../../reference/core/reaction/index.md); photo
reactions carry the `_PHOTON` agent, e.g. `H._PHOTON__H+.e-`) is looked
up in `jaff.db`, and the cross-section arrays are attached to the reaction's
[`photo.xsec`](../working-with-networks/reactions.md).

| Source                                                                                         | Processes                                           | Notes                        |
| ---------------------------------------------------------------------------------------------- | --------------------------------------------------- | ---------------------------- |
| **Leiden** PDR database ([van Dishoeck et al.](https://home.strw.leidenuniv.nl/~ewine/photo/)) | photoabsorption, photodissociation, photoionization | Tabulated cross sections     |
| **NORAD** / OP ([Nahar, OSU](https://norad.astronomy.osu.edu/))                                | photoionization (ground state)                      | Tabulated; 70 ions of 17 elements (H–Fe) |
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
  SymPy-parseable string (symbol `ph_E`, photon energy in **eV**, σ in cm²).

### What lands on the reaction

Which database supplies a photoionization reaction is governed by
[`pi_database`](../configuration/reference.md#networkreactionsphoto)
(`norad` by default, or `verner` / `leiden`, for every reaction in
`[network.reactions.photo]` or per reaction). The
lookup order is the chosen database, then NORAD, then Verner, then Leiden
(duplicates skipped); any fallback is logged once as a summary after loading.
Photodissociation always uses Leiden.

For tabulated sources (Leiden / NORAD), `reaction.photo.xsec` is an
`XsecProps` object carrying the `photon_energy` grid (eV) plus `photo_absorption`
and the single `photodecay` channel (cm² arrays, or `None`);
`process.decay_type` records whether that channel is ionization or
dissociation. The radiation integrator integrates these arrays numerically over
each band.

Verner fits are kept **symbolic**: the object then carries the SymPy `Piecewise`
in `photo.xsec.photodecay_expr` (symbol `ph_E`, eV, σ in cm²) and the arrays are
`None`; the integrator integrates the expression analytically per band. The
database actually used is recorded in `photo.xsec.database`
(`"norad"`, `"verner"` or `"leiden"`, or `"custom"` for a
[`photoCrossSection<N>`](auxiliary-functions.md#custom-cross-sections-photocrosssectionn)
function).

<!-- prettier-ignore -->
!!! note "Threshold energy"
    The `PHOTO, <lo>[, <hi>]` band edges in the network file are only used
    without radiation transport: the rate is then the placeholder function
    `photorates(n, lo, hi)` (with `n` the reaction's 1-based position among the
    photo reactions), for the
    host code to supply. With radiation enabled the edges are ignored; the
    rate is built from the cross section averaged over every configured band,
    so a band contributes wherever the cross section is non-zero.

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

Shielding is opt-in: a `shielding` table in a reaction's `photo` settings
(`[network.reactions.reaction."<serialized>".photo]`), or in
`[network.reactions.photo]` for every photo reaction (see the
[configuration reference](../configuration/reference.md#networkreactionsphotoshielding)).
The reaction **must be a photo-reaction**; the `type` key selects the shielding
function (default `"leiden"`), and a per-reaction table replaces the global one
whole. The serialized key contains `.` separators and so **must be quoted** in
the TOML header.

```toml
# Leiden tabulated line shielding for CO photodissociation
[network.reactions.reaction."CO._PHOTON__C.O".photo.shielding]
type        = "leiden"
radiation   = "ISRF"           # radiation-field subgroup (default "ISRF")
shielded_by = ["self", "H2"]   # shielding species; "self" = the reactant (CO)

# H2 self-shielding via the Hartwig et al. (2015) fit
[network.reactions.reaction."H2._PHOTON__H.H".photo.shielding]
type      = "hg2015"
min_ncol  = 1.0e-35            # optional floors (see below)
min_vdisp = 1.0e-20
```

The settings are checked against the chosen model when the config loads; a
model bound to one reaction (`db1996`, `hg2015`) is only valid in that
reaction's table. The Leiden model keeps its per-reaction table on the reaction;
`jaffgen` writes it to its output directory as `shielding_<reaction>.hdf5`
(`Network.write_shielding_tables(directory)` from Python).

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
| `min_ncol`  | no       | `1e-50` | Floor on the dimensionless column $x = N_\text{H2}/5\times10^{14}\,\text{cm}^{-2}$ in the fit (dimensionless; not a column density). |
| `min_vdisp` | no       | `1e-50` | Floor on the dimensionless Doppler parameter $b_5$ in the fit, keeping $x/b_5$ finite when `vdisp` is 0 (dimensionless; not a velocity). |

<!-- prettier-ignore -->
!!! note "Shielding is exposed programmatically too"
    The reaction's resolved shielding table is `rxn.config.photo.shielding`,
    with the model in `.type` and the remaining keys in `.params`; it is saved
    in `.jaff` files. `jaff.physics.Photochemistry.shielding(reaction, network)`
    returns the symbolic factor; the radiation integrator evaluates it once per
    reaction and multiplies it into every band. Data tables a model builds for
    the generated code (the Leiden `shielding_<reaction>.hdf5`) are kept in
    `rxn.photo.shielding_tables` and written out by
    `Network.write_shielding_tables(directory)`.

---

## Self-consistent photoelectric field (`chi_pe`)

Grain photoelectric heating is driven by far-UV photons in the **photoelectric
band** — from the grain work function (`6 eV`) up to the hydrogen ionisation
edge (`13.6 eV`). Rate and heating expressions parametrise this by `chi_pe`: the
local radiation field in that band, scaled to a reference (Draine/ISRF) field.

The band edges are the `threshold_low` / `threshold_high` keys of the
[`[network.dust.pe]`](../configuration/reference.md#networkdustpe) table
(defaults `6.0` and `13.6` eV).

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
from jaff import Network
from jaff.physics.dust.config import DustConfig
from jaff.physics.photo_reactions.config import RadiationConfig

net = Network(
    "networks/GOW/GOW.jet",
    radiation=RadiationConfig(
        bands=[6.0, 11.2, 13.6, "inf"],   # band edges in eV
        background_field="draine",         # reference field for the scaling
    ),
    dust=DustConfig(),               # enables the dust module
)
```

In a `jaff.toml` or `jaffgen.toml`, add a
[`[network.dust]`](../configuration/reference.md#networkdust) table next to
`[network.radiation]`; its `background_field` key (the reference field) is optional
and defaults to `"draine"`. If either radiation or dust is missing when a rate
references `chi_pe`, generation aborts with a `ParserError`.

Under a reduced speed of light (`c = "c_hat"` in `[network.radiation]`), `chi_pe`
is scaled by `c_reduced / c` so the photoelectric field stays consistent with
the reduced-speed radiation transport (commit `cddf25d`).

Reference the symbol in any rate or `.jfunc` expression exactly like `chi`; JAFF
substitutes the computed expression when the network is built.

---

## Python API

```python
from jaff import Network
from jaff.physics.photo_reactions.config import RadiationConfig

# Enable photochemistry by declaring radiation bands
net = Network(
    "networks/h_photoionization/h_photo.jet",
    radiation=RadiationConfig(
        bands=[13.6, "inf"],   # band edges in eV
        profile_index=0,
        mode="nph",
    ),
)

# One spectral index per band (here two bands: 13.6–20 eV and 20–100 eV)
net_pb = Network(
    "networks/h_photoionization/h_photo.jet",
    radiation=RadiationConfig(bands=[13.6, 20.0, 100.0], profile_index=[0, 1]),
)
[grp.profile_idx for grp in net_pb.radiation.groups]   # [0, 1]

# Inspect photochemical reactions
photo = net.reactions.photo_reactions()
for rxn in photo:
    print(rxn.verbatim)

# Cross sections are attached at load time
rxn = photo[0]
rxn.photo.xsec.photon_energy       # eV grid
rxn.photo.xsec.photodecay          # cm² array (or None)
rxn.photo.xsec_total               # spectrum-averaged σ (cm²); None without radiation
rxn.photo.band_xsecs               # per-band <σ> and fractions (DataFrame)
rxn.rate_coefficient               # Σ_i k_i, with the radiation density symbols
rxn.photo.plot_xsecs()             # visualise σ(E)

# Symbolic radiation ODEs
print(net.sradodes())
```

Cross-section lookup is also exposed directly via the
`jaff.physics.Photochemistry` class. The reactions catalogue owns one,
`net.reactions.photochemistry`, built on first use; constructing one downloads
the cross-section and line-shielding data files on first use (cached
thereafter), so reuse it:

```python
photo = net.reactions.photochemistry
photo.get_xsec(rxn)         # XsecProps from the tabulated databases
photo.get_verner_xsec(rxn)  # analytic Verner σ(E) (sympy) or None
```
