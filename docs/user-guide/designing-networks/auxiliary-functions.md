---
tags:
    - User-guide
    - Network
---

# Auxiliary Function Files

An **auxiliary function file** (`.jfunc`) lets you attach custom symbolic expressions to a network without modifying the original network file. You use it to:

- Supply reaction-specific rate coefficients that cannot be expressed with the standard formula.
- Define composite or multi-step reactions whose effective rate coefficient depends on other species densities.
- Provide gas heating and cooling rates that vary with the physical state of the gas.
- Declare global symbolic constants shared across multiple functions.
- Define radiation generation rates per unit photon energy for photo reactions

JAFF automatically looks for a file with the same stem as the network file and a `.jfunc` extension. You can also specify one explicitly:

```python
from jaff import Network

net = Network("networks/GOW/GOW.jet", funcfile="networks/GOW/GOW.jfunc")
```

Although the above examples passes the path to the auxiliary function file explicitly, JAFF automatically detects and parses a `.jfunc` file in the directory which contains the network file. The auxiliary function file **must have the same name as the network file** for this automatic detection to work. So if the network file name is `GOW.jet`, `GOW.jfunc` or `GOW.jet.jfunc` will be detected as an auxiliary functions file.

---

## File Format

A `.jfunc` file contains two kinds of declarations: **global constants** (`@var`) and **named functions** (`@function`).

### Global constants — `@var`

```text
@var name = expression
```

`expression` is a SymPy-compatible expression. It is resolved once and substituted into every function body and variable that references `name`.

```text
@var kB      = 1.380649e-16       # Boltzmann constant  (erg / K)
@var eV      = 1.602177e-12       # 1 eV in erg
@var amu     = 1.66054e-24        # 1 amu in gram
@var TCMB    = 2.73               # CMB temperature at z = 0
@var d2g_solar = 0.01             # Solar dust-to-gas ratio
```

Expressions may reference previously declared `@var` names:

```text hl_lines="3"
@var sigma_SB = 5.670374419e-5
@var c_light  = 2.99892458e10
@var a_R      = 4 * sigma_SB / c_light    # radiation constant
```

### Functions — `@function`

```text
@function name(arg1, arg2, ...)
    # optional per-argument documentation comments
    local_var = expression
    ...
    @return final_expression
```

- Lines inside a function block are **local variable assignments** evaluated in order.
- The block ends with an `@return` statement whose expression is the function's symbolic value.
- Comments (`#`) inside a block may document individual arguments. A comment is attached to an argument's metadata **only when its first word is one of the declared argument names** (`# tgas Gas temperature in K`); any other comment is ignored.
- Long expressions can be wrapped with a trailing backslash `\`.
- `expressions` can contain sympy functions like `Piecewise`

Functions may call any other function defined in the same file; the order of the definitions does not matter.

### Referencing species densities

Functions usually depend on the number densities of species in the network. You
do **not** wire these up by hand — you reference them through reserved symbol
names (typically as the function's arguments), and JAFF substitutes the matching
species density when the network is built. The substitution understands the
following conventions:

| Symbol            | Resolves to                                                            |
| ----------------- | ---------------------------------------------------------------------- |
| `ntot`            | Total number density (sum over **all** species)                        |
| `n_e`             | Density of the electron (`e-`)                                         |
| `n_<element>_nuc` | Element-nucleus sum, e.g. `n_H_nuc`, `n_He_nuc`, `n_C_nuc` (see below) |
| `n_<species>`     | Density of the named species (see suffix rules below)                  |

For the general `n_<species>` form, the part after `n_` is matched against the
network's species by its j/k-normalized identifier:

- `"+"` is written `j` and `"-"` is written `k` — there is no `0` neutral
  marker. So `n_Cj` → `C+`, `n_Hejj` → `He++`, `n_Hk` → `H-`, and plain `n_H` →
  `H` (the single species, not a sum).
- `n_e` → the electron (`e-`).
- any other `n_<species>` resolves to that exact species' density, e.g. `n_CO` → `CO`.

For the element-nucleus form, `n_<element>_nuc` gives the atom-count-weighted
sum of that element's density over every species bearing it (including ions),
e.g. `n_H_nuc`, `n_He_nuc`. A charge suffix on a nucleus sum (`n_Hj_nuc`) is
meaningless and raises `ParserError`, as does an unknown element (`n_Xx_nuc`)
or a species that does not exist in the network (e.g. `n_Ne` in a network
without neon).

Resolution is an exact reverse lookup, not a trailing-character guess, so
multiply-ionized species and element symbols that end in a charge letter (e.g.
Potassium `K`, tin `Sn`) are unambiguous.

The non-density physical symbols are passed through unchanged and supplied by
the solver at runtime. These are the same canonical symbols available in rate
expressions — `tgas`, `av`, `crate`, `chi`, `tdust`, and so on — listed
in [Rate Expression Variables](network-formats.md#rate-expression-variables).
Any other symbol, such as the dust-to-gas ratio `d2g` or `gradv` (the velocity
gradient used by the CO cooling fit) in the GOW example below, is left as a free
variable that the host code supplies at run time.

<!-- prettier-ignore -->
!!! note "`n_<element>_nuc` expansion and `expand_nuclei`"
    By default JAFF expands `n_H_nuc`, `n_He_nuc`, and other `n_<element>_nuc`
    tokens into explicit sums over the element-bearing species. Constructing
    the network with `Network(..., expand_nuclei=False)` instead keeps
    `n_<element>_nuc` as a standalone free symbol (e.g. `nh_nuc`).

### Referencing other rate coefficients — `rc_<N>`

Inside a custom `@function` body you can reference the auto-generated rate
coefficient of another reaction with the reserved symbol `rc_<N>`, where `N` is
the reaction's 0-based **file-order number** — the same number used by
`chemRate<N>` (see [Reaction numbers](#reaction-numbers)). JAFF resolves `rc_<N>`
to that reaction's `rate_coefficient` (looked up with
`network.reactions.by_source_index(N)`, not `network.reactions[N]`), so you can
build one reaction's rate or energy term out of another's without restating its
Arrhenius expression.

```text
@function chemRate7(tgas)
    # Reuse the coefficient of reaction 3 (e.g. a shared temperature fit)
    @return 0.5 * rc_3
```

If no reaction carries number `N`, loading raises a `ParserError`. A malformed
index (`rc_` not followed by an integer, e.g. `rc_x`) only logs an error naming
the offending symbol and the `.jfunc` file; the symbol is left in the expression
as a free variable.

---

## Reserved Function Names

JAFF recognises seven reserved naming conventions:

| Pattern              | Role                                                                                                                                                                                                                                         |
| -------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `chemRate<N>`        | Custom rate coefficient for reaction number _N_ (0-based, file order). Replaces the Arrhenius expression in the network file.                                                                                                                             |
| `deltaE<N>`          | Change in gas thermal energy per occurrence of reaction _N_ (erg). Positive = exothermic.                                                                                                                                                    |
| `photoCrossSection<N>` | Photodecay cross section (cm²) of photo reaction number _N_ as a function of photon energy `ph_E` (eV). Replaces the database cross section (see [Custom cross sections](#custom-cross-sections-photocrosssectionn)). |
| `spectralYield<N>`   | Radiation energy (erg) reaction _N_ adds to the field per event per unit photon energy `ph_E` (eV), integrated over solid angle (see [Spectral yield](#spectral-yield-spectralyieldn)). |
| `globalEmissivity`   | **Volume** emissivity: radiation energy emitted per unit volume, time and photon energy `ph_E` (erg cm⁻³ s⁻¹ eV⁻¹) by processes not tied to a reaction, integrated over solid angle (see [Global emissivity](#global-emissivity-globalemissivity)). |
| `globalAbsorptionCoefficient` | Absorption coefficient of the medium: fraction of the radiation absorbed per unit path length (cm⁻¹) as a function of photon energy `ph_E` (eV) (see [Global absorption coefficient](#global-absorption-coefficient-globalabsorptioncoefficient)). |
| `heatingCoolingRate` | Net gas heating minus cooling rate from all non-chemical radiative processes (erg s⁻¹ cm⁻³).                                                                                                                                                 |
| Any other name       | Helper function — substituted symbolically into the bodies that call it.                                                                                                                                                                     |

<!-- prettier-ignore -->
!!! warning "Custom function naming"
    Once parsed, the functions are stored in lowercase formats to make them case-insensitive. This is to make sure that they are compatible with languages that are case-insensetive and must be kept in mind while naming functions.

    Note that this applies to **function names only**. `@var` constant names are **case-sensitive** — `kB` and `kb` are different symbols — so a constant must be referenced with the exact case it was declared with.

### Reaction numbers

The `<N>` in `chemRate<N>`, `deltaE<N>`, `photoCrossSection<N>`,
`spectralYield<N>` and `rc_<N>` is the reaction's 0-based **file-order number**: its
position among the reaction lines of the network file (`Reaction.index`). It counts
every reaction, not only photo reactions. When several rows merge into one reaction
(multiple temperature ranges, or duplicates), the merged reaction keeps the number of
its first row and the later numbers are skipped, so `<N>` is not the position in
`network.reactions`; use `network.reactions.by_source_index(N)` to look it up.

---

## Complete Example — GOW Network

The GOW (Gong, Ostriker & Wolfire 2017) network ships with a detailed `.jfunc` file. Excerpts below illustrate the full feature set.

### Global constants

```text
@var kB         = 1.380649e-16
@var c_light    = 2.99892458e10
@var sigma_SB   = 5.670374419e-5
@var a_R        = 4 * sigma_SB / c_light
@var eV         = 1.602177e-12
@var d2g_solar  = 0.01
@var TCMB       = 2.73
@var H2_OPR     = 3.0
@var f_oH2      = H2_OPR / (H2_OPR + 1)
@var f_pH2      = 1 - f_oH2
```

### Helper functions

```text
@function kcr_H_fac(n_H_nuc, n_H, n_H2)
    # n_H_nuc  Total H nucleus density (cm^-3)
    # n_H      Neutral H atom density (cm^-3)
    # n_H2     H2 molecule density (cm^-3)
    @return 2.3 * (n_H2/n_H_nuc) + 1.5 * (n_H/n_H_nuc)

@function kgr_gong(tgas, chi, av, n_e, c0, c1, c2, c3, c4, c5, c6)
    # tgas Gas temperature (K)
    # chi  Radiation field (Draine 1978 units)
    # av   V-band extinction (magnitudes)
    # n_e  Free electron density (cm^-3)
    # c0 c1 c2 c3 c4 c5 c6  Fitting coefficients
    psi  = 1.7 * chi * exp(-1.87 * av) * sqrt(tgas) / n_e
    logT = log(tgas)
    k_gr = 1e-14 * c0 / (1 + c1 * psi**c2 * \
        (1 + c3 * tgas**c4 * psi**(-c5 - c6 * logT)))
    @return k_gr
```

### Custom rate coefficients — `chemRate<N>`

Reactions whose rate cannot be written as a simple Arrhenius expression get a `chemRate<N>` function where `N` is the 0-based file-order [reaction number](#reaction-numbers). The arguments must match the symbolic free variables used in the expression.

```text
@function chemRate0(crate, n_H_nuc, n_H, n_H2)
    # Reaction 0: H + CR -> H+ + e-
    # crate  Primary ionisation rate per H nucleon (s^-1)
    @return kcr_H_fac(n_H_nuc, n_H, n_H2) * crate

@function chemRate14(d2g, n_H_nuc, n_H)
    # Reaction 14: H + H -> H2  (grain-assisted)
    # True rate = 3e-17 * n_H * n_H_nuc * (d2g / d2g_solar)
    # JAFF multiplies by the reactant densities automatically,
    # so we derive k by dividing out the extra n_H factor.
    @return 3.0e-17 * (d2g / d2g_solar) * (n_H_nuc / n_H)

@function chemRate15(d2g, tgas, chi, av, n_H_nuc, n_e)
    # Reaction 15: H+ + e- -> H  (grain-assisted)
    c0 = 12.25
    c1 = 8.074e-6
    c2 = 1.378
    c3 = 5.087e2
    c4 = 1.586e-2
    c5 = 0.4723
    c6 = 1.102e-5
    @return kgr_gong(tgas, chi, av, n_e, c0, c1, c2, c3, c4, c5, c6) * \
        (d2g / d2g_solar) * n_H_nuc / n_e
```

With radiation transport enabled, a photo reaction with a `chemRate<N>` keeps that rate but is **not coupled** to the radiation bands (it removes no photons), because a custom rate does not depend on the transported field. JAFF logs a warning naming each such reaction. To couple it, supply a [`photoCrossSection<N>`](#custom-cross-sections-photocrosssectionn) instead, or set [`override_custom_rate`](../configuration/reference.md#networkreactionsphoto) to use the database cross sections without editing the `.jfunc`. A reaction cannot have both a `chemRate<N>` and a `photoCrossSection<N>`.

### Thermal energy change — `deltaE<N>`

Functions labeled by `deltaE<N>` are used for representing the internal energy change per reaction. `N` is the 0-based file-order [reaction number](#reaction-numbers). The expression returned by these functions are multiplied with the reaction rate and added together to obtain the rate of internal energy change.

$$
\frac{dE}{dt} = \sum_{i} R_i \Delta E_i
$$

where $R_i$ represents the rate of the $i^{th}$ reaction and $\Delta E_i$ represents the change in internal energy per reaction

```text
@function deltaE40()
    # Reaction 40: H2 + H -> H + H + H  (collisional dissociation)
    @return -4.48 * eV

@function deltaE13(tgas, n_H_nuc, n_H, n_H2, chi, av)
    # Reaction 13: H2 + photon -> H + H
    # Combined FUV-pumping + kinetic-energy heating (Visser+ 2018)
    k_photo = 5.7e-11 * chi * exp(-4.18 * av) * fshield_H2(ncol_H2)
    f = 1 / (1 + ncrH2(tgas, n_H_nuc, n_H, n_H2, k_photo) / n_H_nuc)
    @return (0.4 + 8 * 2 * f) * eV
```

### Custom cross sections — `photoCrossSection<N>`

A `photoCrossSection<N>` function replaces the database cross section of photo reaction number _N_ (see [Reaction numbers](#reaction-numbers)). Its body is the photodecay cross section in cm² as a function of the photon-energy symbol `ph_E` in eV (a bare `E` is Euler's number to SymPy):

```text
@function photoCrossSection5()
    # Reaction 5: flat 1e-18 cm^2 above 11.26 eV
    @return Piecewise((1e-18, ph_E >= 11.26), (0, True))
```

The reaction stores it as `photo.xsec` with `database = "custom"` and no photo-absorption cross section. With radiation enabled it is averaged over each band and over the whole spectrum exactly like a Verner fit, and the band averages set the photo rate. A reaction cannot have both a `chemRate<N>` and a `photoCrossSection<N>` (`ParserError`), nor a per-reaction `pi_database` while its custom cross section is in use.

Setting [`override_custom_cross_section`](../configuration/reference.md#networkreactionsphoto) `= true` (for every photo reaction in `[network.reactions.photo]`, or for one in its own `photo` table) ignores the function and uses the database cross section instead. Like the global flag, a per-reaction value has an effect only with radiation.

### Spectral yield — `spectralYield<N>`

A `spectralYield<N>` function gives the radiation that reaction _N_ emits: the energy (erg) added to the radiation field per reaction event per unit photon energy, as a function of the photon-energy symbol `ph_E` (eV), so its units are erg eV⁻¹. Supply it already **integrated over solid angle**: the emission is treated as isotropic. Any reaction may have one, photo or not, e.g. radiative recombination:

```text
@function spectralYield7()
    # Reaction 7: H+ + e- -> H, a recombination line near 13.6 eV
    @return Piecewise((2.2e-11 / 0.2, (ph_E >= 13.6) & (ph_E <= 13.8)), (0, True))
```

With radiation transport enabled, JAFF integrates the yield over each radiation group and stores the result on `rxn.rad.group_yields`, in the units of the group density per event:

| radiation `mode`      | Group yield                                  | Units  |
| --------------------- | -------------------------------------------- | ------ |
| `"u"`                 | $\int_\text{group} Y(E)\,dE$                  | erg    |
| `"nph"` (default)     | $\int_\text{group} Y(E)\,/\,E\;dE$ (`E` in erg) | photons |

Each group yield times the reaction rate is added to that group's **density** equation only; isotropic emission carries no net flux, so the flux equations are unchanged. The yield may use other symbols (e.g. `tgas`); such a yield is integrated analytically, so those symbols stay in the generated equations, and it must have a closed-form integral over every group (else `ValueError`). The same holds for a `photoCrossSection<N>`. Energy emitted outside the group range is dropped. Without radiation the yield is stored on `rxn.rad.spectral_yield` but has no effect.

### Heating / cooling function

The heating cooling function is used to add any non-chemical heating and cooling rates to the reaction network. The returned `expression` must be in terms of rate of internal energy change in units of $erg\ s^{-1}\ cm^{-3}$

```text
@function heatingCoolingRate(chi, av, d2g, tgas, n_H_nuc, n_H, n_H2, n_Cj, n_C, n_O, n_CO, n_e, gradv)
    @return heating_grainPE(chi, av, d2g, tgas, n_H_nuc, n_e) \
        - cooling_LyA(n_H, n_e, tgas) \
        - cooling_H2(n_H2, n_H, n_Hj, n_He_nuc, n_e, tgas) \
        - cooling_Cplus(n_Cj, n_H, n_H2, n_e, tgas) \
        - cooling_C0(n_C, n_H, n_H2, n_e, tgas) \
        - cooling_O0(n_O, n_H, n_H2, n_e, tgas) \
        - cooling_CO(n_CO, n_H, n_H2, n_e, tgas, gradv, n_H_nuc) \
        - cooling_dust_coll(chi, av, d2g, tgas, n_H_nuc) \
        - cooling_dust_rec(chi, av, d2g, tgas, n_e)
```

### Global emissivity — `globalEmissivity`

The `globalEmissivity` function is the network-wide counterpart of `spectralYield<N>`: radiation emitted by processes that are not one reaction (e.g. free-free or dust emission). Its body is the **volume emissivity**: the energy emitted per unit volume, per unit time, per unit photon energy, as a function of the photon-energy symbol `ph_E` (eV), in erg cm⁻³ s⁻¹ eV⁻¹, already **integrated over solid angle** (isotropic emission). There is one per network, with no `<N>`:

<!-- prettier-ignore -->
!!! warning "Which emissivity?"
    "Emissivity" means different things in different sources. JAFF expects the
    volume emissivity per unit **photon energy in eV**, integrated over solid angle.
    Convert other conventions first (2.418 × 10¹⁴ Hz eV⁻¹ is $d\nu/dE = 1/h$):

    | Convention                                  | Usual units                       | JAFF `globalEmissivity`                                 |
    | ------------------------------------------- | --------------------------------- | ------------------------------------------------------- |
    | Emission coefficient $j_\nu$, per steradian | erg cm⁻³ s⁻¹ sr⁻¹ Hz⁻¹             | $4\pi\, j_\nu \times 2.418\times10^{14}$               |
    | Angle-integrated $\varepsilon_\nu = 4\pi j_\nu$ | erg cm⁻³ s⁻¹ Hz⁻¹                  | $\varepsilon_\nu \times 2.418\times10^{14}$            |
    | Mass (specific) emissivity                  | erg g⁻¹ s⁻¹ Hz⁻¹                   | $\rho\,\varepsilon_\text{mass} \times 2.418\times10^{14}$ ($\rho$ in g cm⁻³) |
    | Per reaction event                          | erg eV⁻¹                          | not this: use [`spectralYield<N>`](#spectral-yield-spectralyieldn) |
    | Dimensionless (0–1, grey body)               | —                                 | a surface property, not this quantity                   |

    The example below starts from Rybicki & Lightman's angle-integrated free-free
    $\varepsilon_\nu$ (erg cm⁻³ s⁻¹ Hz⁻¹), hence its factor `2.418e14`.

```text
@function globalEmissivity(n_e, n_Hj, tgas)
    # Free-free emission of an ionised hydrogen plasma (Gaunt factor 1)
    @return 6.8e-38 * 2.418e14 * n_e * n_Hj * tgas**-0.5 * exp(-ph_E / (8.617e-5 * tgas))
```

With radiation transport enabled, JAFF integrates it over each radiation group, with the same units rule as the spectral yield (erg cm⁻³ s⁻¹ in `"u"` mode, photons cm⁻³ s⁻¹ in `"nph"` mode), and adds each group's value to that group's **density** equation only. The spectrum is stored on `net.radiation.emissivity`, the group values on `group.emissivity` (with `n_X` symbols resolved to the species densities), and their sum on `net.radiation.total_emissivity`. Other symbols such as `tgas` stay symbolic, which needs a closed-form integral over every group (else `ValueError`). Without radiation it has no effect.

### Global absorption coefficient — `globalAbsorptionCoefficient`

The `globalAbsorptionCoefficient` function gives absorption by the medium itself, beyond the photo reactions (e.g. a grey dust opacity). Its body is the **absorption coefficient** $\alpha$: the fraction of the radiation absorbed per unit path length, in **cm⁻¹**, as a function of the photon-energy symbol `ph_E` (eV). Unlike the emissivity, $\alpha$ is a value *at* each photon energy, not a density per eV. There is one per network, with no `<N>`:

```text
@function globalAbsorptionCoefficient(n_H_nuc)
    # Grey dust absorption: 1e-21 cm^2 per H nucleus at every photon energy
    @return 1e-21 * n_H_nuc
```

<!-- prettier-ignore -->
!!! note "Converting to an absorption coefficient"
    | Given                                  | Units      | $\alpha$                 |
    | -------------------------------------- | ---------- | ------------------------ |
    | Cross section $\sigma$ per absorber    | cm²        | $n\,\sigma$ ($n$ in cm⁻³) |
    | Mass opacity $\kappa$                  | cm² g⁻¹    | $\rho\,\kappa$ ($\rho$ in g cm⁻³) |
    | Optical depth per path length          | cm⁻¹       | $d\tau/ds$               |

With radiation transport enabled, JAFF averages $\alpha$ over each radiation group, weighted by the group's assumed spectrum, and removes the absorbed radiation from that group's density and flux equations:

$$
\frac{d\,den_i}{dt} \supset -\,c\,\langle\alpha\rangle_i\, den_i,
\qquad
\frac{d\,F_i}{dt} \supset -\,c\,\langle\alpha\rangle_i\, F_i,
\qquad
\langle\alpha\rangle_i = \frac{\int_\text{group} \alpha(E)\, w(E)\, dE}{\int_\text{group} w(E)\, dE}
$$

with $w(E) = E\,n(E)$ in `"u"` mode (the energy density loses each photon's energy) and $w(E) = n(E)$ in `"nph"` mode (the photon density loses photons). The averages are stored on `group.absorption_coefficient` (with `n_X` symbols resolved to the species densities) and the spectrum on `net.radiation.absorption_coefficient`. Absorption removes radiation from every direction in proportion to the intensity, so it also attenuates the flux $F_i$ (`rflux[i]`), with the same $\langle\alpha\rangle_i$. Set `absorption_attenuates_flux = false` (in [`[network.radiation]`](../configuration/reference.md#networkradiation), or `radiation={..., "absorption_attenuates_flux": False}` from Python) to apply it to the density equations only. Other symbols such as `tgas` stay symbolic, which needs a closed-form integral over every group (else `ValueError`). Without radiation it has no effect.

<!-- prettier-ignore -->
!!! warning "Dust module"
    With the [dust module](../configuration/reference.md#networkdust) enabled, dust absorption is already applied to the radiation field. Do not put a dust opacity in `globalAbsorptionCoefficient` as well, or dust absorbs twice.

---

## Radiation Moment Equations

When a network is loaded with radiation transport enabled (a [`[network.radiation]`](../configuration/reference.md#networkradiation) table in its config, or `radiation=` passed to `Network`), JAFF builds a set of **radiation moment equations** alongside the chemical ODEs. Photo-reactions remove radiation from each band they absorb in.

### Radiation density variable — `radeden` / `photden`

Each radiation band carries one density unknown in the generated equations. Its symbolic name depends on how the field is tracked, set by the radiation `mode` (`"nph"` or `"u"`):

| radiation `mode`      | Density variable | Quantity tracked         | Units    |
| --------------------- | ---------------- | ------------------------ | -------- |
| `"u"`                 | `radeden[i]`     | Radiation energy density | erg cm⁻³ |
| `"nph"` (default)     | `photden[i]`     | Photon number density    | cm⁻³     |

Here `i` indexes the radiation band. The photo-reaction rate coefficient for band $i$ is $k_i = c \cdot den_i \cdot {\langle \sigma \rangle}_i$, where ${\langle \sigma \rangle}_i$ is the band-averaged photoionisation cross section and `c` is the speed of light.

Reactions with a [`spectralYield<N>`](#spectral-yield-spectralyieldn) and the network [`globalEmissivity`](#global-emissivity-globalemissivity) add their emission to the density equations, and the network [`globalAbsorptionCoefficient`](#global-absorption-coefficient-globalabsorptioncoefficient) removes the absorbed radiation from them. The radiation moment equations themselves are available through [`Network.sradodes()`](../../reference/core/network/sradodes.md).

---

## Continuation Lines

Break long expressions across multiple lines with a trailing backslash:

```text
@function myRate(tgas, ne)
    term1 = 3.4e-8 * exp(-7.62 / tgas)
    term2 = 6.97e-9 * exp(-1.38 / tgas) + \
            1.31e-7 * exp(-26.6 / tgas) + \
            1.51e-4 * exp(-8110 / tgas)
    @return (term1 + tgas**-1.5 * term2) / ne
```

---

## Supported SymPy Functions

Inside `@function` bodies and `@var` expressions you can use any SymPy built-in:

`exp`, `log`, `sqrt`, `sin`, `cos`, `tan`, `Piecewise`, `Min`, `Max`, `Abs`, `Array`, and all standard python arithmetic operators (`+`, `-`, `*`, `/`, `**`).
