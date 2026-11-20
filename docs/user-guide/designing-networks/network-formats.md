---
tags:
    - User-guide
    - Network
---

# Network Formats

JAFF can currently parse five file formats that are standard in the astrochemical modelling community. Format detection is automatic — JAFF inspects the file content and picks the correct parser without requiring a format flag.

<!-- prettier-ignore -->
!!! tip "You never declare the format"
    Detection is per line, not per file. JAFF classifies each reaction line on its
    own, so a single file can mix formats freely (see [Using these formats in
    JAFF](#using-these-formats-in-jaff)). The file extension is ignored.

## Supported Formats

| Format      | Origin                                     | Paper                                                               |
| ----------- | ------------------------------------------ | ------------------------------------------------------------------- |
| **KIDA**    | Kinetic Database for Astrochemistry        | [A&A 689, A63 (2024)](https://doi.org/10.1051/0004-6361/202450606)  |
| **UDFA**    | UMIST Database for Astrochemistry (Rate22) | [A&A 682, A109 (2024)](https://doi.org/10.1051/0004-6361/202346908) |
| **PRIZMO**  | Protoplanetary disk photochemistry code    | [MNRAS 494, 4471 (2020)](https://doi.org/10.1093/mnras/staa971)     |
| **KROME**   | Astrophysical chemistry & microphysics lib | [MNRAS 439, 2386 (2014)](https://doi.org/10.1093/mnras/stu114)      |
| **UCLCHEM** | Gas-grain astrochemical Python code        | [AJ 154, 38 (2017)](https://doi.org/10.3847/1538-3881/aa773f)       |

---

## Rate Expression Variables

All rate expressions are parsed as SymPy expressions. The following physical symbols are canonical across all formats (they are the fixed symbols of `NetworkSymbols`):

| Symbol  | Description                                      | Units         |
| ------- | ------------------------------------------------ | ------------- |
| `tgas`  | Gas temperature                                  | K             |
| `tdust` | Dust grain temperature                           | K             |
| `av`    | Visual extinction                                | magnitudes    |
| `crate` | Primary cosmic-ray ionisation rate per H nucleus | s⁻¹           |
| `chi`   | Radiation field strength (Draine 1978 units)     | dimensionless |
| `chi_pe`| Local field scaled to the photoelectric band (6–13.6 eV), from self-consistent radiation transport + dust. See [Photochemistry](photochemistry.md#self-consistent-photoelectric-field-chi_pe) | dimensionless |
| `Zd`    | Grain charge                                     | dimensionless |
| `vdisp` | Velocity dispersion                              | cm s⁻¹        |
| `ph_E`  | Photon energy (variable of cross sections and spectral profiles) | eV |

The following convenience symbols are replaced by density expressions after loading (names are matched case-insensitively):

| Symbol             | Replaced by                                                          |
| ------------------ | -------------------------------------------------------------------- |
| `ntot`             | Total number density (sum of all species densities), cm⁻³            |
| `n_<species>`      | That species' number density, e.g. `n_H2`, `n_e`; charges are written `j` (+) and `k` (−), e.g. `n_Hj` for H⁺ |
| `n_<element>_nuc`  | Nucleus number density of an element, e.g. `n_H_nuc`, cm⁻³          |
| `rc_<N>`           | Rate coefficient of reaction number `<N>` in the network file (see [Auxiliary Functions](auxiliary-functions.md)) |

Any other symbol (for example a dust-to-gas ratio `d2g`) is left as a free variable that you supply at run time.

Format-specific shorthand variables (e.g. `t32`, `te`, `invtgas` from KROME files) are automatically rewritten to the canonical symbols above during parsing.

---

## KIDA Format

**Source:** [kida.astrochem-tools.org/](https://kida.astrochem-tools.org/)

KIDA (Kinetic Database for Astrochemistry) distributes networks as **fixed-width** files: each field occupies a fixed character range, not just a whitespace gap. Each data line encodes one reaction with its Arrhenius parameters, uncertainty estimate, formula index, and temperature range.

<!-- prettier-ignore -->
!!! warning "Columns are fixed-width, not just whitespace-separated"
    The reactant and product blocks are read by character position. Reactants
    occupy the first 34 columns and products the next 57; species names are
    padded with spaces to fit. Keep the column alignment of an existing KIDA
    file intact — adding or dropping spaces shifts the fields and misreads the
    reaction.

### Column layout

```text
R1 R2 R3 | P1 P2 P3 P4 P5 | α  β  γ  F  g  err  itype  Tmin  Tmax  frml  ID  ni  r
```

Reactants occupy columns 1–34 and products columns 35–91; the remaining fields are
whitespace-separated. JAFF reads `α`, `β`, `γ`, `itype`, `Tmin`, `Tmax` and `frml`;
`F`, `g`, `err`, `ID`, `ni` and `r` are skipped.

`frml` selects the rate formula. Only formulas 1–5 produce a rate expression:

| `frml` | Formula                     | Rate                                              |
| ------ | --------------------------- | ------------------------------------------------- |
| 1      | Cosmic-ray ionisation       | `α · crate` (times `2 n_H2/n_H_nuc` when a reactant is `CRP`) |
| 2      | Photodissociation (Draine)  | `α · chi · exp(−γ · av)`                          |
| 3      | Kooij                       | `α (tgas/300)^β exp(−γ/tgas)`                     |
| 4      | ionpol1                     | `α β (0.62 + 0.4767 γ √(300/tgas))`              |
| 5      | ionpol2                     | `α β (1 + 0.0967 γ √(300/tgas) + γ² · 300/(10.526 tgas))` |

Any other `frml` (e.g. `6` Troe fall-off, `7` custom) gives a rate of `0.0`; supply the
rate yourself through a `chemRateCoeff<N>` function in an
[auxiliary function file](auxiliary-functions.md), as `networks/GOW/GOW.jfunc` does.

### Example

A comprehensive example of the KIDA format implementation can be found in `networks/GOW/GOW.jet`

```text
! comment lines start with '!'
! alpha, beta, gamma: Arrhenius coefficients
! itype: reaction type; frml: rate-formula index (see KIDA documentation)
! Tmin, Tmax: valid temperature range

H          CR                     H+         e-                                            0.000e+00  0.000e+00  0.000e+00 2.00e+00 0.00e+00 logn  1  -9999   9999  7     1 1  1
H2         CR                     H2+        e-                                            0.000e+00  0.000e+00  0.000e+00 2.00e+00 0.00e+00 logn  1  -9999   9999  7     2 1  1
He         CR                     He+        e-                                            1.100e+00  0.000e+00  0.000e+00 2.00e+00 0.00e+00 logn  1  -9999   9999  1     3 1  1
```

### Temperature bounds

KIDA writes `-9999` and `9999` for "no limit". JAFF maps them to an open bound
(`None`):

- a `Tmin` of zero or below (e.g. `-9999`) means no lower bound;
- a `Tmax` of exactly `9999` means no upper bound, **unless** another range of
  the same reaction starts above this row's `Tmin`: then `9999` is a real bound
  (e.g. `3000–9999` followed by `10000–20000`);
- every other `Tmax`, including values above `9999` (`15000`, `30000`), is a
  real bound.

---

## UDFA Format

**Source:** [umistdatabase.uk](https://umistdatabase.uk/)

The UMIST Database for Astrochemistry (Rate22) uses a colon-delimited format with a leading integer ID and reaction-type tag. Up to two reactants and four products are supported; unused slots are left empty between consecutive colons.

### Column layout

```text
ID:type:R1:R2:P1:P2:P3:P4:NTR:α:β:γ:Tmin:Tmax:<extra>...[:α:β:γ:Tmin:Tmax:<extra>...]
```

| Field       | Description                                                   |
| ----------- | ------------------------------------------------------------- |
| `ID`        | Integer reaction index                                        |
| `type`      | Reaction class code (e.g. `AD`, `IN`, `RA`, `DR`, `CP`, `CR`, `PH`) |
| `R1–R2`     | Reactants (empty if fewer than 2)                             |
| `P1–P4`     | Products (empty if fewer than 4)                              |
| `NTR`       | Number of temperature ranges (parameter blocks) that follow   |
| `α β γ`     | Arrhenius parameters of one block                             |
| `Tmin Tmax` | Valid temperature range (K) of one block                      |
| `<extra>`   | Per-block metadata (accuracy, references, …); ignored         |

Each of the `NTR` parameter blocks becomes one temperature segment of the reaction, so
a reaction with several ranges is written on a single line. All blocks must have the
same number of fields; double-quoted metadata may contain `:`.

The rate expression depends on `type`: `CP` (direct cosmic-ray ionisation) gives
`α/ζ₀ · crate`, `CR` (cosmic-ray-induced photoreaction) gives
`α γ/(1−ω) (tgas/300)^β / ζ₀ · crate` with `ζ₀ = 1.36×10⁻¹⁷ s⁻¹` and `ω = 0.5`, `PH`
gives `α exp(−γ · av)`, and every other type uses the Kooij form
`α (tgas/300)^β exp(−γ/tgas)`.

A `Tmin` of zero or below means no lower bound, and a `Tmax` of `41000` K or
above (UMIST's upper limit) means no upper bound; both become `None`.

### Example

A comprehensive example of the UDFA format implementation can be found in `networks/rate22_final/rate22_final.rates.jet` (that file currently fails to load as a whole; see [Predefined Networks](predefined-networks.md)). In the lines below the `1` after the products is `NTR`.

```text
1:AD:C-:C:C2:e-:::1:5.00e-10:0.00:0.0:10:41000:L:C:"10.1086/190665":"Prasad and Huntress 1980":
2:AD:C-:CH2:C2H2:e-:::1:5.00e-10:0.00:0.0:10:41000:L:C:"10.1086/190665":"Prasad and Huntress 1980":
3:AD:C-:CH:C2H:e-:::1:5.00e-10:0.00:0.0:10:41000:L:C:"10.1086/190665":"Prasad and Huntress 1980":
```

---

## PRIZMO Format

**Source:** [jaff-chemistry.prizmo](https://github.com/jaff-chemistry/prizmo)

PRIZMO networks are Fortran-flavoured text files with an optional `VARIABLES{}` block at the top. Inside the block, shorthand aliases are defined and then used in the rate expressions on subsequent reaction lines. Comment lines start with `#` or `!`.

### `VARIABLES{}` block

```text
VARIABLES{
    variable_name = expression
    ...
}
```

Variables can reference each other and the standard physical symbols (`tgas`, `av`, etc.). Inside the block, Fortran double-precision literals (`1d0`, `3d2`) and `**` are converted to SymPy automatically; an expression that still fails to parse is skipped with a warning. Assignments outside the block are an error.

### Reaction line syntax

```text
Reactant1 [+ Reactant2 ...] -> Product1 [+ Product2 ...]    [tmin, tmax]    rate_expression
```

The square-bracketed temperature range is mandatory; leave it empty (`[]`) or leave
either bound empty (`[, 300]`) for no limit. A `tmin` of zero or below and a `tmax` of
`1e8` or above also mean no limit. Species are separated by ` + ` (spaces around
the `+`).

<!-- prettier-ignore -->
!!! warning "No Fortran literals in reaction rates"
    The `d` exponent conversion only applies inside `VARIABLES{}`. In a reaction's
    rate expression write `1.2e-10`, not `1.2d-10` — the latter fails to parse. The
    bracketed temperature bounds accept either form.

### Photoreaction syntax

Photo reactions are specified with the `PHOTO` keyword (case-insensitive) as the rate.
It accepts an optional photon-energy band in eV, in any of these forms:

```text
H -> H+ + E    []    PHOTO, 13.60
H -> H+ + E    []    PHOTO, 13.60, 100.0
H -> H+ + E    []    photo(13.60, 100.0)
H -> H+ + E    []    PHOTO
```

The first value is the lower band edge (the photon energy threshold) and the second the
upper edge; a missing edge defaults to `1e99`. A `_PHOTON` pseudo-reactant is added to
the reaction. See [Photochemistry](photochemistry.md) for how the band is used.

### Example

```text
# custom variables can be defined this way
VARIABLES{
    invt  = 1d0 / Tgas
    t32   = Tgas / 3d2
    sqrtt = sqrt(Tgas)
}

H + O -> OH        [10, 1e4]    1.2e-10 * sqrtt
CO -> C + O        []           1.0e-10 * exp(-3.53 * av)
H -> H+ + E        []           PHOTO, 13.60
```

---

## KROME Format

**Source:** [kromepackage.org](https://kromepackage.org/)

KROME files are comma-separated reaction lines whose column semantics are set by an optional `@format:` header. Global Fortran-style variable aliases can be defined with `@var:`. Other KROME declaratives (e.g. `@common:`, `@noTabNext`) are ignored. Comment lines start with `#` or `!`.

### `@format:` header

```text
@format:idx,R,R,R,P,P,P,P,tmin,tmax,rate
```

Tokens: `idx` (index), `R` (reactant), `P` (product), `tmin`/`tmax` (temperature bounds, K), `rate` (reaction rate coefficient). Tokens are case-insensitive. A header applies to every following reaction line until the next header; lines before any header use the default layout shown above (`idx,R,R,R,P,P,P,P,tmin,tmax,rate`). Each line must have exactly as many reactant and product slots as the active layout declares; empty slots are left blank between commas.

Temperature bounds may be written `NONE` (no limit) and may carry Fortran comparison prefixes such as `.LE.` or `>`, which are stripped. Fortran double-precision literals (`1d4`) are accepted in bounds and rates.

### Built-in temperature shortcuts

As in KROME, a rate may use these names without defining them (names are
case-insensitive; `Tgas` is the gas temperature in K):

| Name | Value | Name | Value |
| --- | --- | --- | --- |
| `logT` | `log10(Tgas)` | `sqrTgas` | `sqrt(Tgas)` |
| `lnT` | `ln(Tgas)` | `invsqrT` | `1 / sqrt(Tgas)` |
| `Te` | `Tgas * 8.617343e-5` (eV) | `sqrT32` | `sqrt(T32)` |
| `lnTe` | `ln(Te)` | `invsqrT32` | `1 / sqrt(T32)` |
| `T32`, `t3` | `Tgas / 300` | `Tgas2`, `Tgas3`, `Tgas4` | `Tgas**2`, `Tgas**3`, `Tgas**4` |
| `t4` | `Tgas / 1e4` | `T0` | `288` (K) |
| `invT`, `invTgas` | `1 / Tgas` | `T02`, `T03`, `T04` | `T0**2`, `T0**3`, `T0**4` |
| `invT32` | `1 / T32` | `T0inv` | `1 / T0` |
| `invTe` | `1 / Te` | | |

An `@var` of the same name replaces the built-in value.

### Variable aliases

```text
@var:te      = tgas * 8.617343e-5
@var:invtgas = 1e0 / tgas
@var:t32     = tgas / 3e2
```

<!-- prettier-ignore -->
!!! warning "No inline comments on `@var` lines"
    A trailing comment such as `@var:x = tgas / 3e2   ! note` makes the expression
    unparseable: the variable is skipped with a warning and `x` is left as a free
    symbol in every rate that uses it. Put comments on their own line.

### Example

A comprehensive example can be found at `networks/COthin/react_COthin.jet`

```text
#rates for CO network similar to Glover+2010
@var:Hnuclei = get_Hnuclei(n(:))
@common: user_crate,user_Av,user_Tdust

@var: Te = Tgas*8.617343d-5
@var: invT = 1d0/Tgas

@format:idx,R,R,R,P,P,P,Tmin,Tmax,rate
9,He,H+,,He+,H,,NONE,1d4,1.26d-9*Tgas**(-.75)*exp(-127500d0*invT)
10,He,H+,,He+,H,,1d4,NONE,4d-37*Tgas**4.74
```

The two lines cover adjacent temperature ranges of the same reaction, so JAFF merges
them into one reaction with a piecewise rate.

---

## UCLCHEM Format

**Source:** [uclchem.github.io](https://uclchem.github.io/)

UCLCHEM networks are **comma-separated**, produced by the UCLCHEM Python tool. Every reaction always carries three reactant and four product slots; empty slots are filled with the `NAN` sentinel (this `NAN` is what identifies the format). Gas-phase species have no prefix; ice-surface species are prefixed with `#`, bulk-ice species with `@`.

### Column layout

```text
R1,R2,R3,P1,P2,P3,P4,α,β,γ,Tmin,Tmax,reduced_mass,extrapolate
```

| Field          | Description                                            |
| -------------- | ------------------------------------------------------ |
| `R1–R3`        | Reactants; unused slots are `NAN`                      |
| `P1–P4`        | Products; unused slots are `NAN`                       |
| `α β γ`        | Arrhenius parameters                                   |
| `Tmin Tmax`    | Valid temperature range (K)                            |
| `reduced_mass` | Reduced mass of the reactants                          |
| `extrapolate`  | `True`/`False` — in UCLCHEM, allow rate use outside `[Tmin, Tmax]` |

UCLCHEM switches a rate off only when the gas temperature falls outside
`[Tmin, Tmax]`, and gives the reactions it generates itself (freeze-out,
desorption, ice swaps, …) the placeholder range `0`–`10000` K. A `Tmin` of zero or
below therefore never applies, and JAFF treats it as no lower bound (`None`).

<!-- prettier-ignore -->
!!! warning "`extrapolate` is currently not honoured"
    Due to a known parser issue the flag is never read as `True`: every row keeps
    its own `[Tmin, Tmax]` whatever its `extrapolate` value. The `reduced_mass`
    column is not used either.

### Example

A comprehensive example can be found at `networks/uclchem_small_gas/uclchem_small_gas_network.jet`

```text
!
! Small gas network from UCLCHEM v3.5.1
! Species notation: @ = bulk ice, # = surface ice, no prefix = gas phase
!

H,CRP,NAN,H+,E-,NAN,NAN,5.98e-18,0.0,0.0,10.0,41000.0,0.0,False
C,PHOTON,NAN,C+,E-,NAN,NAN,3.5e-10,0.0,3.76,10.0,41000.0,0.0,False
C,CO+,NAN,CO,C+,NAN,NAN,1.1e-10,0.0,0.0,10.0,41000.0,0.0,False
@C,BULKSWAP,NAN,#C,NAN,NAN,NAN,1.0,0.0,0.0,0.0,10000.0,0.0,False
```

### Mechanism keywords

UCLCHEM encodes the reaction *mechanism* as a keyword in a reactant slot
(`CRP`, `CRPHOT`, `PHOTON`, `FREEZE`, `THERM`, `DESCR`, `DEUVCR`, `DESOH2`,
`LH`, `LHDES`, `ER`, `ERDES`, `H2FORM`, `BULKSWAP`, `SURFSWAP`). JAFF reads this
keyword to (a) set the reaction [`type`](../working-with-networks/reactions.md#reaction-types)
and (b) select the rate-coefficient formula. The radiation and cosmic-ray keywords
`CRP`, `CRPHOT` and `PHOTON` are kept as the pseudo-species `_CRP`, `_CRPHOT` and
`_PHOTON` (e.g. `H,CRP,…` becomes `H + _CRP -> H+ + e-`); the grain-mechanism
keywords (`FREEZE`, `THERM`, `LH`, `BULKSWAP`, …) are dropped from the species
list. Reactions with no keyword are two-body gas-phase reactions using the
Kooij/Arrhenius form.

Not every mechanism has a closed-form, per-line rate. The gas-phase and grain
non-thermal mechanisms (`CRP`, `CRPHOT`, `PHOTON`, `FREEZE`, `DESCR`, `DEUVCR`,
two-body) are emitted as rate expressions; the surface-diffusion and
bulk/surface-swap mechanisms (`THERM`, `DESOH2`, `LH`, `LHDES`, `ER`, `ERDES`,
`H2FORM`, `BULKSWAP`, `SURFSWAP`) need machinery JAFF does not yet expose (dust
temperature, a surface competition formula, or coupling to another reaction's
rate) and are emitted as `0.0`. The full per-mechanism formula table, source
references, and the free symbols each rate expects are documented in the parser's
`README.md` (`src/jaff/core/parsers/network/_formats/uclchem/`).

!!! note "Multiple temperature ranges"
    When the same reaction (same species **and** mechanism) is listed over
    several disjoint `[Tmin, Tmax]` ranges, JAFF merges the rows into one
    reaction with a piecewise rate — see
    [Temperature ranges and piecewise rates](../working-with-networks/reactions.md#temperature-ranges-and-piecewise-rates).
    Overlapping ranges raise a `ValueError`; the same reaction repeated over the
    exact same range is a duplicate, resolved by the
    [duplicate policy](../configuration/reference.md#networkreactions).

---

## Using these formats in JAFF

Since these formats are simple text files, JAFF can read any of these formats irrespective of the file name extension. However, it is recommended to use the `.jet` extension for JAFF networks for the sole purpose of consistency. Network formats in JAFF can also be combined in a single file and JAFF will be able to correctly detect and parse that format.

### How a line's format is recognised

Each line is checked on its own, in this order:

1. **PRIZMO**: the reaction arrow `->` and the `[tmin, tmax]` bracket, before any comma.
2. **UDFA**: the line starts with its index and a colon (`326:CE:...`).
3. **KROME**: any other line with a comma (a `,NAN` field makes it UCLCHEM instead).
4. **UCLCHEM**, then **KIDA** (fixed-width columns).

So a KROME rate may contain `:` (e.g. `get_Hnuclei(n(:))`) and a KROME comment may
contain `->` without the line being taken for UDFA or PRIZMO. A line no format
recognises but that has `->` or `:` is reported as an invalid PRIZMO or UDFA
reaction (e.g. a PRIZMO line missing its `[tmin, tmax]`), not skipped.

### Example

An example of combined reaction formats can be found at `networks/demos/demo1.jet`

```text
VARIABLES{
    invt = 1d0 / Tgas
    t32 = Tgas / 3d2
}

H+ + E -> H                    [0, 300]   3.61e-12*t32**(-0.75)

# KIDA format
N2         CR                     N          N                                             5.000e+00  0.000e+00  0.000e+00 1.25e+00 0.00e+00 logn  1  -9999   9999  1     1 1  1

# UMIST format
326:CE:CO:N2+:N2:CO+:::1:7.40e-11:0.00:0.0:10:41000:M:A:"10.1063/1.438893"::

# KROME format (default)
16,H2,E,,H,H,E,,NONE,NONE,5.6e-11*exp(-102124e0*invT)*Tgas**0.5

# KROME format (custom)
@format:R,R,P,P,rate
N2+,E,N,N,1.7e-7*(Tgas/3e2)**(-0.3)
```
