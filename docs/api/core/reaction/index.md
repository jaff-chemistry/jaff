---
tags:
    - Api
    - Reaction
---

# Reaction

`jaff.core.reaction.Reaction`

The `Reaction` class represents a single chemical reaction, holding its reactants, products, symbolic rate expression, valid temperature bounds, and associated energy and radiation terms.

## Constructor

`#!python Reaction(net, reactants, products, rate_coefficient, tmin, tmax, original_string, index, dE=None, t_cutoff="clip", type="unknown", rad=None, photo=None)`

**Parameters**

**net** : _Network_
: Network the reaction belongs to. Supplies the `.jfunc` auxiliary functions, the radiation field, the cross-section lookup and the conservation-check mode (`net.spec.errors`).

**reactants** : _list\[Specie\]_
: Reactant species.

**products** : _list\[Specie\]_
: Product species.

**rate_coefficient** : _sympy.Expr_
: Symbolic rate-coefficient expression.

**tmin** : _float or None_
: Minimum valid temperature in Kelvin.

**tmax** : _float or None_
: Maximum valid temperature in Kelvin.

**original_string** : _str_
: Raw reaction string from the network file.

**index** : _int_
: File-side source number of the reaction (its row in the network file); also keys the `.jfunc` functions `chemRate<N>`, `deltaE<N>`.

**dE** : _sympy.Expr or None, optional_
: Internal energy released or absorbed per reaction event (erg). A `.jfunc` `deltaE<N>` function takes precedence; `None` (default) without one means zero.

**t_cutoff** : _str, optional_
: Out-of-range behaviour for the rate segment(s): `"clip"` (default) holds the boundary value, `"extrapolate"` evaluates the expression as-is. See [`jaff.toml` temperature cutoffs](../../../user-guide/configuration/index.md).

**type** : _str, optional_
: Reaction type concluded by the network-format parser, stored verbatim. Generic gas-phase values are `"photo"`, `"cosmic_ray"`, `"3_body"`, `"unknown"`; grain/ice networks add surface-mechanism types (`"freeze"`, `"desorption_thermal"`, `"desorption_cr"`, `"desorption_uvcr"`, `"desorption_h2"`, `"eley_rideal"`, `"eley_rideal_desorption"`, `"langmuir_hinshelwood"`, `"langmuir_hinshelwood_desorption"`, `"h2_formation"`, `"bulk_swap"`, `"surface_swap"`). Default `"unknown"`.

**photo** : _PhotoReactionProps or None, optional_
: Photo-reaction state (see the `photo` attribute). `Network` passes one for every photo reaction. Default `None`.

With a photo reaction, the constructor also looks up its cross sections and, when the network has a radiation field, computes its per-band rate coefficients.

## Attributes

| Attribute             | Type            | Description                                                                                                                                      |
| --------------------- | --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| `reactants`           | `Species`       | Ordered species catalogue of reactant species                                                                                                    |
| `products`            | `Species`       | Ordered species catalogue of product species                                                                                                     |
| `rate_coefficient`    | `sympy.Expr`    | SymPy expression for the rate coefficient (units depend on reaction order; typically cm³ s⁻¹ for two-body reactions). Reassigning it refreshes `rate` |
| `rate`                | `sympy.Expr`    | Reaction rate per unit volume: `rate_coefficient` times the number density of each core reactant (cached)                                         |
| `tmin`                | `float or None` | Minimum gas temperature at which the rate is valid (K). `None` means no lower bound                                                              |
| `tmax`                | `float or None` | Maximum gas temperature at which the rate is valid (K). `None` means no upper bound                                                              |
| `t_cutoff`            | `str`           | Out-of-range behaviour (`"clip"` / `"extrapolate"`) for the rate segment(s)                                                                       |
| `rate_segments`       | `RateSegments`  | Per-temperature-range rate pieces; one segment for a single-range reaction, several when merged across disjoint ranges, collapsed into `rate_coefficient` |
| `dE`                  | `sympy.Basic`   | SymPy expression for the energy released per reaction event (erg)                                                                                |
| `verbatim`            | `str`           | Human-readable `"R1 + R2 -> P1 + P2"` form                                                                                                       |
| `serialized`          | `str`           | Canonical name-level form `"<sorted_reactants>__<sorted_products>"`                                                                              |
| `serialized_exploded` | `str`           | Like `serialized` but built from atom-level serialized forms of each species (isomer-insensitive)                                                |
| `index`               | `int`           | File-side source number (`source_index`): the file-order row number. Stable per parsed row but **gapped** when rows merge (`duplicate_policy` or a temperature-range merge), so it is **not** a valid array subscript. Used for aux lookups (`chemRateN`) and `Reactions.by_source_index()` |
| `catalogue_index`     | `int`           | Dense `0..count-1` position in the parent `Reactions` catalogue, assigned by `Reactions.add()`. Safe as an array subscript in generated code (what the `reaction_idx` template property emits). `-1` until catalogued |
| `type`                | `str`           | Reaction type concluded by the parser (verbatim): gas-phase `"photo"`/`"cosmic_ray"`/`"3_body"`/`"unknown"`, or a grain surface-mechanism type (`"freeze"`, `"desorption_*"`, `"eley_rideal*"`, `"langmuir_hinshelwood*"`, `"h2_formation"`, `"bulk_swap"`, `"surface_swap"`) |
| `photo`               | `PhotoReactionProps or None` | Photo-reaction state; `None` for non-photo reactions. Holds `xsec` (the reaction's [`XsecProps`](../../physics/index.md#xsecprops) cross sections, or `None` when no database has an entry; `xsec.total` is the cross section averaged over the whole spectrum, cm², `None` without radiation), `pi_database` (per-reaction photoionization database override), `override_custom_rate` (per-reaction override of `RadiationProps.override_custom_rate`, or `None`), `override_custom_cross_section` (per-reaction override of `RadiationProps.override_custom_cross_section`, or `None`), `custom_rad_rate` (`True` when the rate came from a `.jfunc` aux function rather than the cross sections), `rgrouped` (one entry per radiation band with that band's `k`, `xsec` and `xsec_frac`, set when a radiation field is configured, `None` for a `chemRate<N>` reaction; see [`band_xsecs`](band_xsecs.md)) and `shielding` (`ShieldingProps` with the shielding model `type` and its `props`, or `None`) |
| `rad`                 | `ReactionRadiationProps` | Radiation the reaction emits: `spectral_yield` (from a `.jfunc` `spectralYield<N>`, erg eV⁻¹ per event as a function of `ph_E`, or `None`) and `group_yields` (its per-group integrals, set when a radiation field is configured; see [Spectral yield](../../../user-guide/designing-networks/auxiliary-functions.md#spectral-yield-spectralyieldn)) |
