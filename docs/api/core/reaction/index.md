---
tags:
    - Api
    - Reaction
---

# Reaction

`jaff.core.reaction.Reaction`

The `Reaction` class represents a single chemical reaction, holding its reactants, products, symbolic rate expression, valid temperature bounds, and associated energy and radiation terms.

## Constructor

`#!python Reaction(reactants, products, rate, tmin, tmax, dE, dRad, original_string, index, t_cutoff="clip", type="unknown", errors=False)`

**Parameters**

**reactants** : _list\[Specie\]_
: Reactant species.

**products** : _list\[Specie\]_
: Product species.

**rate** : _sympy.Expr_
: Symbolic rate expression.

**tmin** : _float or None_
: Minimum valid temperature in Kelvin.

**tmax** : _float or None_
: Maximum valid temperature in Kelvin.

**dE** : _sympy.Basic_
: SymPy expression for the internal energy released or absorbed per reaction event (erg).

**dRad** : _sympy.Basic_
: Extra radiation energy density per unit energy (eV) that may be produced/absorbed by the reaction. The energy symbol must be `E` for the quantity to be integrated during cross-section calculation

**original_string** : _str_
: Raw reaction string from the network file.

**index** : _int_
: Position in the network reaction list.

**t_cutoff** : _str, optional_
: Out-of-range behaviour for the rate segment(s): `"clip"` (default) holds the boundary value, `"extrapolate"` evaluates the expression as-is. See [`jaff.toml` temperature cutoffs](../../../user-guide/working-with-networks/jaff-toml.md).

**type** : _str, optional_
: Reaction type concluded by the network-format parser, stored verbatim. Generic gas-phase values are `"photo"`, `"cosmic_ray"`, `"3_body"`, `"unknown"`; grain/ice networks add surface-mechanism types (`"freeze"`, `"desorption_thermal"`, `"desorption_cr"`, `"desorption_uvcr"`, `"desorption_h2"`, `"eley_rideal"`, `"eley_rideal_desorption"`, `"langmuir_hinshelwood"`, `"langmuir_hinshelwood_desorption"`, `"h2_formation"`, `"bulk_swap"`, `"surface_swap"`). Default `"unknown"`.

**errors** : _bool, optional_
: If `True`, terminate the process on mass or charge conservation violations instead of merely logging a warning. Default `False`.

## Attributes

| Attribute             | Type            | Description                                                                                                                                      |
| --------------------- | --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| `reactants`           | `Species`       | Ordered species catalogue of reactant species                                                                                                    |
| `products`            | `Species`       | Ordered species catalogue of product species                                                                                                     |
| `rate`                | `sympy.Expr`    | SymPy expression for the rate coefficient (units depend on reaction order; typically cm³ s⁻¹ for two-body reactions)                             |
| `tmin`                | `float or None` | Minimum gas temperature at which the rate is valid (K). `None` means no lower bound                                                              |
| `tmax`                | `float or None` | Maximum gas temperature at which the rate is valid (K). `None` means no upper bound                                                              |
| `t_cutoff`            | `str`           | Out-of-range behaviour (`"clip"` / `"extrapolate"`) for the rate segment(s)                                                                       |
| `rate_segments`       | `RateSegments`  | Per-temperature-range rate pieces; one segment for a single-range reaction, several when merged across disjoint ranges, collapsed into `rate`     |
| `dE`                  | `sympy.Basic`   | SymPy expression for the energy released per reaction event (erg)                                                                                |
| `dRad`                | `sympy.Basic`   | SymPy expression for the extra photon absorption/emission rate contribution to the radiation moment equations                                    |
| `verbatim`            | `str`           | Human-readable `"R1 + R2 -> P1 + P2"` form                                                                                                       |
| `serialized`          | `str`           | Canonical name-level form `"<sorted_reactants>__<sorted_products>"`                                                                              |
| `serialized_exploded` | `str`           | Like `serialized` but built from atom-level serialized forms of each species (isomer-insensitive)                                                |
| `index`               | `int`           | File-side source number (`source_index`): the file-order row number. Stable per parsed row but **gapped** when rows merge (`duplicate_policy` or a temperature-range merge), so it is **not** a valid array subscript. Used for aux lookups (`chemRateN`) and `Reactions.by_source_index()` |
| `catalogue_index`     | `int`           | Dense `0..count-1` position in the parent `Reactions` catalogue, assigned by `Reactions.add()`. Safe as an array subscript in generated code (what the `reaction_idx` template property emits). `-1` until catalogued |
| `type`                | `str`           | Reaction type concluded by the parser (verbatim): gas-phase `"photo"`/`"cosmic_ray"`/`"3_body"`/`"unknown"`, or a grain surface-mechanism type (`"freeze"`, `"desorption_*"`, `"eley_rideal*"`, `"langmuir_hinshelwood*"`, `"h2_formation"`, `"bulk_swap"`, `"surface_swap"`) |
| `custom_rad_rate`     | `bool`          | `True` when the radiation rate was supplied via a `.jfunc` aux function rather than computed from cross-sections                                 |
| `xsecs_dict`          | `XsecsProps or None` | Photo cross-section data for the reaction's single decay channel. Holds `units`, `_equations` (`pa` photo-absorption flag and `decay_type`, either `"ionization"` or `"dissociation"`), `photon_energy` (eV), the `photo_absorption` / `photodecay` arrays (cm², or `None` where absent; both `None` for Verner), `database` (`"norad"` / `"verner"` / `"leiden"`: the source used) and `photodecay_expr` (Verner SymPy expression in `E`, eV; else `None`). `None` for non-photo reactions |
| `rad_groups`          | `list[RadiationGroup]` | Back-references to the radiation bands this reaction contributes to, populated when a radiation field is configured (empty otherwise). See the [`band_xsecs`](band_xsecs.md) property for the band-averaged cross sections |
