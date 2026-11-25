---
tags:
    - Development
icon: phosphor/cube
---

# Writing Plugins

A plugin is a Python class in the repo-root `plugins/` folder that registers itself
under a name; a user selects it from a `jaff.toml` / `jaffgen.toml` table whose `type`
is that name (see [Plugins](../user-guide/configuration/plugins.md) for the config
side and how plugins combine with the `.jfunc` file). Every kind follows one pattern:

1. a settings model, a subclass of the kind's `...Params` class (a strict
   [config model](adding-configuration.md));
2. the plugin class, a subclass of the kind's `...Function` base, with a `name`, a
   `config_model()` returning the settings model and one method returning a SymPy
   expression;
3. the `register` decorator of the kind;
4. an import of the class in its folder's `__init__.py`.

## Kinds

| Kind | Import from | Base / settings | Method | Returns | Folder |
| --- | --- | --- | --- | --- | --- |
| Shielding | `jaff.physics.photo_reactions.shielding` | `ShieldingFunction` / `ShieldingParams` | `get_shielding(reaction)` | `ShieldingResult` (dimensionless factor) | `plugins/reactions/<r>/shielding/`, `plugins/physics/radiation/shielding/` |
| Energy change | `jaff.physics.thermodynamics.energy_change` | `EnergyChangeFunction` / `EnergyChangeParams` | `energy_change(reaction)` | erg per reaction event | `plugins/reactions/<r>/thermodynamics/`, `plugins/physics/thermodynamics/` |
| Heating / cooling rate | `jaff.physics.thermodynamics.heating_cooling_rate` | `HeatingCoolingRateFunction` / `HeatingCoolingRateParams` | `heating_cooling_rate(net, params)` | erg cm⁻³ s⁻¹ | `plugins/physics/thermodynamics/` |
| Emissivity | `jaff.physics.photo_reactions.emissivity` | `EmissivityFunction` / `EmissivityParams` | `emissivity(net, params)` | erg cm⁻³ s⁻¹ eV⁻¹ in `ph_E` | `plugins/physics/radiation/emissivity/` |
| Absorption coefficient | `jaff.physics.photo_reactions.absorption_coefficient` | `AbsorptionCoefficientFunction` / `AbsorptionCoefficientParams` | `absorption_coefficient(net, params)` | cm⁻¹ in `ph_E` | `plugins/physics/radiation/absorption_coefficient/` |

`<r>` is the sanitised serialized reaction (`.` → `_`, `+` → `j`, `-` → `k`; e.g.
`H2._PHOTON__H.H` → `H2__PHOTON__H_H`). Shielding has its own page,
[Adding a Custom Shielding Function](adding-shielding-functions.md).

**Per-reaction kinds** (shielding, energy change) set `reaction`: the serialized
reaction the plugin is bound to, or `None` for one any reaction may select. A
reaction's own plugin wins over a general one of the same name. The method receives
the reaction and reads its settings with `self.params(reaction)`.

**Network-wide kinds** (heating / cooling rate, emissivity, absorption coefficient)
belong to no reaction (`register` rejects a `reaction`). One type may be listed several
times with different settings, so the method receives the entry's settings as
`params`.

## How plugins are found

Each kind has a registry, filled on its first lookup: it imports the whole `plugins/`
folder as the package `jaff.plugins` (`jaff.common.plugins.load_plugins`), every folder
at any depth whose name does not start with `_`. Only a folder's `__init__.py` is
imported, so:

- every folder on the way needs an `__init__.py`;
- the `__init__.py` must import the plugin's module (`from .line_cooling import
  LineCooling`), or its `@register` never runs;
- helpers go in a folder starting with `_` (e.g. `_utils/`), which is not imported
  as a plugin; import jaff absolutely (`from jaff... import ...`) and your own helpers
  relatively.

## Expressions

A plugin returns a [SymPy](https://www.sympy.org) expression in the network's symbols
(`net.symbols`, or `reaction.net.symbols`):

| Symbol | Meaning |
| --- | --- |
| `net.symbols.tgas` | Gas temperature (K) |
| `net.symbols.ndens[i]` | Number density of species *i* (cm⁻³); `i = net.species.from_name("CO").index` |
| `net.symbols.ntot`, `net.symbols.n_hnuc` | Total and H-nuclei number densities (cm⁻³) |
| `net.symbols.av`, `chi`, `crate` | Visual extinction, UV field, cosmic-ray ionisation rate |
| `net.symbols.ph_E` | Photon energy (eV): the variable of emissivity and absorption spectra |

Any other free symbol becomes an input of the generated code. Emissivity and
absorption spectra are functions of `ph_E`: JAFF integrates (emissivity) or averages
(absorption) them over each radiation band, so the expression must have a closed-form
integral over the band edges.

## Example: energy change of a reaction

Bound to `H + H -> H2`, in `plugins/reactions/H_H__H2/thermodynamics/`:

```python
from typing import cast

from sympy import Expr, Float

from jaff.physics.thermodynamics.energy_change import (
    EnergyChangeFunction,
    EnergyChangeParams,
    register,
)
from jaff.physics.units import eV_to_erg


class FormationHeatingParams(EnergyChangeParams):
    """Settings of the H2 formation heating.

    Attributes
    ----------
    fraction : float
        Fraction of the 4.48 eV binding energy that heats the gas.
    """

    fraction: float = 1.0 / 3.0


@register
class FormationHeating(EnergyChangeFunction):
    """Part of the H2 binding energy heats the gas (``type = "formation_heating"``)."""

    name = "formation_heating"
    reaction = "H.H__H2"

    @staticmethod
    def config_model() -> type[FormationHeatingParams]:
        return FormationHeatingParams

    def energy_change(self, reaction) -> Expr:
        params = cast(FormationHeatingParams, self.params(reaction))

        return Float(params.fraction * eV_to_erg(4.48))
```

```toml
[network.reactions.reaction."H.H__H2".thermodynamics]
energy_change = { type = "formation_heating", fraction = 0.33 }
```

A `deltaE<N>` for the reaction in the `.jfunc` file takes preference over the plugin.

## Example: a heating / cooling rate

In `plugins/physics/thermodynamics/line_cooling/`:

```python
from typing import cast

from sympy import Expr, exp

from jaff.physics.thermodynamics.heating_cooling_rate import (
    HeatingCoolingRateFunction,
    HeatingCoolingRateParams,
    register,
)


class LineCoolingParams(HeatingCoolingRateParams):
    """Settings of a two-level line cooling.

    Attributes
    ----------
    species : str
        The cooling species.
    coefficient : float
        Cooling coefficient (erg cm^3 s^-1).
    energy : float
        Upper-level energy (K).
    """

    species: str
    coefficient: float
    energy: float


@register
class LineCooling(HeatingCoolingRateFunction):
    """Two-level line cooling of one species (``type = "line_cooling"``)."""

    name = "line_cooling"

    @staticmethod
    def config_model() -> type[LineCoolingParams]:
        return LineCoolingParams

    def heating_cooling_rate(self, net, params) -> Expr:
        params = cast(LineCoolingParams, params)
        sym = net.symbols
        n = sym.ndens[net.species.from_name(params.species).index]

        return -params.coefficient * n * sym.ntot * exp(-params.energy / sym.tgas)
```

```toml
[network.thermodynamics]
heating_cooling_rate = [
    { type = "line_cooling", species = "CO", coefficient = 1e-27, energy = 5.5 },
    { type = "line_cooling", species = "O", coefficient = 1e-26, energy = 228.0 },
]
```

The rates of all entries are added to the `.jfunc` `heatingCoolingRate`.

## Example: emissivity and absorption spectra

In `plugins/physics/radiation/emissivity/` and
`plugins/physics/radiation/absorption_coefficient/`:

```python
from typing import cast

from sympy import Expr, Float

from jaff.physics.photo_reactions.absorption_coefficient import (
    AbsorptionCoefficientFunction,
    AbsorptionCoefficientParams,
)
from jaff.physics.photo_reactions.absorption_coefficient import (
    register as register_absorption,
)
from jaff.physics.photo_reactions.emissivity import EmissivityFunction, register


@register
class FlatEmission(EmissivityFunction):
    """A flat emission spectrum, 1e-25 erg cm^-3 s^-1 eV^-1 (``type = "flat"``)."""

    name = "flat"

    def emissivity(self, net, params) -> Expr:
        return Float(1e-25)


class DustExtinctionParams(AbsorptionCoefficientParams):
    """Settings of the dust extinction.

    Attributes
    ----------
    cross_section : float
        Dust cross section per H nucleus at 10 eV (cm^2).
    """

    cross_section: float = 1e-21


@register_absorption
class DustExtinction(AbsorptionCoefficientFunction):
    """Dust absorption rising linearly with photon energy (``type = "dust_extinction"``)."""

    name = "dust_extinction"

    @staticmethod
    def config_model() -> type[DustExtinctionParams]:
        return DustExtinctionParams

    def absorption_coefficient(self, net, params) -> Expr:
        params = cast(DustExtinctionParams, params)
        sym = net.symbols

        return params.cross_section * sym.n_hnuc * sym.ph_E / 10.0
```

```toml
[network.radiation]
bands = [6.0, 13.6]
emissivity = [{ type = "flat" }]
absorption_coefficient = [{ type = "dust_extinction" }]
```

The spectra are added to the `.jfunc` `globalEmissivity` /
`globalAbsorptionCoefficient` before the totals are integrated / averaged per band.

<!-- prettier-ignore -->
!!! warning "Dust module"
    With the dust module enabled, dust absorption is already applied to the radiation
    field: a dust absorption plugin on top of it absorbs twice.

## Settings models

A `...Params` subclass is a regular [config model](adding-configuration.md): strict
types, unknown keys rejected, and its numpydoc `Attributes` become the key
descriptions in the editor schemas. Give each plugin its own settings class. Plugins
registered under one `name` (the same type bound to several reactions) must share one
`config_model()`: a table is validated by its `type` before the reaction is known
(`register` raises otherwise). A plugin without settings keeps the default
`config_model()` and accepts no keys.

## Testing

Test a plugin through a `Network` that selects it. In JAFF's own tests, test-only
plugins never go into the real registries (they would leak into later tests, e.g. the
committed-schema comparison): `tests/helpers/energy_change_plugins.py`,
`tests/helpers/heating_cooling_rate_plugins.py` and `tests/helpers/radiation_plugins.py` install them
into a fresh registry per test with `install(monkeypatch)`.

## Checklist

- [x] Settings class (`...Params` subclass, numpydoc `Attributes`) returned by
      `config_model()`, or none
- [x] Plugin class subclasses the kind's base, sets `name` (and `reaction` for a
      per-reaction kind) and is decorated with the kind's `register`
- [x] Its method returns an expression in the network's symbols, in the kind's units
- [x] Placed in the kind's folder; every folder has an `__init__.py`, the last one
      imports the class; no folder on the way starts with `_`
- [x] Schemas regenerated (`uv run --extra dev python -m jaff.config.json_schema`) so
      editors complete the new type
- [x] `uv run --extra dev ty check plugins` is clean

## See Also

- [Plugins](../user-guide/configuration/plugins.md): selecting plugins and how they
  combine with the `.jfunc` file
- [Adding a Custom Shielding Function](adding-shielding-functions.md)
- [Auxiliary functions](../user-guide/designing-networks/auxiliary-functions.md)
- [Adding configuration](adding-configuration.md)
