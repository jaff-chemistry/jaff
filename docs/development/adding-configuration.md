---
tags:
    - Development
    - Configuration
---

# Adding Configuration

Every configurable part of JAFF owns its settings: a **config model** in a
`config.py` next to the code it configures. The models form a tree that mirrors
the TOML files: a parent model lists its children's models as fields. From that
tree JAFF validates `jaff.toml` / `jaffgen.toml`, generates the JSON Schemas
editors read and the [configuration reference](../user-guide/configuration/reference.md),
so a new setting is written **once**, in its model.

```text
src/jaff/
├── config/base.py                       # ConfigModel, LowerStr, PathLike, ...
├── core/network/config.py               # NetworkConfig ([network]), JaffTomlConfig
├── core/reaction/config.py              # ReactionsConfig, ReactionConfig, PhotoConfig
├── physics/photo_reactions/config.py    # RadiationConfig ([network.radiation])
├── physics/photo_reactions/shielding/config.py
├── physics/dust/config.py               # DustConfig, PhotoelectricEmissionConfig
├── physics/thermodynamics/config.py     # ThermodynamicsConfig
├── physics/thermodynamics/eos/config.py # the EOS models
└── cli/jaffgen/config.py                # JaffgenConfig ([jaffgen]), JaffgenTomlConfig
```

## Write the model

Subclass `ConfigModel` (`jaff.config.base`). It rejects unknown keys, is strict
about types (a TOML string is never turned into a number or a boolean) and is
frozen once validated. Document every key in a numpydoc `Attributes` section:
the entries become the key descriptions in the schema and the reference, so
write them for someone editing a TOML file.

```python
# src/jaff/physics/cosmic_rays/config.py
"""Config of the cosmic-ray module (``[network.cosmic_rays]``)."""

from __future__ import annotations

from typing import Annotated, Literal

from pydantic import BeforeValidator, Field

from ...config.base import ConfigModel, PathLike, lower

Spectrum = Annotated[Literal["power_law", "tabulated"], BeforeValidator(lower)]


class CosmicRaysConfig(ConfigModel):
    """Cosmic rays (``[network.cosmic_rays]``); the table enables the module.

    Attributes
    ----------
    ionisation_rate : float
        Primary ionisation rate per H2 in s^-1, > 0.
    spectrum : str
        Spectrum shape: ``"power_law"`` or ``"tabulated"`` (case-insensitive).
    table : Path or None
        Tabulated spectrum, used with ``spectrum = "tabulated"``.
    """

    ionisation_rate: float = Field(1.0e-17, gt=0.0)
    spectrum: Spectrum = "power_law"
    table: PathLike | None = None
```

| Building block | Use it for |
| --- | --- |
| `Annotated[Literal[...], BeforeValidator(lower)]` | A choice, case-insensitive and stored lower-case. |
| `LowerStr` | Free text compared case-insensitively. |
| `PathLike` | A path; a relative one is resolved against the file it is written in. |
| `Field(default, gt=..., ...)` | Bounds; they appear in the schema too. |
| `X \| None = None` | An optional key; TOML has no null, so *unset* means absent. |
| `Dict[str, Model]` | A table of named entries, like `reaction."<serialized>"`. |

## Attach it to its parent

Add the model as a field of the model of the enclosing table, e.g. in
`NetworkConfig`:

```python
cosmic_rays: CosmicRaysConfig | None = None   # the table enables the module
```

Use `| None = None` when the table's presence switches the module on (like
`radiation` and `dust`); use a default instance (`= CosmicRaysConfig()`) when
the module is always there (like `thermodynamics`). Document the field in the
parent's `Attributes` too. There is no registry: the parent's import of the
child is what puts it in the tree.

Then register the new module in two lists:

- `CONFIG_MODULES` in `tests/config/test_config_modules.py`, which imports every
  `config.py` in a fresh interpreter to catch import cycles;
- `CONFIG_MODULES` in `src/jaff/_utils/api_docs.py` (module path → nav title), which
  gives its models pages in the API reference. Give the module an `__all__`.

A `config.py` should import only pydantic, `jaff.config` and other `config.py`
modules.

## Use it

Every configured component follows the same pattern:

1. **It takes its own section**: `Thing(network, config)`, where `config` is the
   resolved model of its table, never the whole config.
2. **It keeps it as `self.config`**, read-only, and reads its settings from there.
3. **It passes `config.<child>` down** to the components it builds (as `Dust`
   does with `config.pe` for `PhotoelectricEmission`).
4. **Its other attributes are computed state** only: what it derives from the
   settings, never a copy of them, and it never reads another component's
   settings out of `net.config`.

```python
class CosmicRays:
    def __init__(self, network: Network, config: CosmicRaysConfig):
        self.network = network
        self.config = config  # read-only settings
        self.rate = self._build_rate()  # computed from them

    def _build_rate(self) -> Expr:
        return self.config.ionisation_rate * ...
```

`Network` builds it from its own config, as it does `Radiation`, `Dust` and
`Thermodynamics`:

```python
self.cosmic_rays = (
    CosmicRays(self, self.config.cosmic_rays)
    if self.config.cosmic_rays is not None
    else None
)
```

Per-reaction settings work the same way. The catalogue is
`Reactions(config.reactions)` and owns the `[network.reactions]` settings: every
`Reaction` holds `reaction.config`, its resolved `ReactionConfig`
(`net.reactions.config.resolve(serialized)`: the global table with its own
`reaction."<serialized>"` entry applied over it); `net.reactions.add_rate`
resolves a second rate over the same range by the reaction's
`duplicate_policy`; and `net.reactions.check_config()` reports per-reaction
entries that name no reaction. A photo reaction's `reaction.photo` is
itself a component, `PhotoReaction(reaction, reaction.config.photo)`: it keeps
the photo section as `photo.config`, computes the cross sections, band data and
shielding tables (`photo.compute()`), and offers `band_xsecs` and
`plot_xsecs()`. To tell a key set on
the reaction itself from the global value, use `net.reactions.config.own(serialized)`.

## Where each check belongs

| Check | Where | Example |
| --- | --- | --- |
| A value on its own: type, choice, bound | the field's type | `t_cutoff`, `gamma > 1` |
| Several keys of one table | a `field_validator` / `model_validator` raising `ValueError` | one `profile_index` per band |
| Against the network | where `Network` is built, raising `ParserError` with the key path (`jaff.config.base.dotted`) | a reaction key naming no reaction |

The loader turns validation errors into one `ParserError` listing every problem
by key path. Keep messages short and in TOML terms.

## Tables chosen by a `type`

When a table's keys depend on its `type` (the EOS, shielding), write one model
per type and select it by `type`:

- set `_atomic: ClassVar[bool] = True` on each model: a layer then replaces the
  table whole instead of merging it key by key with another type's keys;
- give each a `type: Annotated[Literal["..."], BeforeValidator(lower)]` field and
  combine them in an `Annotated[Union[...], Discriminator(...)]` (see
  `physics/thermodynamics/eos/config.py`);
- adjust a model's schema in `_extend_json_schema(cls, schema)` (the EOS models
  mark `type` required there).

For an open set of plugins, let each plugin declare its settings instead — the
shielding models return a `ShieldingParams` subclass from `config_model()` (see
[Adding a Custom Shielding Function](adding-shielding-functions.md)).

## `Network` keyword arguments and command-line options

`Network`'s keyword arguments are the keys of the `[network]` table, named and
nested as in the file (`reactions=`, `radiation=`, `dust=`, `thermodynamics=`,
...), and `_keyword_layer` passes them through unchanged. So:

- **a new top-level `[network]` table or key** (a new field of `NetworkConfig`)
  needs a `Network` parameter of the same name (default `None`, meaning not
  given), typed as the model or a `dict`, and documented in `__init__`;
  `tests/core/network/test_network_config.py` fails until it exists;
- **a key inside an existing table** needs nothing: it is reached through that
  table's keyword (`reactions={"photo": {...}}`), a `config=` layer or `--set`.

A command-line shortcut is optional: add the typer option and an entry in
`CLI_KEYS` (`cli/_helper.py`) mapping it to its key path.

## Regenerate and test

```bash
uv run --extra dev python -m jaff.config.json_schema   # editor schemas, src/jaff/config/schemas/
uv run --extra dev python -m jaff.config.reference     # docs/user-guide/configuration/reference.md (for local docs)
uv run --extra dev pytest tests/config
```

`tests/config/test_config_json_schema.py` fails when a committed schema is out of date
with the models. The reference page is not committed: CI generates it before
building the docs.

## Checklist

- [x] `config.py` next to the module, with a `ConfigModel` subclass
- [x] Every key documented in numpydoc `Attributes`, worded for TOML users
- [x] Field added to the parent model (and documented there)
- [x] Module has `__all__` and is listed in both `CONFIG_MODULES`
      (`tests/config/test_config_modules.py`, `src/jaff/_utils/api_docs.py`)
- [x] The consuming class takes `(network, config)`, keeps it as `self.config`
      and stores only computed state; `Network` builds it
- [x] A new top-level `[network]` key has a same-named `Network` keyword
- [x] Checks placed by kind (type / validator / network build)
- [x] Type-selected tables: `_atomic`, one model per type, discriminated union
- [x] Schemas and reference regenerated; tests for defaults, bad values and
      unknown keys

## See Also

- [Configuration overview](../user-guide/configuration/index.md)
- [Adding a Custom Shielding Function](adding-shielding-functions.md)
- [Writing API Documentation](api-documentation.md) — how the `Attributes` text is rendered
- [Code Style Guide](code-style.md)
