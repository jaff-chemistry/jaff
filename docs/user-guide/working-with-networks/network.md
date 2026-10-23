---
tags:
    - User-guide
    - Network
---

# Network

`Network` is the entry point for every JAFF workflow and the hub the other
objects hang off. One constructor call reads a network file, parses it into
typed [`Species`](species.md), [`Reactions`](reactions.md), and
[`Elements`](elements.md) catalogues, rewrites every rate into a common symbolic
form, validates the chemistry, and exposes the symbolic ODE/flux/radiation
expressions that downstream code generation consumes.

```python
from jaff import Network

net = Network("networks/h_photoionization/h_photo.jet")
```

---

## What loading does

A single `Network(...)` call runs a fixed pipeline. Knowing the phases explains
where every attribute comes from:

1. **Read** — the file is resolved and, if it is a binary `.jaff`, deserialized
   directly; otherwise the text parser auto-detects the format (KROME, PRIZMO,
   UDFA, KIDA, UCLChem, …).
2. **Parse & assemble** — reactions are read one by one. Each new species name
   becomes a [`Specie`](species.md) appended to `net.species`; each line becomes
   a [`Reaction`](reactions.md) appended to `net.reactions`. Any accompanying
   `.jfunc` auxiliary file is loaded here and its custom rates, heating/cooling,
   and radiation terms are wired in.
3. **Standardize symbols** — every rate and energy expression is rewritten so
   that number densities use the common [`nden[i]`](#number-densities-the-nden-representation)
   representation (shorthands like `n_X`, `n_X_nuc`, `ntot` are expanded).
4. **Validate** — mass/charge conservation, sinks/sources, recombinations,
   isomers, and duplicate reactions are [checked](#validation).
5. **Finalize** — integer [stoichiometry matrices](#key-attributes) are built
   and the unique [`Elements`](elements.md) set is derived from the species.

After this, the network is a fully assembled, queryable model.

---

## Constructor

```python
Network(
    fname,
    config=None,
    strict=None,
    label=None,
    funcfile=None,
    duplicate_policy=None,
    expand_nuclei=None,
    radiation_props=None,
    dust_props=None,
    use_proxy_photoreaction=None,
    eos_props=None,
    use_jaff_toml=True,
)
```

Settings come from the network's [configuration](../configuration/index.md):
the `jaff.toml` next to the network file, then the `config` layers, then the
keyword arguments. Every keyword argument defaults to `None` — *not given* — so
it only overrides the config when you pass it. The validated result is
`net.config`.

| Parameter                 | Type                                 | Description                                                                                                                         |
| ------------------------- | ------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------- |
| `fname`                   | `str or Path`                        | Path to the network file (required), or a predefined network name; `.jaff` files are loaded as binary                              |
| `config`                  | path, `dict`, `JaffTomlConfig`, or a list of them | Config layers over the network's `jaff.toml`; a `dict` is shaped like a `jaff.toml` (`{"network": {...}}`)            |
| `strict`                  | `bool`                               | Treat validation problems (conservation, sinks/sources, ...) as fatal (`network.strict`; default `False`)                          |
| `label`                   | `str`                                | Human-readable network name (`network.label`; default the file stem)                                                               |
| `funcfile`                | `bool or str or Path`                | `.jfunc` auxiliary file; `True` scans the network dir, `False` skips loading (`network.funcfile`; default `True`)                  |
| `duplicate_policy`        | `str`                                | `preserve-first`, `preserve-last` or `error` (`network.reactions.duplicate_policy`; default `preserve-first`)                      |
| `expand_nuclei`           | `bool`                               | Expand `n_<element>_nuc` into the element-nucleus sum (`network.expand_nuclei`; default `True`)                                    |
| `radiation_props`         | `RadiationConfig` or `dict`          | Radiation field (`network.radiation`); enables radiation transport and photochemistry                                              |
| `dust_props`              | `DustConfig` or `dict`               | Dust (`network.dust`); enables the dust module                                                                                      |
| `use_proxy_photoreaction` | `bool`                               | Look up cross-sections under the proxy key (`CHx`→`CH`, `OHx`→`OH`) (`network.reactions.photo.use_proxy_photoreaction`; default `True`) |
| `eos_props`               | EOS config model or `dict`           | Equation of state (`network.thermodynamics.eos`); replaces a configured EOS whole                                                   |
| `use_jaff_toml`           | `bool`                               | Stack the `jaff.toml` next to the network file under the `config` layers (default `True`)                                          |

### Examples

```python
from jaff import Network
from jaff.physics.photo_reactions.config import RadiationConfig

# Minimal load (with the network's jaff.toml, if any)
net = Network("networks/GOW/GOW.jet")

# Strict validation — abort on the first conservation failure
net = Network("networks/GOW/GOW.jet", strict=True)

# Custom label and explicit auxiliary-function file
net = Network(
    "networks/GOW/GOW.jet",
    label="GOW-2017",
    funcfile="networks/GOW/GOW.jfunc",
)

# Config layers: a file, then a dict
net = Network(
    "networks/GOW/GOW.jet",
    config=["tuned.toml", {"network": {"reactions": {"t_cutoff": "clip"}}}],
)

# With radiation bands (needed for photochemistry during code generation)
net = Network(
    "networks/h_photoionization/h_photo.jet",
    radiation_props=RadiationConfig(bands=[13.6, "inf"], mode="nph"),
)
```

---

## Key attributes

| Attribute         | Type              | Description                                                                                                         |
| ----------------- | ----------------- | ------------------------------------------------------------------------------------------------------------------- |
| `filename`        | `Path`            | Absolute path to the source file                                                                                    |
| `label`           | `str`             | Network name                                                                                                        |
| `config`          | `NetworkConfig`   | The validated [configuration](../configuration/index.md) the network was built from                                |
| `spec`            | `NetworkSpec`     | Normalized construction params (resolved `fname`, `aux_funcs`)                                                      |
| `species`         | `Species`         | Ordered [`Species`](species.md) catalogue                                                                           |
| `reactions`       | `Reactions`       | Ordered [`Reactions`](reactions.md) catalogue                                                                       |
| `elements`        | `Elements`        | [`Elements`](elements.md) catalogue derived from the species                                                        |
| `reactant_matrix` | `np.ndarray`      | Integer stoichiometry, shape `(n_reactions, n_species)` — reactant counts                                           |
| `product_matrix`  | `np.ndarray`      | Integer stoichiometry, same shape — product counts                                                                  |
| `thermodynamics`  | `Thermodynamics`  | Thermal equations: `eos`, `dEdt_chemical`, `dEdt_extra`, `dEdt_tot`, `dTdt_chemical`, `dTdt_extra`, `dTdt_tot`; see [thermodynamics](../../reference/physics/thermodynamics/index.md) |
| `radiation`       | `Radiation\|None` | Radiation field object; `None` when the config has no `[network.radiation]` table                                   |
| `mass_dict`       | `dict`            | Element mass dictionary used to build the species                                                                   |
| `symbols`         | `NetworkSymbols`  | Canonical symbols, densities (`ndens`, `ntot`, `rho`, `n_hnuc`), introspection sets and `standardize`; see [symbols](../../reference/core/networksymbols/index.md) |

```python
net = Network("networks/h_photoionization/h_photo.jet")

net.label                   # 'h_photo'
net.species.count           # 3
net.reactions.count         # 2
net.reactant_matrix.shape   # (2, 3)   — 2 reactions × 3 species
```

The three catalogues are the same typed collections
documented on [Species](species.md), [Reactions](reactions.md) and
[Elements](elements.md) — anything there works straight off the network:

```python
net.species.names()                 # ['H', 'H+', 'e-']
net.species.charged("name")         # ['H+', 'e-']
net.reactions.verbatim()            # ['H -> H+ + e-', 'H+ + e- -> H']
net.reactions.photo_reactions()     # [ReactionObject(H -> H+ + e-)]
net.elements.symbols()              # ['H']
```

### Stoichiometry matrices

`reactant_matrix[i, j]` counts how many times species _j_ appears as a reactant
in reaction _i_ (and `product_matrix` likewise for products). Row order is
reaction index; column order is species index.

```python
net.reactant_matrix
# [[1 0 0]      reaction 0: H -> H+ + e-   consumes 1× H (species 0)
#  [0 1 1]]     reaction 1: H+ + e- -> H   consumes 1× H+ and 1× e-

net.product_matrix
# [[0 1 1]
#  [1 0 0]]
```

The net stoichiometry `product_matrix - reactant_matrix` is the change in each
species per reaction — the backbone of the ODE system below.

---

## Number densities: the `nden` representation

Inside JAFF every species number density is a reference into one SymPy
`MatrixSymbol("nden", n_species, 1)`. Species _X_ is `nden[idx_X, 0]`, where the
index is the species' position in `net.species`. This is why rate and ODE
expressions print terms like `nden[1, 0]` rather than a named density — it keeps
the whole network on one indexable vector the code generator can emit as an
array.

During loading (phase 3), convenience shorthands in rate expressions are
expanded into this form:

| Shorthand                                      | Expands to                                                   |
| ---------------------------------------------- | ------------------------------------------------------------ |
| `ntot`                                         | sum of `nden[i, 0]` over **all core** species                |
| `n_<element>_nuc` (e.g. `n_H_nuc`, `n_He_nuc`) | sum over element-bearing species, **weighted** by atom count |
| `n_e`                                          | electron density                                             |
| `n_<species>` (e.g. `n_CO`, `n_Hj`, `n_Hk`)    | `nden[idx_X, 0]` for that one species (`j`→`+`, `k`→`-`)     |

For example, a cosmic-ray ionization rate written with the `n_H_nuc` shorthand
comes out as a weighted `nden` sum (note `2*nden[3, 0]` — H₂ contributes two H
atoms):

```python
g = Network("networks/GOW/GOW.jet")
g.reactions.with_type("cosmic_ray")[0].rate_coefficient
# crate*(1.5*nden[0, 0]/(nden[0, 0] + nden[1, 0] + 2*nden[3, 0] + ...) + ...)
```

Pass `expand_nuclei=False` to keep `n_<element>_nuc` tokens as opaque free
symbols (e.g. `nh_nuc`) instead of expanding them — useful when an external
driver supplies those densities directly.

---

## Symbolic expressions

Once loaded, the network exposes its dynamics as SymPy expressions.

### Species ODEs

`Network.sodes()` returns one expression per species — the net rate of change of
its number density (cm⁻³ s⁻¹), summed over every reaction it takes part in:

<!-- prettier-ignore -->
$$ \frac{dn}{dt} = \sum_{i} (-1)^a f_i \quad \begin{cases} a = 0 & n \text{ is a product of reaction } i \\ a = 1 & n \text{ is a reactant of reaction } i \end{cases} $$

where $f_i$ is the flux of the $i^{\text{th}}$ reaction.

```python
for specie, expr in zip(net.species, net.sodes()):
    print(f"d[{specie}]/dt = {expr}")

# d[H]/dt  = 1.659e-10*nden[1, 0]*nden[2, 0]/tgas**0.7 - photorates(1, 13.6, 1.0e+99)*nden[0, 0]
# d[H+]/dt = -1.659e-10*nden[1, 0]*nden[2, 0]/tgas**0.7 + photorates(1, 13.6, 1.0e+99)*nden[0, 0]
# d[e-]/dt = -1.659e-10*nden[1, 0]*nden[2, 0]/tgas**0.7 + photorates(1, 13.6, 1.0e+99)*nden[0, 0]
```

### Fluxes and radiation moments

```python
net.sfluxes()       # list[Expr] — per-reaction flux  k_i * nden[r1] * nden[r2] ...
net.sradodes(0)     # list[Expr] — radiation moment ODEs (order 0), if radiation configured
```

A flux is just the reaction rate times its reactant densities; the species ODEs
are signed sums of these. See [`sfluxes`](../../reference/core/network/sfluxes.md) and
[`sradodes`](../../reference/core/network/sradodes.md) for details.

---

## Temperature cutoffs

Each reaction's rate is only valid within a temperature range `[Tmin, Tmax]`.
Whether `tgas` is clamped to that range (`clip`, the default) or left to
extrapolate beyond it is the `t_cutoff` setting, for every reaction in
`[network.reactions]` or for one reaction in its own table — see the
[configuration reference](../configuration/reference.md#networkreactions).

---

## Validation

JAFF validates the network during loading. With `strict=False` (default) each
problem is logged as a warning; with `strict=True` the process exits on the
first one. The individual checks can also be re-run by hand — each takes an
`errors` flag (exit when it finds a problem), logs its findings and returns a
report whose `.ok` is `True` when nothing was found:

```python
net.check_sink_sources(errors=False).sinks     # species never produced
net.check_recombinations(errors=False).missing  # cations lacking an electron recombination
net.check_isomers(errors=False).groups          # species sharing one atomic composition
```

Duplicate reactions (same species, mechanism/`type`, and temperature range) are
resolved at load time by the network's
[`duplicate_policy`](../configuration/reference.md#networkreactions) rather than by a separate
check.

<!-- prettier-ignore -->
!!! note "The checks report, they don't return"
    These methods write warnings through the logger and return `None`. To act on
    the results programmatically, inspect the catalogues directly (e.g. compare
    `reactant_matrix` / `product_matrix` columns for sinks and sources) rather
    than expecting a returned list.

Mass and charge conservation are validated per reaction at construction time —
see [Reactions](reactions.md#conservation-is-checked-at-construction).

---

## Comparing networks

`compare_reactions` and `compare_species` diff two loaded networks, logging what
is shared and what is unique to each. `verbosity=1` (default) also prints the
full sets:

```python
net2 = Network("networks/GOW/GOW.jet")

net.compare_reactions(net2)             # reactions in one but not the other
net.compare_species(net2, verbosity=0)  # counts only, no per-item listing
```

Comparison is by [serialized form](reactions.md#reaction-identity-two-serialized-forms),
so isomer/duplicate handling is consistent with the rest of JAFF.

---

## Export and caching

### Rate tables

Reaction rates can be tabulated against temperature and written to HDF5 or text
for use by external solvers:

```python
net.to_hdf5("rates.hdf5", T_min=10, T_max=1e4, nT=200, err_tol=1e-3)
net.to_txt("rates.txt", T_min=10, T_max=1e4, nT=200)
net.write_table("rates.hdf5", format="auto")  # general form; format from extension (.hdf5/.hdf/.txt)
```

`err_tol` drives adaptive temperature sampling (set it to `None` for a fixed
grid). See [`to_hdf5`](../../reference/core/network/to_hdf5.md) and
[`to_txt`](../../reference/core/network/to_txt.md).

### JAFF binary format

`to_jaff` serializes an already-loaded network to a gzip-compressed JSON file
(conventional extension `.jaff`). Re-loading skips parsing and the expensive
SymPy assembly, so it is the recommended cache for large networks:

```python
net = Network("networks/GOW/GOW.jet")
net.to_jaff("gow.jaff")        # save

net2 = Network("gow.jaff")     # reload — same species, reactions, and ODEs
```

<!-- prettier-ignore -->
!!! warning "Networks with undefined functions can't be serialized"
    If any rate contains an undefined SymPy function — most commonly an
    unresolved `photorates(...)` (a photo-reaction loaded without radiation)
    or a custom `interp(...)` — `to_jaff` raises `ValueError`, because such
    expressions cannot be round-tripped through the JSON schema. Resolve the
    radiation/interpolation first, or keep those networks in their text form.
