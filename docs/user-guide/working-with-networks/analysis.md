---
tags:
    - User-guide
    - Network
---

# Analysis mode

A network's rate expressions are symbolic: they hold parameters (`tgas`, `av`,
`crate`, column densities, ...), species number densities (`nden[i]`) and
sometimes unknown functions (functions the `.jfunc` file calls but does not
define, such as GOW's CO shielding table `fshield_co_interp2d`). Analysis mode,
`net.analysis`, holds numeric values for all of them in one place. Two things
read them:

- the evaluate methods of reactions, which give rate coefficients *k* and rates
  𝓡 as numbers ([Evaluating rate coefficients and rates](reactions.md#evaluating-rate-coefficients-and-rates));
- the network figures, which colour reactions by *k* or 𝓡 and list the values
  in a legend ([Rate colouring](visualization.md#rate-colouring)).

```python
from jaff import Network

net = Network("GOW")
net.analysis.enable(
    parameters={"tgas": 20.0, "av": {"value": 10.0, "label": r"$A_V$"}, "chi": 1.0, "crate": 2e-16},
    densities={"H": 1.0, "H2": {"value": 100.0, "label": True}},
    functions={"fshield_co_interp2d": 7.3e-5},
)
net.analysis.parameters.values()
# {'av': 10.0, 'chi': 1.0, 'crate': 2e-16, 'tgas': 20.0}
```

---

## Enabling and disabling

`net.analysis` exists on every network but is off until enabled
(`net.analysis.enabled` is `False` and its tables are `None`).

- `enable(parameters=..., densities=..., functions=...)` builds the tables from
  the network when off, then sets the given entries. Calling it again merges
  the new values over the stored ones.
- `disable()` drops the tables and their values.
- `pop()` clears every entry of every table; analysis mode stays on.
- `get_spec()` returns a read-only snapshot of the three tables (an
  `AnalysisSpec`), or `None` when off.

Values live only on the network object: they are not part of `net.config` and
are not saved.

---

## The tables

`enable()` builds three tables, each with a fixed set of names taken from the
network:

| table        | names                                              | a value is                                                                                   |
| ------------ | -------------------------------------------------- | -------------------------------------------------------------------------------------------- |
| `parameters` | the symbols of the rate coefficients (not `nden`)  | a finite number                                                                              |
| `densities`  | every species, by name                             | a finite number ≥ 0, in cm⁻³; `enable` also takes an array in species order                  |
| `functions`  | the unknown functions of the rate coefficients     | a number, a sympy expression of the parameters, or a callable without arguments returning either |

```python
net.analysis.parameters.names
# ('av', 'chi', 'crate', 'd2g', 'ncol_C', 'ncol_CO', 'ncol_H2', 'tgas', 'vdisp')
net.analysis.functions.names
# ('fshield_co_interp2d',)
```

A function takes one value wherever it is called: its arguments in the network
are not used. An expression may use only parameters, and neither an expression
nor a callable's result may contain an unknown function.

Each table has the same methods:

| method                          | does                                                                                     |
| ------------------------------- | ---------------------------------------------------------------------------------------- |
| `set(name, value, label=None)`  | set a value (and a label, see below)                                                     |
| `update(mapping)`               | set several: per name a value, or `{"value": ..., "label": ...}`                         |
| `get(name)` / `get()`           | the `Entry(value, label)` of a name / a copy of every entry                              |
| `values()`                      | the value of every name that has one                                                     |
| `pop(name)` / `pop()`           | clear one entry and return the old one / clear the whole table                           |
| `missing(names)`                | the names without a value                                                                |

```python
net.analysis.parameters.get("av")  # Entry(value=10.0, label='$A_V$')
net.analysis.parameters.set("tgass", 1.0)
# KeyError: "'tgass' is not a parameter of the network; did you mean tgas?"
net.analysis.parameters.pop("tgas")  # Entry(value=20.0, label=None)
```

A name the network does not use raises `KeyError`, with close matches; a value
that is not allowed raises `ValueError`.

---

## Labels

Each entry has an optional label. Labelled entries are listed in the values
legend of rate-coloured figures; unlabelled ones are used but not shown.

- A `{}` in the label takes the value: `r"$T_\mathrm{gas} = {}$ K"` reads
  *T*<sub>gas</sub> = 20 K. Inside maths the value goes in as it is, outside it in
  its own maths.
- A label without `{}` reads `label = value`: `r"$A_V$"` gives *A*<sub>V</sub> = 10.
- `label=True` gives the table's default: `$tgas = {}$` for a parameter,
  `$n(<species LaTeX>) = {}$` for a density, `name = {}` for a function.
- `label=False` removes the label; leaving `label` out keeps the current one.

Values show 3 significant figures: plain from 0.01 to 1000, else `m × 10ⁿ`.

```python
net.analysis.densities.get("H2")  # Entry(value=100.0, label='$n({\\rm H_{2}}) = {}$')
```

---

## Evaluation

The evaluate methods (`evaluate_rate_coefficient()`, `evaluate_rate()` and the
bulk `net.reactions.evaluate_rate_coefficients()` / `evaluate_rates()`) take
`parameters`, `densities` and `functions` too. With analysis mode on, each
argument is merged over the matching table, so a call needs only what differs;
with it off, only the given values are used. Missing values raise one
`ValueError` naming all of them:

```python
rxn = net.reactions[24]  # CO + H3+ -> H2 + HCO+
rxn.evaluate_rate_coefficient()  # 1.7e-09
rxn.evaluate_rate()
# ValueError: Cannot evaluate: no density for CO, H3+
```
