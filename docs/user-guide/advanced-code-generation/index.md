---
tags:
    - User-guide
    - Code-generation
icon: phosphor/cpu
---

# Advanced Code Generation

JAFF exposes a programmatic interface for code generation that bypasses the
template-file workflow: [`Codegen`](codegen.md) takes a parsed `Network` and
returns target-language code directly in Python, as strings for rates, fluxes,
ODEs, Jacobians and more. Use it when you need fine-grained control over
individual output blocks.

For a complete, runnable solver project, use the built-in templates through
[`jaffgen --template <name>`](../code-generation/jaffgen.md) instead
(`python_solve_ivp`, `fortran_dlsodes`, `kokkos_ode`, `microphysics`, ...).

---

<div class="grid cards" markdown>

- :phosphor-gear:{ .sm .middle } **Codegen**

    ***

    Low-level engine for generating rates, fluxes, ODEs, RHS vectors, Jacobians, and radiation ODEs as formatted strings for any of the seven supported languages.

    [:octicons-arrow-right-24: Codegen](codegen.md)

</div>

---

## Example

```python
from jaff import Network
from jaff.codegen import Codegen

net = Network("networks/GOW/GOW.jet")
cg = Codegen(net, lang="python")
print(cg.get_rates_str())
print(cg.get_ode_str())
```
