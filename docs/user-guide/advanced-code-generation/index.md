---
tags:
    - User-guide
    - Code-generation
icon: phosphor/cpu
---

# Advanced Code Generation

JAFF exposes two programmatic interfaces for code generation that bypass the template-file workflow entirely. Both accept a parsed `Network` and produce target-language code directly in Python.

| Interface               | Use when                                                                                     |
| ----------------------- | -------------------------------------------------------------------------------------------- |
| [`Builder`](builder.md) | You want a complete, runnable solver project in one call using a built-in template plugin    |
| [`Codegen`](codegen.md) | You need fine-grained control over individual output strings: rates, fluxes, ODEs, Jacobians |

---

<div class="grid cards" markdown>

- :phosphor-cube:{ .sm .middle } **Builder**

    ***

    One-call code generation using named plugin templates. Generates complete solver projects (Python `solve_ivp`, Fortran DLSODES, C++ VODE via the `kokkos_ode` template, AMReX Microphysics) from a loaded network.

    [:octicons-arrow-right-24: Builder](builder.md)

- :phosphor-gear:{ .sm .middle } **Codegen**

    ***

    Low-level engine for generating rates, fluxes, ODEs, RHS vectors, Jacobians, and radiation ODEs as formatted strings for any of the seven supported languages.

    [:octicons-arrow-right-24: Codegen](codegen.md)

</div>

---

## Quick Comparison

```python
from jaff import Network
net = Network("networks/GOW/GOW.jet")

# Builder — complete project in one call
from jaff.codegen import Builder
b = Builder(net)
b.build(template="python_solve_ivp", output_dir="output/")

# Codegen — individual expression strings
from jaff.codegen import Codegen
cg = Codegen(net, lang="python")
print(cg.get_rates_str())
print(cg.get_ode_str())
```
