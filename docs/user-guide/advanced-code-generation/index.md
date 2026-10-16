---
tags:
    - User-guide
    - Code-generation
icon: phosphor/cpu
---

# Advanced Code Generation

JAFF exposes a programmatic interface for code generation that bypasses the template-file workflow entirely. It accepts a parsed `Network` and produces target-language code directly in Python.

| Interface               | Use when                                                                                     |
| ----------------------- | -------------------------------------------------------------------------------------------- |
| [`Codegen`](codegen.md) | You need fine-grained control over individual output strings: rates, fluxes, ODEs, Jacobians |

---

<div class="grid cards" markdown>

- :phosphor-gear:{ .sm .middle } **Codegen**

    ***

    Low-level engine for generating rates, fluxes, ODEs, RHS vectors, Jacobians, and radiation ODEs as formatted strings for any of the seven supported languages.

    [:octicons-arrow-right-24: Codegen](codegen.md)

</div>

---

## Quick Example

```python
from jaff import Network
net = Network("networks/GOW/GOW.jet")

# Codegen — individual expression strings
from jaff.codegen import Codegen
cg = Codegen(net, lang="python")
print(cg.get_rates_str())
print(cg.get_ode_str())
```
