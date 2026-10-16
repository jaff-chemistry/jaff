---
tags:
    - Api
    - Code-generation
icon: phosphor/brackets-curly
---

# jaff.codegen

The `jaff.codegen` subpackage contains the classes for generating solver source code from a chemical reaction network.

## Classes

| Class                             | Description                                                       |
| --------------------------------- | ----------------------------------------------------------------- |
| [`Codegen`](codegen/index.md)           | Low-level multi-language code generator (rates, ODEs, Jacobian)   |
| [`Preprocessor`](preprocessor/index.md) | Pragma-based template preprocessor for file-level code generation |

## Quick Start

```python
from jaff import Network
from jaff.codegen import Codegen

net = Network("networks/react_COthin")

# Generate individual code blocks
cg = Codegen(network=net, lang="c++")
rates = cg.get_rates_str(use_cse=True)
odes  = cg.get_ode_str(use_cse=True)
jac   = cg.get_jacobian_str(use_cse=True)
```
