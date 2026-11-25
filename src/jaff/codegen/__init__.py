# ABOUTME: Code generation from a reaction network: Codegen and $$JAFF templates;
# ABOUTME: re-exports the public classes and typed dicts.
"""Generate solver source code from a chemical reaction network.

``jaff.codegen`` turns a loaded :class:`~jaff.core.network.Network` into
source code for rate coefficients, fluxes, ODE right-hand sides and
Jacobians, in C, C++, Fortran, Python, Rust, Julia or R. It offers two
levels of entry point:

- :class:`Codegen` is the low-level generator: each ``get_*_str`` method
  returns one block of code as a string and each ``get_indexed_*`` method
  returns the expressions as an :class:`IndexedReturn`. :class:`Language`
  holds the per-language syntax it prints with.
- :class:`TemplateParser` renders a ``$$JAFF`` template, an ordinary source
  file with ``$$JAFF`` directives that expand to generated code; it is the
  engine behind the ``jaffgen`` command (template syntax is described in the
  user guide's *Template Syntax* page). The built-in templates
  (``python_solve_ivp``, ``fortran_dlsodes``, ``kokkos_ode``, ``microphysics``
  ...) live in ``templates/<name>/`` at the repository root and are rendered with
  ``jaffgen --template <name>``.

Examples
--------
Generate individual code blocks with :class:`Codegen`:

>>> from jaff import Network
>>> from jaff.codegen import Codegen
>>> net = Network("COthin")
>>> cg = Codegen(network=net, lang="c++")
>>> rates = cg.get_rates_str(use_cse=True)
>>> odes = cg.get_ode_str(use_cse=True)
>>> jac = cg.get_jacobian_str(use_cse=True)
>>> odes.splitlines()[0]
'const double cse0 = std::log(tgas);'

Render a ``$$JAFF`` template. Given ``nspec.cpp``:

```cpp
$$JAFF SUB nspec
const int n = $nspec$;
$$JAFF END
```

the directives expand against the network and are echoed as comments:

>>> from pathlib import Path
>>> from jaff.codegen import TemplateParser
>>> print(TemplateParser(net, Path("nspec.cpp")).parse_file(), end="")
// $$JAFF SUB nspec
const int n = 37;
// $$JAFF END
"""

from ._languages import Language
from ._typing import CseGroup, ExtrasDict, IndexedReturn
from .codegen import Codegen
from .templating import TemplateFunctions, TemplateParser, jaffgen_function

__all__ = [
    "Codegen",
    "Language",
    "IndexedReturn",
    "CseGroup",
    "TemplateParser",
    "TemplateFunctions",
    "jaffgen_function",
    "ExtrasDict",
]
