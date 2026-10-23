# ABOUTME: Code generation from a reaction network: Codegen, templates, preprocessor
# ABOUTME: and the plugin Builder; re-exports the public classes and typed dicts.
"""Generate solver source code from a chemical reaction network.

``jaff.codegen`` turns a loaded :class:`~jaff.core.network.Network` into
source code for rate coefficients, fluxes, ODE right-hand sides and
Jacobians, in C, C++, Fortran, Python, Rust, Julia or R. It offers three
levels of entry point:

- :class:`Codegen` is the low-level generator: each ``get_*_str`` method
  returns one block of code as a string and each ``get_indexed_*`` method
  returns the expressions as an :class:`IndexedReturn`. :class:`Language`
  holds the per-language syntax it prints with.
- :class:`TemplateParser` renders a ``$$JAFF`` template, an ordinary source
  file with ``$$JAFF`` directives that expand to generated code; it is the
  engine behind the ``jaffgen`` command (template syntax is described in the
  user guide's *Template Syntax* page).
- :class:`Builder` builds a complete solver project from a named plugin
  template (``python_solve_ivp``, ``fortran_dlsodes``, ``kokkos_ode`` or
  ``microphysics``); the plugins fill ``PREPROCESS_<KEY>`` marker blocks of
  their template files with :class:`Preprocessor`.

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
'const double cse0 = nden[24]*nden[9];'

Build a ready-to-run ``scipy.integrate.solve_ivp`` solver into ``output/``:

>>> from jaff.codegen import Builder
>>> path = Builder(net).build(template="python_solve_ivp", output_dir="output/")
Building network with template: python_solve_ivp
Network built successfully using template 'python_solve_ivp'.
Output files are located in: output/

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
from ._typing import ExtrasDict, IndexedReturn
from .builder import Builder
from .codegen import Codegen
from .preprocessor import Preprocessor
from .templating import TemplateParser

__all__ = [
    "Builder",
    "Codegen",
    "Language",
    "IndexedReturn",
    "Preprocessor",
    "TemplateParser",
    "ExtrasDict",
]
