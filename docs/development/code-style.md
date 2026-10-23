---
tags:
    - Development
icon: phosphor/chart-bar
---

# Code Style

## Python Version

JAFF requires Python 3.11 and higher (`requires-python = ">=3.11"` in `pyproject.toml`).

Start every module with `from __future__ import annotations`, then use **built-in
generics** (`list`, `dict`, `tuple`) and **PEP 604 unions** (`X | None`) — not the
legacy `typing` aliases.

```python
# Good - built-in generics and PEP 604 unions
from __future__ import annotations

def process(items: list[str]) -> dict[str, int]:
    pass

def lookup(key: str) -> int | None:
    pass

# Avoid - legacy typing aliases
from typing import List, Dict, Optional

def process(items: List[str]) -> Optional[Dict[str, int]]:
    pass
```

## Code Formatting

Use Ruff for formatting and linting, and `ty` for type checking; all three are in
the `dev` extra. The Ruff config lives in `pyproject.toml` under `[tool.ruff]`:

- **Line length:** 90 characters
- **Quote style:** double quotes
- **Indent:** 4 spaces

The tools read their config automatically — no need to pass flags:

```bash
# Check formatting without modifying
uv run --extra dev ruff format --check src/jaff/core/network/network.py

# Format a file you changed
uv run --extra dev ruff format src/jaff/core/network/network.py

# Lint
uv run --extra dev ruff check src/jaff/core/network/network.py

# Type check (ruff does not catch type errors)
uv run --extra dev ty check src/jaff/core/network/network.py
```

Format and lint **only the files you touched**, and only reformat a file that was
already formatted before your change (check with `ruff format --check` first).
`tests/` is not lint-clean, and a repo-wide `ruff format` or `ruff check --fix`
rewrites unrelated files. `ty` reports some pre-existing errors; make sure your
change adds none.

Ruff's lint selection is `E4`, `E7`, `E9` and `F` (pyflakes): it does not sort
imports. Order them by hand: `from __future__ import annotations`, then the standard
library, third-party packages and local (`jaff`) imports, each group separated by a
blank line.

## File Headers

Each code file starts with a two-line `ABOUTME` comment, before the module
docstring:

```python
# ABOUTME: <short description>
# ABOUTME: <optional second line>
```

Not every existing module has one yet and nothing checks it; add the header to new
files and to files you edit.

## Naming Conventions

Follow PEP 8 naming throughout.

### Variables and Functions

Use `snake_case`:

```python
# Good
user_name = "Alice"
def calculate_rate(temperature):
    pass

# Avoid
userName = "Alice"
def CalculateRate(temperature):
    pass
```

### Classes

Use `PascalCase`:

```python
# Good
class Network:
    pass

class CodeGenerator:
    pass

# Avoid
class network:
    pass
```

### Constants

Use `UPPER_SNAKE_CASE`:

```python
# Good
MAX_ITERATIONS = 1000
DEFAULT_TEMPERATURE = 300.0

# Avoid
maxIterations = 1000
```

### Private Members

Prefix with underscore:

```python
class Network:
    def __init__(self):
        self._internal_data = []  # Private
        self.public_data = []     # Public

    def _helper_method(self):  # Private
        pass

    def public_method(self):   # Public
        pass
```

## Type Hints

Annotate every public function — both parameters and return type.

### Always Use Type Hints

Untyped signatures hide intent and disable static checking. Always declare types.

```python
# Good
def compute_rate(temperature: float, alpha: float) -> float:
    return alpha * temperature

# Avoid
def compute_rate(temperature, alpha):
    return alpha * temperature
```

### Import Types

Use built-in generics and `X | None`. Only import from `typing` for names that
have no built-in form (e.g. `Any`, `TypedDict`, `cast`, `TYPE_CHECKING`).
Older modules still use `List`/`Dict`/`Optional`; convert them when you touch them.

```python
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

def process_species(
    names: list[str],
    masses: np.ndarray,
    options: dict[str, Any] | None = None,
) -> tuple[list[str], np.ndarray]:
    pass
```

### Complex Types

Name recurring or compound types with aliases to keep signatures readable.

```python
from __future__ import annotations

# Type aliases for clarity
SpeciesDict = dict[str, int]
RateExpression = str | float

def load_network(
    filename: str,
    species_map: SpeciesDict | None = None,
) -> Network:
    pass
```

## Docstrings

Every public module, class, and function gets a docstring.

### NumPy Style

Use **NumPy-style** (numpydoc) docstrings with underlined sections
(`Parameters`, `Returns`, `Raises`, `Notes`, `Examples`). Cross-reference other
objects with Sphinx roles — `:class:`, `:meth:`, `:func:` — so the rendered API
docs link correctly; a leading `~` (`:class:`~jaff.core.network.Network``) shows only
the last name. Write inline code as ``` ``x`` ```, use fenced code blocks (never a
reST `::` literal block), and put `>>>` lines only in the `Examples` section.

```python
def species_masses(network: Network) -> np.ndarray:
    """Return the mass of every species in a network.

    The masses follow the order of
    :attr:`~jaff.core.network.Network.species`.

    Parameters
    ----------
    network : Network
        Chemical reaction network.

    Returns
    -------
    np.ndarray
        Species masses in g.

    Raises
    ------
    ValueError
        If the network has no species.

    Examples
    --------
    >>> net = Network("COthin")
    >>> masses = species_masses(net)
    >>> print(masses[0])
    1.673773e-24
    """
    if not network.species:
        raise ValueError("Network has no species")

    return np.array([sp.mass for sp in network.species])
```

### Class Docstrings

Describe the class in the **class** docstring and document the constructor
arguments in the `__init__` docstring under a `Parameters` section, as
:class:`~jaff.codegen.codegen.Codegen` and :class:`~jaff.core.network.Network` do.
The API page generator builds a class page's "Constructor" section from the
`__init__` docstring (falling back to the class docstring when `__init__` has none),
and lists the class docstring's `Attributes` section as a table.

```python
class Codegen:
    """Generate rates, fluxes, ODEs, and Jacobians from a Network in multiple languages.

    Examples
    --------
    >>> net = Network("COthin")
    >>> cg = Codegen(network=net, lang="c++")
    >>> rates = cg.get_rates_str(use_cse=True)
    """

    def __init__(self, network: Network, lang: str = "c++") -> None:
        """Bind a network to a target language.

        Parameters
        ----------
        network : Network
            Parsed chemical reaction network.
        lang : str, optional
            Target language alias (``"c++"``, ``"c"``, ``"fortran"``,
            ``"python"``, ...).  Default is ``"c++"``.

        Raises
        ------
        InvalidLanguageError
            If *lang* is not a supported language.
        """
```

### Module Docstrings

Open each module with a one-line summary followed by a short description.
Reference the key classes the module exposes with `:class:` roles.

```python
"""Chemical reaction network module.

Exposes the :class:`~jaff.core.network.Network` class for loading and managing
chemical reaction networks from various file formats.

Examples
--------
>>> from jaff import Network
>>> net = Network("COthin")
"""
```

## Error Handling

Fail loudly with specific exceptions and descriptive messages.

### Specific Exceptions

Catch the narrowest exception that applies — never a bare `except:`.

```python
# Good
try:
    net = Network(filename)
except FileNotFoundError:
    print(f"File not found: {filename}")
except ValueError as e:
    print(f"Invalid network format: {e}")

# Avoid
try:
    net = Network(filename)
except:  # Too broad
    pass
```

### Custom Exceptions

Define a domain exception hierarchy and raise it from a meaningful base.

```python
class NetworkError(Exception):
    """Base exception for network errors."""
    pass

class ParseError(NetworkError):
    """Error parsing network file."""
    pass

def load_network(filename: str) -> Network:
    try:
        return Network(filename)
    except ValueError as e:
        raise ParseError(f"Failed to parse {filename}: {e}")
```

### Error Messages

State what was expected and what was actually received.

```python
# Good - descriptive
raise ValueError(f"Temperature must be positive, got {temp}")

# Avoid - vague
raise ValueError("Bad temperature")
```

<!--## Functions

### Single Responsibility

```python
# Good - single purpose
def calculate_rate(temperature: float, alpha: float, beta: float) -> float:
    return alpha * (temperature / 300) ** beta

def write_to_file(content: str, filename: str) -> None:
    with open(filename, 'w') as f:
        f.write(content)

# Avoid - multiple responsibilities
def calculate_and_save_rate(temperature, alpha, beta, filename):
    rate = alpha * (temperature / 300) ** beta
    with open(filename, 'w') as f:
        f.write(str(rate))
```

### Function Length

Keep functions short and focused (generally < 50 lines):

```python
# Good
def process_network(filename: str) -> Network:
    """Process network file."""
    net = load_network(filename)
    validate_network(net)
    return net

# If too long, split into smaller functions
```

### Arguments

```python
# Good - clear parameters
def compute_ode(
    species: list[str],
    reactions: list[Reaction],
    use_cse: bool = True
) -> str:
    pass

# Avoid - too many positional arguments
def compute_ode(species, reactions, cse, optimize, cache, debug):
    pass
```

## Classes

### Keep Classes Focused

```python
# Good - single responsibility
class Network:
    """Represents a chemical reaction network."""

    def load(self, filename: str) -> None:
        pass

    def validate(self) -> bool:
        pass

# Avoid - doing too much
class NetworkManager:
    """Does everything."""
    pass
```

### Properties vs Methods

```python
class Network:
    @property
    def num_species(self) -> int:
        """Number of species (property - no computation)."""
        return len(self._species)

    def calculate_rates(self, temperature: float) -> np.ndarray:
        """Calculate rates (method - requires computation)."""
        return self._compute_rates(temperature)
```
-->

## Testing Style

See the [Testing Guide](testing.md) for the full workflow; the rules below are
style only.

### Test Function Names

Name tests after the behaviour they verify, not `test1`/`test2`.

```python
# Good - descriptive
def test_network_loads_krome_format():
    pass

def test_codegen_raises_error_for_invalid_language():
    pass

# Avoid
def test1():
    pass
```

### Test Organization

Group related tests in a `Test<Subject>` class.

```python
class TestNetwork:
    """Tests for Network class."""

    def test_load_from_file(self):
        """Test loading network from file."""
        pass

    def test_validate_species(self):
        """Test species validation."""
        pass
```

### Assertions

Assert concrete values, not bare truthiness.

```python
# Good - clear assertions
assert len(net.species) == 35
assert net.label == "react_COthin"
assert_almost_equal(rate, 1.2e-10, decimal=15)

# Avoid - unclear
assert x
```

## Performance

Prefer the idiomatic, fast construct when it costs no clarity.

### List Comprehensions

Build lists with comprehensions instead of `append` loops.

```python
# Good - fast
species_names = [s.name for s in network.species]

# Slower
species_names = []
for s in network.species:
    species_names.append(s.name)
```

### String Formatting

Use f-strings over concatenation.

```python
# Good - f-strings (fast, readable)
message = f"Network has {n} species"

# Avoid - slow
message = "Network has " + str(n) + " species"
```

## Best Practices

General Python habits the codebase follows.

### Use Context Managers

Manage resources with `with` so they always close, even on error.

```python
# Good
with open(filename, 'w') as f:
    f.write(content)

# Avoid
f = open(filename, 'w')
f.write(content)
f.close()
```

### Use Pathlib

Use `pathlib.Path` for filesystem paths instead of `os.path`.

```python
from pathlib import Path

# Good
path = Path("networks") / "react.dat"
if path.exists():
    content = path.read_text()

# Avoid
import os
path = os.path.join("networks", "react.dat")
if os.path.exists(path):
    with open(path) as f:
        content = f.read()
```

### Early Returns

Return early to handle edge cases first and avoid deep nesting.

```python
# Good - early return
def process(value: int | None) -> int:
    if value is None:
        return 0

    if value < 0:
        return -1

    return value * 2

# Avoid - nested if
def process(value):
    if value is not None:
        if value >= 0:
            return value * 2
        else:
            return -1
    else:
        return 0
```

## Code Review Checklist

Before submitting code:

- [x] Changed files formatted with ruff (`ruff format --check`)
- [x] Imports ordered: standard library, third-party, local
- [x] No linting errors (`ruff check`) on changed files
- [x] No new type errors (`ty check`) on changed files
- [x] `ABOUTME` header on new files
- [x] Type hints added
- [x] Docstrings written
- [x] Tests added/updated
- [x] No commented-out code
- [x] No hardcoded paths
- [x] Error handling in place

## See Also

- [Testing Guide](testing.md)
- [PEP 8](https://pep8.org/)
- [numpydoc style guide](https://numpydoc.readthedocs.io/en/latest/format.html)
