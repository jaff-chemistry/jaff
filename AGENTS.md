# AGENTS.md - Agent Guidelines for JAFF

JAFF (Just Another Fancy Format) is an astrochemical network parser that parses chemical reaction network files and generates code in multiple languages (C++, C, Fortran, Python).

## Build, Lint, and Test Commands

### Installation

```bash
uv pip install --system -e ".[dev]"    # Install with dev dependencies
```

### Running Tests

```bash
uv run pytest                           # Run all tests
uv run pytest -v                        # Verbose output
uv run pytest tests/test_network.py     # Run specific test file
uv run pytest tests/test_network.py::test_load_network  # Run single test
uv run pytest -k "network"              # Run tests matching pattern
uv run pytest --cov=jaff                # Run with coverage
```

### Code Formatting & Linting

```bash
ruff format src/ tests/                 # Format code
ruff format --check src/                # Check formatting without modifying
ruff check src/                         # Lint code
ruff check --fix src/                   # Fix linting issues
```

## Code Style Guidelines

### Python Version

- Support Python 3.11+
- Use `List`, `Dict` from `typing` (not built-in generics)

### Formatting (Ruff)

- Line length: 90 characters, Indent width: 4 spaces
- Quote style: double quotes, Format docstrings with code blocks

### Naming Conventions

- **Variables/functions**: `snake_case` (e.g., `user_name`, `calculate_rate`)
- **Classes**: `PascalCase` (e.g., `Network`, `CodeGenerator`)
- **Constants**: `UPPER_SNAKE_CASE` (e.g., `MAX_ITERATIONS`)
- **Private members**: prefix with underscore (e.g., `self._internal_data`)

### Type Hints

- Always use type hints for function signatures
- Import from `typing`: `List`, `Dict`, `Optional`, `Union`, `Tuple`, `Any`
- Use `pathlib.Path` for file paths

```python
from typing import List, Dict, Optional, Union, Tuple, Any
from pathlib import Path
import numpy as np

def process_species(
    names: List[str],
    masses: np.ndarray,
    options: Optional[Dict[str, Any]] = None
) -> Tuple[List[str], np.ndarray]:
    pass
```

### Docstrings (Google Style)

```python
def compute_rates(network: Network, temperature: float) -> np.ndarray:
    """Compute reaction rate coefficients.

    This function calculates rate coefficients for all reactions
    in the network at the specified temperature.

    Args:
        network: Chemical reaction network
        temperature: Gas temperature in Kelvin

    Returns:
        Array of rate coefficients in cm³/s

    Raises:
        ValueError: If temperature is negative
    """
```

### Error Handling

- Use specific exceptions (not bare `except:`)
- Provide descriptive error messages with context

```python
# Good
try:
    net = Network(filename)
except FileNotFoundError:
    raise FileNotFoundError(filename)
except ValueError as e:
    raise ValueError(f"Invalid network format: {e}")
```

### Best Practices

- Use f-strings for string formatting
- Use context managers (`with open(...)`)
- Use `pathlib.Path` for file paths
- Prefer early returns to reduce nesting
- Keep functions focused and under 50 lines
- Use list comprehensions over loops when appropriate
- Test one thing per test function

### Imports Order

1. Standard library, 2. Third-party packages, 3. Local imports

```python
import os
from pathlib import Path
from typing import List, Dict, Optional
import numpy as np
from jaff import Network
```

### Testing Conventions

- Test file naming: `test_<module>.py`
- Test class naming: `Test<ClassName>`
- Test function naming: `test_<description>()`
- Use pytest fixtures in `conftest.py`
- Use `@pytest.mark.parametrize` for parameterized tests

### Special Project Rules

1. **File headers**: Each code file starts with a 2-line comment:
    ```python
    # ABOUTME: <short description>
    # ABOUTME: <optional second line>
    ```
2. **Coverage goals**: >80% overall, >90% for critical modules
3. **Preserve existing comments** unless demonstrably false
4. **TDD approach**: Write tests BEFORE implementing new features

## Directory Structure

```
src/jaff/       # Main source code
tests/          # Test files
  fixtures/     # Test data
examples/       # Example marimo notebooks (run by tests/test_example_notebooks.py)
docs/           # Documentation
```

## Key Modules

- `network.py`: Network class for loading/managing chemical networks
- `codegen.py`: Multi-language code generation
- `file_parser.py`: File format parsing
- `reaction.py`: Reaction representations
- `species.py`: Species representations
