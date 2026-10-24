# AGENTS.md - Agent Guidelines for JAFF

JAFF (Just Another Fancy Format) is an astrochemical network parser that parses chemical reaction network files and generates code in multiple languages (C, C++, Fortran, Python, Rust, Julia, R).

## Build, Lint, and Test Commands

### Installation

```bash
uv sync --extra dev                     # Install with dev dependencies
```

Always run tools as `uv run --extra dev ...`. Plain `uv run` does an inexact sync: it
leaves already-installed dev extras in place but will not install them, so on a fresh
environment `pytest`, `ruff`, `ty` or `zensical` may be missing.
Data tables download on first use; set `JAFF_OFFLINE=1` to forbid network access.

### Running Tests

```bash
uv run --extra dev pytest                # Run all tests
uv run --extra dev pytest -v             # Verbose output
uv run --extra dev pytest tests/core/network                  # One package's tests
uv run --extra dev pytest tests/core/network/test_network.py  # One test file
uv run --extra dev pytest tests/core/network/test_network.py::TestConstruction  # One class
uv run --extra dev pytest -k "network"   # Run tests matching pattern
uv run --extra dev pytest --cov=jaff     # Run with coverage
```

### Code Formatting & Linting

```bash
uv run --extra dev ruff format --check <files>  # Check formatting without modifying
uv run --extra dev ruff format <files>          # Format the files you changed
uv run --extra dev ruff check <files>           # Lint code
uv run --extra dev ty check <files>             # Type check (ruff misses type errors)
```

Only format and fix the files you touched, and only if they were already formatted:
`tests/` is not lint-clean and a repo-wide `ruff format` / `--fix` rewrites unrelated
files. Pre-existing `ty` errors exist; make sure you add none.

## Code Style Guidelines

### Python Version

- Support Python 3.11+
- Start modules with `from __future__ import annotations`; use built-in generics
  (`list[str]`, `dict[str, Any]`) and `X | None`, not `List`/`Dict`/`Optional`

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
- Built-in generics and PEP 604 unions; import only what has no built-in form
  (`Any`, `Callable`, `Mapping`, `Sequence`, ...) from `typing`
- Use `pathlib.Path` for file paths
- Older modules still use `List`/`Dict`/`Optional`; convert them when you touch them

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

### Docstrings (numpydoc)

Docstrings are the only source of the API reference (`docs/reference/`, generated
by `uv run --extra dev python -m jaff._utils.api_docs`) and of the config key
descriptions (schemas and `docs/user-guide/configuration/reference.md`). Use numpydoc
sections; document config keys in an `Attributes` section, worded for someone
editing a TOML file.

```python
def compute_rates(network: Network, temperature: float) -> np.ndarray:
    """Compute reaction rate coefficients.

    This function calculates rate coefficients for all reactions
    in the network at the specified temperature.

    Parameters
    ----------
    network : Network
        Chemical reaction network.
    temperature : float
        Gas temperature in Kelvin.

    Returns
    -------
    np.ndarray
        Rate coefficients in cm^3 s^-1.

    Raises
    ------
    ValueError
        If temperature is negative.
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
from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import numpy as np

from jaff import Network
```

### Testing Conventions

- Test file naming: `test_<module>.py`
- Test class naming: `Test<ClassName>`
- Test function naming: `test_<description>()`
- Use pytest fixtures in `conftest.py`
- Use `@pytest.mark.parametrize` for parameterized tests
- Compare floats with `pytest.approx` (also numbers inside snapshots)

## Tests

`tests/` mirrors `src/jaff/`: put a test in the folder of the package whose code it
tests (`tests/core/network/`, `tests/physics/dust/`, `tests/codegen/templating/`, ...).
Cross-cutting tests go in `tests/e2e/` (golden code, generated solvers, bundled
networks, snapshots) or `tests/meta/` (public API, API docs, imports, notebooks). A
module's config tests sit with the module; `tests/config/` holds only the loader,
schema and reference machinery. Every folder has an `__init__.py`.

- **Data paths**: import from `tests._paths` (`FIXTURES_DIR`, `SNAPSHOTS_DIR`,
  `GOLDEN_DIR`, `REPO_ROOT`); never build paths from `Path(__file__)`.
- **Shared helpers** live at `tests/` root and are imported as
  `from tests.xsec_props import ...` / `from tests.codegen_render import ...`.
- **Fixtures** (`tests/conftest.py`): `make_network(lines, name=..., **kwargs)` builds a
  `Network` from inline reaction lines in `tmp_path`; `fixtures_dir`,
  `sample_network`; an autouse fixture silences the banner and logging.
- **Regenerating references**: `JAFF_UPDATE_GOLDEN=1` (golden code in `tests/golden/`),
  `JAFF_UPDATE_SNAPSHOTS=1` (snapshots); review the diff before keeping it.
- **Reproducing a rendered-code failure**: `JAFF_GOLDEN_SEED=<seed>` (the seed is in
  the failure message).

## Configuration

`jaff.toml` and `jaffgen.toml` share one schema built from pydantic models. Each
configurable module owns a `config.py` next to its code (e.g.
`physics/dust/config.py`), and a parent model lists its children's models as fields,
so the TOML keys mirror the class structure (`[network.dust.pe]`). One merge rule
(`config/loader.py`: tables key by key, `eos` / shielding tables whole, one
validation at the end) builds every config. `Network` owns `[network]`: defaults <
the `jaff.toml` next to the network (always applied: the network's own defaults) <
`config=` (one source: a TOML file, a `dict` of the `[network]` keys or a
`NetworkConfig`) < keyword args. jaffgen layers `jaffgen.toml` < CLI options <
`--set key.path=value` over that.

`Network(...)`'s keyword arguments are the `[network]` keys, named and nested as in
the file: `strict`, `label`, `funcfile`, `expand_nuclei`, and the tables `reactions`,
`radiation`, `dust`, `thermodynamics` (each a `dict` or its config model), e.g.
`Network("GOW", reactions={"photo": {"use_proxy_photoreaction": True}},
thermodynamics={"eos": {"gamma": 1.4}})`. Settings are read back from `net.config`
(same structure); don't copy config values onto other attributes.

To add a setting, follow `docs/development/adding-configuration.md`. In short:

1. Add the field to the module's `ConfigModel` subclass (`jaff.config.base`; strict,
   frozen, unknown keys rejected) and document it in the class's numpydoc
   `Attributes` (that text becomes the schema and reference description).
2. A new `config.py` gets an `__all__`, a field on its parent model, an entry in
   `CONFIG_MODULES` (`tests/config/test_config_modules.py`) and one in
   `CONFIG_MODULES` of `src/jaff/_utils/api_docs.py`.
3. The class doing the work takes its own section (`Thing(network, config)`), never
   the whole config, keeps it read-only as `self.config`, passes `config.<child>` to
   the components it builds, and holds only computed state otherwise: no copies of
   settings, no reads of another component's settings from `net.config`. Reactions
   do the same: the catalogue is `Reactions(config.reactions)` and owns the reaction
   settings (`net.reactions.config.resolve(s)` gives a reaction's `reaction.config`,
   `.own(s)` the keys set on it alone; `add_rate` applies `duplicate_policy`;
   `check_config()` checks per-reaction entries), and a
   photo reaction's `reaction.photo` is a `PhotoReaction(reaction, config.photo)`
   component holding its photo settings, computed photo state and the photo logic
   (`compute()`, `band_xsecs`, `plot_xsecs()`); every reaction's `reaction.rad` is a
   `ReactionRadiation(reaction)` component (its `.jfunc` spectral yield, read not
   copied, and the computed `group_yields`). A reaction component lists every
   expression it holds in `expressions()` and rewrites them in
   `substitute(mapping)`: `Reactions.link()` resolves `rc_<N>` through them and
   `tests/e2e/test_no_shorthand_survives.py` checks them. `.jfunc` functions are
   read through `net.functions` (`value(name)` / `call(name, args)`) or
   `reaction.function(kind)`, already standardized. Checks against the network
   raise `ParserError` with the key path.
4. `Network`'s keyword arguments are exactly the `[network]` keys: a new top-level
   field of `NetworkConfig` needs a keyword of the same name (a test enforces it);
   nothing else maps keywords. Optional CLI shortcut: an entry in `CLI_KEYS`
   (`src/jaff/cli/_helper.py`).
5. Regenerate the committed editor schemas:
   `uv run --extra dev python -m jaff.config.json_schema`
   (`tests/config/test_config_json_schema.py` fails when they are stale).

There are no aliases for old keys: a renamed key is a hard break.

## Documentation

- Docs are built with Zensical: `uv run --extra dev zensical build --strict`.
- The API reference (`docs/reference/`) and the configuration reference
  (`docs/user-guide/configuration/reference.md`) are generated from docstrings and
  are **not committed** (gitignored; CI generates them). Generate them before a local
  build or serve:

  ```bash
  uv run --extra dev python -m jaff._utils.api_docs   # API pages + API nav in zensical.toml
  uv run --extra dev python -m jaff.config.reference  # configuration reference
  ```

- The API nav section of `zensical.toml` *is* committed; regenerate it when exports
  change (`tests/meta/test_api_docs.py` checks it).
- Documented packages are `PACKAGES` and `CONFIG_MODULES` in `src/jaff/_utils/api_docs.py`;
  their `__all__` (strings only) decides what gets a page.
- Write docstring markup as reST-style ``code``, :class:`~a.b.C` (tilde = last name
  only), fenced code blocks (never `::` / `.. math::`), and `>>>` only in `Examples`.
- How docstrings become pages: `docs/development/api-documentation.md`. In short:
  - a class page = description, Constructor (from the `__init__` docstring's
    `Parameters`), the class `Attributes` table, then `Notes`, `Examples`, and a
    page per public method/property the class itself defines (inherited ones get none);
  - each member's first docstring line is its row text in the tables, so make it
    stand alone;
  - one `Parameters` entry per name (grouped `a, b : float` gives both the first type);
  - a constant's description is the string literal right after its assignment;
  - config models' `Attributes` also become the JSON-schema descriptions and the
    configuration reference: write them for TOML users.
- `Examples` are not run in CI (Network prints a banner): run them and paste real output.

## Special Project Rules

1. **File headers**: Each code file starts with a 2-line comment:
    ```python
    # ABOUTME: <short description>
    # ABOUTME: <optional second line>
    ```
    Not every existing module has one yet and nothing checks it; add the header to
    new files and to files you edit.
2. **Coverage goals**: >80% overall, >90% for critical modules
3. **Preserve existing comments** unless demonstrably false
4. **TDD approach**: Write tests BEFORE implementing new features

## Directory Structure

```
src/jaff/       # Main source code
tests/          # Test files, mirroring src/jaff/ (core/, physics/, codegen/, ...)
                # plus e2e/ and meta/; paths via tests/_paths.py
  fixtures/     # Test data
examples/       # Example marimo notebooks (run by tests/meta/test_example_notebooks.py)
docs/           # Documentation
```

## Key Modules

- `core/network/network.py`: `Network`, loading/managing a chemical network
- `core/network/ode_system.py`: `OdeSystem` (`net.ode_system(...)`), the state vector and
  rhs a solver integrates; `Codegen` prints it rather than assembling equations
- `core/reaction/`, `core/species/`, `core/elements/`: reactions, species, elements
- `core/parsers/`: network file formats and `.jfunc` auxiliary functions
- `physics/`: radiation and photo-reactions, dust, thermodynamics/EOS, constants
- `codegen/`: `Codegen`, `Language`s, the `$$JAFF` template engine (`templating/`)
- `config/`: config base models, loader, JSON schemas, reference generator
- `cli/`: `jaffgen` and `jaffx`
- `drivers/`, `io/`: HDF5/SQLite/TOML/CSV access, data downloads, logging
