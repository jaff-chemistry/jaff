---
tags:
    - Development
    - Installation
icon: phosphor/download-simple
---

# Installation

This guide covers a **development** install (editable, with dev tooling). For a
plain user install, see the [Getting Started installation guide](../getting-started/installation.md).

## Prerequisites

- `python >= 3.11`
- Git
- [uv](https://docs.astral.sh/uv/) (recommended) - a fast Python package and project manager; the repository is locked with it

Check your Python version:

```bash
python --version
```

## Development Installation

### Step 1: Clone the Repository

```bash
git clone https://github.com/jaff-chemistry/jaff.git
cd jaff
```

### Step 2: Install in Development Mode

The repository is locked with `uv` (`uv.lock`), and CI installs from that lock.
Install JAFF in editable mode with all development dependencies:

=== "uv (Recommended)"

    ```bash
    uv sync --extra dev
    ```

    This creates `.venv/` if needed and installs the exact locked versions, with JAFF
    itself in editable mode. You do not need to activate the environment: run tools
    with `uv run --extra dev ...`. Plain `uv run` will not install the dev extras
    (it leaves them in place if they are already installed), so keep the
    `--extra dev`.

=== "pip"

    ```bash
    python -m venv .venv
    source .venv/bin/activate     # Linux/macOS (.venv\Scripts\activate on Windows)
    pip install -e ".[dev]"
    ```

    This resolves the dependencies afresh instead of using `uv.lock`, so versions may
    differ from CI. With an activated environment, drop the `uv run --extra dev`
    prefix from the commands in these guides.

Editable mode means changes you make to the source code are picked up without
reinstalling.

## Verifying Your Development Setup

After installation, verify everything is working:

=== "Package import"

    ```bash
    uv run --extra dev python -c "import jaff; print(jaff.__version__)"
    ```

=== "CLI"

    ```bash
    uv run --extra dev jaffgen --help
    uv run --extra dev jaffx --help
    ```

=== "Tests"

    ```bash
    uv run --extra dev pytest
    ```

=== "Code style"

    ```bash
    uv run --extra dev ruff check src/
    uv run --extra dev ruff format --check src/
    ```

    `tests/` is not lint-clean, so check `src/` (or just the files you changed)
    rather than `.`. Type checking with `ty` is part of linting; see the
    [Code Style Guide](code-style.md#code-formatting).

## Dependencies

### Core Dependencies

These are automatically installed with JAFF:

| Package      | Minimum Version | Purpose                                  |
| ------------ | --------------- | ---------------------------------------- |
| `numpy`      | ≥2.0.0          | Numerical computations                   |
| `scipy`      | ≥1.13.0         | Scientific computing                     |
| `sympy`      | ≥1.14.0         | Symbolic mathematics                     |
| `pandas`     | ≥2.3.3          | Tabular data handling                    |
| `astropy`    | —               | Physical constants and units             |
| `matplotlib` | ≥3.9.4          | Plotting                                 |
| `h5py`       | ≥3.9.0          | HDF5 file support                        |
| `pooch`      | ≥1.9.0          | Downloading data files                   |
| `lark`       | ≥1.3.1          | Template-language grammar                |
| `pydantic`   | ≥2.14.1         | Configuration models and validation      |
| `tomli-w`    | ≥1.2.0          | Writing TOML                             |
| `typer`      | ≥0.27.0         | Command-line interfaces (`jaffgen`, `jaffx`) |
| `rich`       | ≥15.0.0         | Rich terminal output                     |
| `pygments`   | ≥2.19.2         | Syntax highlighting                      |
| `ipykernel`  | ≥7.2.0          | Jupyter kernel support                   |
| `marimo`     | ≥0.23.8         | Interactive notebooks                    |

### Development Dependencies

When installing with the `dev` extra, the following additional packages are installed:

#### Testing and Linting Tools

| Package            | Minimum Version | Purpose                          |
| ------------------ | --------------- | -------------------------------- |
| `pytest`           | ≥7.0            | Testing framework                |
| `pytest-cov`       | —               | Code coverage reporting          |
| `ruff`             | —               | Fast Python linter and formatter |
| `ty`               | —               | Type checker                     |
| `check-jsonschema` | —               | JSON schema validation           |
| `kaleido`          | —               | Static image export for plots    |

#### Documentation Tools

| Package                | Minimum Version | Purpose                             |
| ---------------------- | --------------- | ----------------------------------- |
| `zensical`             | —               | Documentation static-site generator |
| `griffelib`            | ≥2.3            | Reads docstrings for the API pages  |
| `pymdown-extensions`   | ≥10.5           | Advanced Markdown extensions        |

(`pygments` is also listed in the `dev` extra, without a version bound; the core
requirement above already installs it.)

## Updating Dependencies

After pulling new changes, bring your environment back in line with the lock file:

=== "uv"

    ```bash
    git pull
    uv sync --extra dev
    ```

=== "pip"

    ```bash
    git pull
    pip install -e ".[dev]" --upgrade
    ```

If you add or change a dependency in `pyproject.toml`, run `uv lock` and commit the
updated `uv.lock` with it: CI installs with `uv sync --locked`, which fails when the
lock file is out of date.

## Building Documentation

JAFF uses [Zensical](https://zensical.org/) as a static site generator for the documentation.

The API reference (`docs/reference/`) and the configuration reference
(`docs/user-guide/configuration/reference.md`) are generated from the docstrings
and are not committed; CI generates them before every build. Generate them
before serving or building locally:

```bash
uv run --extra dev python -m jaff._utils.api_docs   # API pages + the API nav in zensical.toml
uv run --extra dev python -m jaff.config.reference  # configuration reference page

# Serve with live reload (recommended during development)
uv run --extra dev zensical serve

# Build static site into site/
uv run --extra dev zensical build

# Build with strict link validation (what CI runs)
uv run --extra dev zensical build --strict
```

The API nav section of `zensical.toml` is committed. If `jaff._utils.api_docs`
changed it, commit the change: CI fails when the committed nav differs from the
generated one.

When serving, open `http://localhost:8000`. The site rebuilds automatically on save.

There is no manual deploy step: the docs workflow (`.github/workflows/docs.yml`)
builds the site on every push and pull request, and publishes it to GitHub Pages on
pushes to `main`.

## Working with uv (Recommended)

[uv](https://github.com/astral-sh/uv) is a fast Python package installer and resolver written in Rust, significantly faster than pip.

### Install uv

=== "macOS / Linux"

    ```bash
    curl -LsSf https://astral.sh/uv/install.sh | sh
    ```

=== "Windows"

    ```powershell
    powershell -c "irm https://astral.sh/uv/install.ps1 | iex"
    ```

=== "pip"

    ```bash
    pip install uv
    ```

### Common uv Commands

```bash
# Create/update .venv from uv.lock, with the dev tools
uv sync --extra dev

# Add a dependency (updates pyproject.toml and uv.lock)
uv add <package_name>

# Re-lock after editing pyproject.toml by hand
uv lock

# Run a command inside the project environment
uv run --extra dev python script.py
```

## IDE Setup

=== "VS Code"

    Recommended extensions:

    - Python
    - Pylance
    - Ruff
    - Even Better TOML

    Settings (`.vscode/settings.json`):

    ```json
    {
        "python.defaultInterpreterPath": ".venv/bin/python",
        "python.testing.pytestEnabled": true,
        "python.testing.unittestEnabled": false,
        "editor.formatOnSave": true,
        "python.linting.enabled": true,
        "ruff.format.args": ["--config", "pyproject.toml"]
    }
    ```

=== "Zed"

    Settings (`.zed/settings.json`):

    ```json
    {
        "languages": {
            "Python": {
                "language_servers": ["pyright"],
                "format_on_save": true
            }
        },

        "python": {
            "interpreter": ".venv/bin/python"
        },

        "formatter": {
            "external": {
                "command": "ruff",
                "arguments": ["format", "--config", "pyproject.toml", "-"]
            }
        },

        "lint": {
            "external": {
                "command": "ruff",
                "arguments": ["check", "--config", "pyproject.toml", "-"]
            }
        }
    }
    ```

## Next Steps

- **Running tests** — see the [Testing Guide](testing.md).
- **Code style and linting** — see the [Code Style Guide](code-style.md).
- **Submitting changes** — see the [Contributing Guide](contributing.md).
