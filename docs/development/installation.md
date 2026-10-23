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
- (Optional) [uv](https://docs.astral.sh/uv/) - A fast Python package installer and resolver

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

### Step 2: Create a Virtual Environment

It's strongly recommended to use a virtual environment for development:

=== "venv"

    ```bash
    # Create virtual environment
    python -m venv .venv

    # Activate (Linux/macOS)
    source .venv/bin/activate

    # Activate (Windows)
    .venv\Scripts\activate
    ```

=== "uv (Recommended)"

    ```bash
    # Create virtual environment with uv (much faster)
    uv venv

    # Activate
    source .venv/bin/activate  # Linux/macOS
    .venv\Scripts\activate     # Windows
    ```

=== "conda"

    ```bash
    # Create conda environment
    conda create -n .venv python=3.11

    # Activate
    conda activate .venv
    ```

### Step 3: Install in Development Mode

Install JAFF in editable mode with all development dependencies:

=== "pip"

    ```bash
    pip install -e ".[dev]"
    ```

=== "uv (Recommended)"

    ```bash
    uv pip install -e ".[dev]"
    ```

The `-e` flag installs the package in "editable" mode, meaning changes you make to the source code are immediately reflected without reinstalling.

## Verifying Your Development Setup

After installation, verify everything is working:

=== "Package import"

    ```bash
    python -c "import jaff; print(jaff.__version__)"
    ```

=== "CLI"

    ```bash
    jaffgen --help
    jaffx --help
    ```

=== "Tests"

    ```bash
    pytest
    ```

=== "Code style"

    ```bash
    ruff check .
    ```

## Dependencies

### Core Dependencies

These are automatically installed with JAFF:

| Package      | Minimum Version | Purpose                |
| ------------ | --------------- | ---------------------- |
| `numpy`      | ≥2.0.0          | Numerical computations |
| `scipy`      | ≥1.13.0         | Scientific computing   |
| `sympy`      | ≥1.14.0         | Symbolic mathematics   |
| `pandas`     | ≥2.3.3          | Tabular data handling  |
| `matplotlib` | ≥3.9.4          | Plotting               |
| `h5py`       | ≥3.9.0          | HDF5 file support      |
| `rich`       | ≥15.0.0         | Rich terminal output   |
| `pooch`      | ≥1.9.0          | Downloading data files |
| `pygments`   | ≥2.19.2         | Syntax highlighting    |
| `ipykernel`  | ≥7.2.0          | Jupyter kernel support |
| `marimo`     | ≥0.23.8         | Interactive notebooks  |

### Development Dependencies

When installing with `[dev]`, the following additional packages are installed:

#### Testing Tools

| Package            | Minimum Version | Purpose                          |
| ------------------ | --------------- | -------------------------------- |
| `pytest`           | ≥7.0            | Testing framework                |
| `pytest-cov`       | —               | Code coverage reporting          |
| `ruff`             | —               | Fast Python linter and formatter |
| `check-jsonschema` | —               | JSON schema validation           |

#### Documentation Tools

| Package                | Minimum Version | Purpose                             |
| ---------------------- | --------------- | ----------------------------------- |
| `zensical`             | —               | Documentation static-site generator |
| `griffelib`            | ≥2.3            | Reads docstrings for the API pages  |
| `pymdown-extensions`   | ≥10.5           | Advanced Markdown extensions        |
| `pygments`             | ≥2.19.2         | Syntax highlighting                 |

## Updating Dependencies

After pulling new changes, refresh your editable install:

=== "pip"

    ```bash
    git pull
    pip install -e ".[dev]" --upgrade
    ```

=== "uv"

    ```bash
    git pull
    # uv automatically resolves to the latest compatible versions
    uv pip install -e ".[dev]" --upgrade
    ```

## Building Documentation

JAFF uses [Zensical](https://zensical.org/) as a static site generator for the documentation.

The API reference (`docs/reference/`) and the configuration reference
(`docs/user-guide/configuration/reference.md`) are generated from the docstrings
and are not committed; CI generates them before every build. Generate them
before serving or building locally:

```bash
python -m jaff._utils.api_docs      # API pages + the API nav in zensical.toml
python -m jaff.config.reference     # configuration reference page

# Serve with live reload (recommended during development)
zensical serve

# Build static site into site/
zensical build

# Build with strict link validation
zensical build --strict

# Deploy to GitHub Pages (maintainers only)
zensical gh-deploy
```

When serving, open `http://127.0.0.1:8000`. The site rebuilds automatically on save.

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
# Create virtual environment
uv venv

# Install a package
uv pip install <package_name>

# Install JAFF in editable mode
uv pip install -e ".[dev]"

# Run a command inside the virtual environment
uv run python script.py
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
