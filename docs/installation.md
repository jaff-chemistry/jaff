# Installation

## Requirements

- Python 3.9 or higher
- Runtime dependencies: numpy, sympy, scipy, tqdm, h5py

## From Source

Clone the repository and install:

```bash
git clone https://github.com/tgrassi/jaff.git
cd jaff
pip install .
```

## Development Installation

For development work, install in editable mode with development dependencies:

```bash
git clone https://github.com/tgrassi/jaff.git
cd jaff
pip install -e ".[dev]"
```

This will install additional tools for development:
- pytest (testing framework)
- pytest-cov (coverage reporting)
- black (code formatting)
- ruff (linting)

## Using uv (Recommended)

If you have `uv` installed, you can use it for faster dependency resolution:

```bash
uv pip install .
```

Or for development:

```bash
uv pip install -e ".[dev]"
```

## Verification

After installation, you can verify JAFF is working by running:

```bash
jaff --help
```

You should see the command-line help output showing available options.

## Notes for C++ Code Generation

The `kokkos_ode` template uses a header-only integrators dependency. By default, its CMake project fetches this via `FetchContent`. In offline environments, configure CMake with `-DUSE_SYSTEM_INTEGRATORS=ON` and ensure the `integrators` package is installed and discoverable via `find_package`.
