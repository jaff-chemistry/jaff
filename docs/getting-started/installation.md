---
tags:
    - Introduction
    - Installation
icon: phosphor/download-simple
---

# Installation

## Requirements

JAFF requires `python>=3.11`. Check your Python version:

```bash
python --version
```

## Installation

<!-- prettier-ignore -->
!!! tip "Managing virtual environments"
    Although JAFF can be installed in a variety of ways, it is encouraged to use [uv](https://docs.astral.sh/uv/) to install JAFF and manage its virtual environment.

We recommend installing JAFF inside a virtual environment to avoid conflicts with other packages.

### 1. Create and activate a virtual environment

=== "venv"

    ```bash
    # Create virtual environment
    python -m venv .venv

    # Activate (Linux/macOS)
    source .venv/bin/activate

    # Activate (Windows)
    .venv\Scripts\activate
    ```

=== "conda"

    ```bash
    # Create conda environment
    conda create -n jaff python=3.11

    # Activate
    conda activate jaff
    ```

=== "uv"

    ```bash
    # Create virtual environment (faster)
    uv venv

    # Activate (Linux/macOS)
    source .venv/bin/activate

    # Activate (Windows)
    .venv\Scripts\activate
    ```

### 2. Install from source

=== "pip"

    ```bash
    # Clone the repository
    git clone https://github.com/jaff-chemistry/jaff.git
    cd jaff

    # Install the package
    pip install -e .
    ```

=== "uv"

    ```bash
    # Clone the repository
    git clone https://github.com/jaff-chemistry/jaff.git
    cd jaff

    # Install the package
    uv pip install -e .
    ```

## Verifying Installation

After installation, verify that JAFF is working correctly:

```bash
# Check JAFF version
python -c "import jaff; print(jaff.__version__)"

# Test the code generator CLI
jaffgen --help
```

You should see the JAFF code generator help message.

### Available Commands

After installation, JAFF provides the following command-line tools:

- **`jaffx`** - Quick network inspection (`jaffx get`) and export (`jaffx export`),
  for example a rate-coefficient table in HDF5:

```bash
jaffx export hdf5 --network networks/GOW/GOW.jet --file gow.hdf5
```

- **`jaffgen`** - Code generator for chemical reaction networks

```bash
jaffgen --network networks/GOW/GOW.jet --template microphysics --outdir generated/
```

## Data files

Some data tables (photo cross sections, shielding factors, radiation fields and
dust data) are not shipped with the package. They are downloaded on first use
with [pooch](https://www.fatiando.org/pooch/) into the `data/` directory of the
installed `jaff` package and reused from there afterwards, so the first load of
a network that needs them requires network access. Set `JAFF_OFFLINE=1` (any
non-empty value) to disable downloads; a table that is not already present then
cannot be loaded.

## Troubleshooting

### `ImportError`: No module named 'jaff'

Make sure you've activated your virtual environment and that the installation completed successfully.

### NumPy/SciPy Installation Issues

On some systems, if you are using pip, you may need to install `NumPy`, `SciPy`,
and `Astropy` separately:

```bash
# Install scientific stack first
pip install numpy scipy astropy
```

## Updating JAFF

JAFF is installed from source, so update by pulling the latest changes and reinstalling:

=== "pip"

    ```bash
    cd jaff
    git pull
    pip install -e .
    ```

=== "uv"

    ```bash
    cd jaff
    git pull
    uv pip install -e .
    ```

## Uninstalling

To remove JAFF from your system/venv:

=== "pip"

    ```bash
    pip uninstall jaff
    ```

=== "uv"

    ```bash
    uv pip uninstall jaff
    ```

## Next Steps

Now that JAFF is installed:

- **Quick Start**: Follow the [Quick Start Guide](quickstart.md) to run your first network analysis
- **Basic Concepts**: Learn about [chemical networks and reactions](concepts.md)
- **User Guide**: Explore detailed [usage documentation](../user-guide/working-with-networks/index.md)
- **Developer Guide**: If you want to contribute or develop JAFF, see the [Developer Installation Guide](../development/installation.md)
