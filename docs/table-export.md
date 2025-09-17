# Table Export

JAFF provides functionality to export reaction rate coefficient tables in multiple formats, suitable for use in hydrodynamics simulations and other applications.

## Overview

The `write_table()` method allows you to export temperature-dependent rate coefficients to disk in either text or HDF5 format. The tables use adaptive temperature sampling to achieve a specified interpolation accuracy while minimizing table size.

## Basic Usage

```python
from jaff import Network

# Load network
network = Network("networks/react_COthin")

# Export to HDF5 format
network.write_table("rates.hdf5", T_min=10, T_max=1000)

# Export to text format
network.write_table("rates.txt", T_min=10, T_max=1000)
```

## Parameters

### Required Parameters

- `fname` (str): Output filename. The extension determines the format if `format='auto'`

### Temperature Range

- `T_min` (float): Minimum temperature in Kelvin. If None, uses the minimum from all reactions
- `T_max` (float): Maximum temperature in Kelvin. If None, uses the maximum from all reactions

### Sampling Control

- `nT` (int): Initial number of temperature points (default: 64)
- `err_tol` (float): Relative error tolerance for adaptive sampling (default: 0.01)
- `fast_log` (bool): Use fast_log2 sampling instead of natural log (default: False)

### Rate Limits

- `rate_min` (float): Minimum rate for error calculation (default: 1e-30)
- `rate_max` (float): Maximum rate - values above this are clipped (default: 1e100)

### Output Control

- `format` (str): Output format - 'auto', 'txt', or 'hdf5' (default: 'auto')
- `include_all` (bool): Include non-tabulated reactions as NaN (default: False)
- `verbose` (bool): Print progress during adaptive refinement (default: False)

## Adaptive Sampling

The table export uses adaptive temperature sampling to achieve the specified error tolerance:

1. Start with `nT` logarithmically-spaced temperature points
2. For each pair of adjacent points, evaluate the exact rate at the midpoint
3. Compare with logarithmic interpolation between the endpoints
4. If error exceeds `err_tol`, add the midpoint to the table
5. Repeat until all intervals meet the tolerance

The error metric is: `abs((interp - exact) / (exact + rate_min))`

## Output Formats

### Text Format

The text format is compatible with the Quokka hydrodynamics code:

```
# JAFF auto-generated rate coefficient table
# Network name: my_network
# Reactions included
#   (reactants) (products) (reaction type)
#   {'H': 1, 'O': 1} {'OH': 1} 2_body
1                    # 1D table
25                   # Number of reactions
2                    # Log spacing (2) or fast_log spacing (3)
128                  # Number of temperature points
1.000000e+01 1.000000e+03  # T_min T_max
1.234567e-12 2.345678e-11 ... # Rate coefficients
```

### HDF5 Format

The HDF5 format follows the Quokka table standard with a group structure:

```
/reaction_coeff/              # Main group
    data                      # Dataset: (n_reactions, n_temps)
    output_names              # Dataset: reaction descriptions
    output_units              # Dataset: units for each reaction
    @input_names  = ['temperature']
    @input_units  = ['K']
    @xlo          = [T_min]
    @xhi          = [T_max]
    @spacing      = ['log'] or ['fast_log']
```

Note: Output names and units are stored as datasets rather than attributes to avoid HDF5 attribute size limitations when dealing with large reaction networks.

## Fast Logarithm Sampling

When `fast_log=True`, temperature points are sampled uniformly in fast_log2 space rather than natural log space. This can provide better sampling for some applications:

```python
# Standard logarithmic sampling
network.write_table("rates_log.hdf5", fast_log=False)

# Fast logarithm sampling
network.write_table("rates_fastlog.hdf5", fast_log=True)
```

## Filtering Reactions

By default, only reactions that can be tabulated as functions of temperature are included. Reactions depending on other variables (density, radiation field, etc.) are excluded unless `include_all=True`:

```python
# Include all reactions (non-tabulated ones as NaN)
network.write_table("rates_all.hdf5", include_all=True)

# Only tabulated reactions (default)
network.write_table("rates_tabulated.hdf5", include_all=False)
```

## Example: Complete Export

```python
from jaff import Network

# Load network
network = Network("networks/gas_reactions_kida.uva.2024.in")

# Export with custom settings
network.write_table(
    fname="kida_rates.hdf5",
    T_min=5.0,           # Very low temperature
    T_max=5000.0,        # High temperature  
    nT=128,              # Start with 128 points
    err_tol=0.01,        # 1% accuracy
    rate_min=1e-40,      # Very small rates
    rate_max=1e-5,       # Clip large rates
    fast_log=True,       # Use fast_log sampling
    format='hdf5',       # Explicit format
    include_all=False,   # Only tabulated reactions
    verbose=True         # Show progress
)
```

## Notes

- Reactions with temperature limits (`tmin`, `tmax`) are only evaluated within their valid range
- The adaptive sampling ensures accurate interpolation between table points
- For photochemical reactions, use `av=0` when generating tables
- Cosmic ray reactions use `crate=1` for table generation