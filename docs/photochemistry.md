# Photochemistry

JAFF provides support for photochemical reactions including photoionization and photodissociation processes, with automatic loading of cross-section data.

## Overview

The `Photochemistry` module handles:
- Loading cross-section data from files
- Automatic detection of photoionization vs photodissociation
- Energy-dependent cross-section retrieval
- Integration with the reaction network

## Cross-Section Data

### File Format

Cross-section data files should be placed in `src/jaff/data/xsecs/` directory. The files follow the Leiden database format:

```
# Example: CO__C_O.dat
# Photodissociation of CO -> C + O
# wavelength(nm)  dissociation(cm2)
111.8    2.23e-17
112.0    2.31e-17
112.5    2.45e-17
...
```

### File Naming Convention

Files must be named as: `REACTANTS__PRODUCTS.dat`

- Reactants and products separated by double underscore `__`
- Multiple species separated by single underscore `_`
- Species names must match those in your network
- Example: `H2O__H_OH.dat` for H2O → H + OH

### Data Columns

The data file should contain:
- **Wavelength**: in nanometers (nm)
- **Cross-section**: in cm²

The module automatically identifies which column contains the cross-section based on the reaction type (ionization vs dissociation).

## Using Photochemistry

### Basic Usage

```python
from jaff import Network
from jaff.photochemistry import Photochemistry

# Load network with photoreactions
network = Network("networks/test.dat")

# Initialize photochemistry module
photo = Photochemistry()

# Get cross-section data for a specific reaction
reaction = network.get_reaction_by_verbatim("H -> H+ + e-")
xsec_data = photo.get_xsec(reaction)

# Access the data
energies = xsec_data['energy']  # in erg
xsecs = xsec_data['xsecs']      # in cm²
```

### Automatic Loading

The `Photochemistry` class automatically loads all cross-section files from the data directory on initialization:

```python
photo = Photochemistry()
# All .dat files in data/xsecs/ are now loaded
```

## Energy Calculations

The module converts wavelengths to photon energies:

```python
# Energy (erg) = h * c / wavelength
# where wavelength is converted from nm to cm
```

Constants used:
- Speed of light: c = 2.99792458×10¹⁰ cm/s
- Planck constant: h = 6.62607015×10⁻²⁷ erg·s

## Reaction Types

The module automatically determines reaction type based on charge balance:

### Photoionization
Reactions where products have higher charge than reactants:
```
H + photon -> H+ + e-
CO + photon -> CO+ + e-
```

### Photodissociation
Reactions where charge is conserved:
```
H2 + photon -> H + H
H2O + photon -> H + OH
```

## Integration with Networks

When working with networks containing photoreactions:

```python
# Find all photoreactions
photo_reactions = [r for r in network.reactions 
                   if 'photon' in [s.name for s in r.reactants]]

# Get cross-sections for all photoreactions
photo = Photochemistry()
for reaction in photo_reactions:
    try:
        data = photo.get_xsec(reaction)
        print(f"Loaded data for: {reaction.get_verbatim()}")
    except SystemExit:
        print(f"Missing data for: {reaction.get_verbatim()}")
```

## Adding New Cross-Section Data

To add cross-section data for a new photoreaction:

1. **Obtain data**: Get wavelength-dependent cross-sections from literature or databases
2. **Format file**: Create a `.dat` file with the proper format
3. **Name correctly**: Use the `REACTANTS__PRODUCTS.dat` naming scheme
4. **Add header**: Include comment lines describing the reaction and units
5. **Place in directory**: Copy to `src/jaff/data/xsecs/`

Example file creation:
```python
import numpy as np

# Your data
wavelengths = np.array([100, 110, 120, 130])  # nm
cross_sections = np.array([1e-17, 2e-17, 3e-17, 2.5e-17])  # cm²

# Write file (example - modify path as needed)
filename = 'NH3__NH2_H.dat'
with open(filename, 'w') as f:
    f.write("# Photodissociation of NH3\n")
    f.write("# NH3 + photon -> NH2 + H\n")
    f.write("# wavelength(nm)  dissociation(cm2)\n")
    for wl, xs in zip(wavelengths, cross_sections):
        f.write(f"{wl:.1f}    {xs:.2e}\n")
```

## Common Issues

### Missing Data Files

If a photoreaction in your network doesn't have corresponding cross-section data:

```
ERROR: reaction H2O_photon__H_OH+ not found in photochemistry data.
Add the file to the data/xsecs folder as H2O__H_OH+.dat
```

### Reaction Serialization

The module uses "serialized" reaction names with sorted species:
- Input: `H + OH + photon -> H2O`
- Serialized: `H_OH__H2O`

This ensures consistent matching regardless of species order in the network file.

## Best Practices

1. **Verify units**: Ensure wavelengths are in nm and cross-sections in cm²
2. **Check completeness**: Include data across relevant wavelength ranges
3. **Document sources**: Add references in file headers
4. **Test loading**: Verify new data loads correctly before using in calculations
5. **Handle missing data**: Check for reactions without cross-section data