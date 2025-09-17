# Python API

## Core Classes

### Network Class

The main class for loading and working with chemical networks.

```python
from jaff import Network

# Initialize a network (optionally pass an auxiliary functions file)
network = Network(
    "networks/react_COthin",
    errors=False,
    label=None,
    funcfile=None,   # or path to a functions file, or 'none'
)
```

**Parameters:**
- `fname` (str): Path to the network file
- `errors` (bool): If True, exit on validation errors (default: False)
- `label` (str): Custom label for the network (default: filename)
- `funcfile` (str|None): Optional auxiliary functions file. If `None`, JAFF will look for a file named `fname + "_functions"`. If `'none'`, auxiliary function parsing is disabled. Any functions found are substituted into rate expressions before processing.

### Species Access

```python
# Access species list
species_list = network.species

# Get number of species
n_species = len(network.species)

# Get species by name
species = network.get_species_object("CO")

# Get species index
idx = network.get_species_index("CO")
```

### Reaction Access

```python
# Access reactions list
reactions_list = network.reactions

# Get number of reactions
n_reactions = len(network.reactions)

# Get reaction by verbatim string
idx = network.get_reaction_index("H + e- -> H+ + e- + e-")
```

### Rate Table Generation

Generate temperature-dependent rate coefficient tables:

```python
temps, rates = network.get_table(
    T_min=10,           # Minimum temperature (K)
    T_max=1000,         # Maximum temperature (K)  
    nT=64,              # Initial number of temperature points
    err_tol=0.01,       # Relative error tolerance
    rate_min=1e-30,     # Minimum rate for error calculation
    rate_max=1e100,     # Maximum rate (clipped)
    fast_log=False,     # Use fast_log2 sampling (default: False)
    verbose=False       # Print adaptive refinement info
)
```

**Returns:** 
- `temps`: `numpy.ndarray` of temperature values (K)
- `rates`: `numpy.ndarray` with shape `(n_reactions, n_temperatures)`

### Table Export

Export rate tables to disk in text or HDF5 format:

```python
network.write_table(
    fname="rates.hdf5",    # Output filename (.txt or .hdf5/.hdf)
    T_min=10,              # Minimum temperature (K)
    T_max=1000,            # Maximum temperature (K)
    nT=64,                 # Initial number of temperature points
    err_tol=0.01,          # Relative error tolerance
    rate_min=1e-30,        # Minimum rate for error calculation
    rate_max=1e100,        # Maximum rate (clipped)
    fast_log=False,        # Use fast_log2 sampling
    format='auto',         # 'auto', 'txt', or 'hdf5'
    include_all=False,     # Include non-tabulated reactions as NaN
    verbose=False          # Print adaptive refinement info
)
```

The HDF5 format follows the Quokka table standard with proper metadata attributes.

### Network Comparison

Compare two networks:

```python
network1 = Network("networks/react_COthin")
network2 = Network("networks/react_COthin")

# Compare reactions
network1.compare_reactions(network2, verbosity=1)

# Compare species
network1.compare_species(network2, verbosity=1)
```

## Species Properties

Individual species objects have these attributes:

```python
species = network.get_species_object("CO")

print(species.name)        # Species name
print(species.mass)        # Molecular mass
print(species.charge)      # Electric charge
print(species.index)       # Index in species list
print(species.exploded)    # Elemental composition
print(species.latex)       # LaTeX representation
```

## Reaction Properties

Individual reaction objects have these attributes:

```python
reaction = network.reactions[0]

print(reaction.reactants)  # List of reactant Species objects
print(reaction.products)   # List of product Species objects
print(reaction.rate)       # Sympy rate expression
print(reaction.tmin)       # Minimum temperature (if set)
print(reaction.tmax)       # Maximum temperature (if set)
```

## Utility Methods

```python
# Get LaTeX representation of species
latex_str = network.get_latex("CO", dollars=True)

# Get reaction verbatim string
verbatim = network.get_reaction_verbatim(idx)

# Check if networks have common elements
# (methods exist but require careful examination of their implementation)
```

## Primitive Variables

Rate expressions can contain these primitive variables:

- `tgas`: Gas temperature (K)
- `av`: Visual extinction (Draine units)
- `crate`: Cosmic ray ionization rate of H₂ (s⁻¹)
- `ntot`: Total number density (cm⁻³)
- `hnuclei`: H nuclei number density (cm⁻³)
- `d2g`: Dust-to-gas mass ratio

## Code Generation

Generate ODE solver code from a network:

```python
from jaff.builder import Builder

# Create builder instance
builder = Builder(network)

# Generate Python solver code
builder.build(template="python_solve_ivp")

# Generate Fortran solver code
builder.build(template="fortran_dlsodes")

# Generate C++ solver (header-only integrators)
builder.build(template="kokkos_ode")
```

Generated code is placed in the `builds/` directory.

### Analytical Jacobian and CSE

- `Network.get_symbolic_ode_and_jacobian(idx_offset=0, use_cse=True, language="c++")` returns strings for ODE RHS and analytic Jacobian, with optional common subexpression elimination (CSE) and language-specific indexing.
- `Network.get_rates(idx_offset=0, rate_variable="k", language="python|f90|c++", use_cse=True)` generates rate coefficient code for the selected language. CSE significantly reduces redundant computations in C++ codegen.

## Photochemistry

Access photochemical cross-section data:

```python
from jaff.photochemistry import Photochemistry

# Initialize photochemistry module
photo = Photochemistry()

# Get cross-section data for a reaction
# First load a network with photoreactions
from jaff import Network
network = Network("networks/react_COthin")
# Example: access photoreaction data if it exists
# reaction = network.reactions[0]  # Get first reaction
# data = photo.get_xsec(reaction)
# Returns dict with keys: 'energy' (erg), 'xsecs' (cm²)
```

Cross-section data files should be placed in `data/xsecs/` in Leiden database format.

## Fast Logarithm Functions

For performance-optimized temperature sampling:

```python
from jaff.fastlog import fast_log2, inverse_fast_log2

# Fast approximation of log2(x)
y = fast_log2(100.0)

# Inverse function (machine precision)
x = inverse_fast_log2(y)
```
