# Quick Start

## Basic Usage

### Loading a Network

```python
from jaff import Network

# Load a chemical network from file
network = Network("networks/react_COthin")
```

### Accessing Species

```python
# Get total number of species
print(f"Total species: {len(network.species)}")

# Iterate through species
for species in network.species:
    print(f"{species.name}: mass={species.mass:.2f}, charge={species.charge:+d}")

# Get a specific species by name
co_species = network.get_species_object("CO") 
```

### Accessing Reactions

```python
# Get total number of reactions
print(f"Total reactions: {len(network.reactions)}")

# Iterate through reactions
for i, reaction in enumerate(network.reactions):
    print(f"Reaction {i+1}: {reaction}")

# Get reaction by its string representation
reaction_idx = network.get_reaction_index("C + H3+ -> CH+ + H2")
```

### Generating Rate Tables

```python
# Create a temperature-dependent rate table
temp, rates = network.get_table(
    T_min=10,      # Minimum temperature (K)
    T_max=1000,    # Maximum temperature (K) 
    nT=64,         # Number of temperature points
    err_tol=0.01   # Error tolerance for adaptive sampling
)

print(f"Rate table shape: {rates.shape}")  # (n_reactions, n_temperatures)
```

## Command Line Usage

### Basic Network Loading

```bash
# Load and validate a network
jaff networks/react_COthin
```

### List Species and Reactions

```bash
# List all species
jaff networks/react_COthin --list-species

# List all reactions  
jaff networks/react_COthin --list-reactions
```

### Validation Checks

```bash
# Exit on validation errors
jaff networks/react_COthin --errors

# Set a custom network label
jaff networks/react_COthin --label "CO_thin_network"
```

## Network Properties

The Network class automatically performs several validation checks:

- **Sink/Source Detection**: Identifies species that only appear as products (sources) or reactants (sinks)
- **Recombination Checks**: Verifies that positive ions have electron recombination reactions
- **Isomer Detection**: Finds species with identical elemental composition but different names
- **Duplicate Reactions**: Identifies reactions that appear multiple times

All validation results are printed during network loading.