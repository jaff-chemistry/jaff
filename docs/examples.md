# Examples

## Code Generation (C++)

```python
from jaff import Network
from jaff.builder import Builder

network = Network("networks/react_COthin")
Builder(network).build(template="kokkos_ode")
# Use CMake in the output directory to build and run the generated code
```

## Loading Different Network Formats

### KIDA Format
```python
from jaff import Network

# Load a KIDA format network
network = Network("networks/gas_reactions_kida.uva.2024.in")
print(f"Loaded {len(network.reactions)} reactions")
```

### PRIZMO Format  
```python
from jaff import Network

# PRIZMO networks use -> separator and VARIABLES{} blocks
network = Network("networks/react_COthin")

# Check what variables were defined
for reaction in network.reactions:
    print(f"Rate expression: {reaction.rate}")
```

### KROME Format
```python
from jaff import Network

# KROME networks use comma separation with @format headers
network = Network("networks/react_COthin") 

# KROME supports custom shortcuts like t32, te, etc.
print("Variables found:", [str(s) for s in network.reactions[0].rate.free_symbols])
```

## Species Analysis

### Basic Properties
```python
from jaff import Network

network = Network("networks/react_COthin")

# Find all carbon-bearing species
carbon_species = []
for species in network.species:
    if 'C' in species.exploded:
        carbon_species.append(species)
        
print(f"Found {len(carbon_species)} carbon-bearing species")
```

### Mass and Charge Distribution
```python
from jaff import Network
import matplotlib.pyplot as plt
import numpy as np

# Load network
network = Network("networks/react_COthin")

# Plot mass distribution
masses = [s.mass for s in network.species]
plt.hist(masses, bins=50)
plt.xlabel('Molecular Mass (amu)')
plt.ylabel('Number of Species')
plt.title('Species Mass Distribution')
plt.show()

# Plot charge distribution  
charges = [s.charge for s in network.species]
unique_charges, counts = np.unique(charges, return_counts=True)
plt.bar(unique_charges, counts)
plt.xlabel('Charge')
plt.ylabel('Number of Species')
plt.title('Species Charge Distribution')
plt.show()
```

## Reaction Analysis

### Rate Coefficient Tables
```python
from jaff import Network

# Load network
network = Network("networks/react_COthin")

# Generate rate table for temperature range
temp, rates = network.get_table(
    T_min=10,      # 10 K
    T_max=1000,    # 1000 K
    nT=100,        # 100 temperature points
    err_tol=0.01,  # 1% error tolerance
    verbose=True   # Show adaptive refinement
)

print(f"Rate table shape: {rates.shape}")
print(f"Temperature range: {temp[0]:.1f} - {temp[-1]:.1f} K")
```

### Temperature-Dependent Plotting
```python
from jaff import Network
import numpy as np
import matplotlib.pyplot as plt

# Load network
network = Network("networks/react_COthin")

# Get rates for specific reaction
reaction_idx = 2  # Use reaction index 2 since we only have 3 reactions
temp, rates = network.get_table(T_min=10, T_max=1000, nT=100)
reaction_rates = rates[reaction_idx, :]

plt.loglog(temp, reaction_rates)
plt.xlabel('Temperature (K)')
plt.ylabel('Rate Coefficient (cm³/s)')
plt.title(f'Reaction: {network.reactions[reaction_idx]}')
plt.show()
```

## Network Comparison

### Compare Two Networks
```python
from jaff import Network

# Load two different networks
network1 = Network("networks/react_COthin")
network2 = Network("networks/gas_reactions_kida.uva.2024.in")

# Compare reactions (verbosity=1 shows differences)
# network1.compare_reactions(network2, verbosity=1)

# Note: Species comparison can fail if networks have very different species sets
# network1.compare_species(network2, verbosity=1)
```

### Find Common Elements
```python
# Get species names from both networks
species1 = set(s.name for s in network1.species)
species2 = set(s.name for s in network2.species)

# Find overlaps
common_species = species1.intersection(species2)
unique_to_1 = species1.difference(species2)
unique_to_2 = species2.difference(species1)

print(f"Common species: {len(common_species)}")
print(f"Unique to network 1: {len(unique_to_1)}")
print(f"Unique to network 2: {len(unique_to_2)}")
```

## Advanced Usage

### Custom Rate Evaluation
```python
from sympy import symbols

# Get a reaction's symbolic rate expression
reaction = network.reactions[0]
rate_expr = reaction.rate

# Substitute specific values
tgas, av, crate = symbols('tgas av crate')
substituted_rate = rate_expr.subs([
    (tgas, 100),    # 100 K
    (av, 1.0),      # 1 mag extinction  
    (crate, 1.3e-17) # Standard cosmic ray rate
])

print(f"Rate at T=100K: {substituted_rate}")
```

### Export to Other Formats
```python
# Export species list to CSV
import csv

with open('species_list.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Name', 'Mass', 'Charge', 'Formula'])
    
    for species in network.species:
        writer.writerow([
            species.name,
            species.mass,
            species.charge,
            species.exploded
        ])
```

### Working with ODEs
```python
# The network automatically generates ODE structure
print(f"Reactant list shape: {network.rlist.shape}")
print(f"Product list shape: {network.plist.shape}")

# These arrays define the ODE system structure:
# rlist[i, :] = indices of reactants for reaction i
# plist[i, :] = indices of products for reaction i
```
