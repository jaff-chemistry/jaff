---
tags:
    - User-guide
    - Network
icon: phosphor/graph
---

# Working with Networks

Once a network is loaded, JAFF exposes a rich set of attributes and methods for inspecting species, reactions, and elements — and for exporting the network in various formats.

```python
from jaff import Network

net = Network("networks/h_photoionization/h_photo.jet")

print(f"Label:     {net.label}")
print(f"Species:   {net.species.count}")
print(f"Reactions: {net.reactions.count}")
print(f"Elements:  {net.elements.count}")
```

<div class="grid cards" markdown>

- :phosphor-flask:{ .sm .middle } **Network**

    ***

    Top-level `Network` class: constructor options, key attributes, export methods, and stoichiometry matrices.

    [:octicons-arrow-right-24: Network](network.md)

- :phosphor-sliders:{ .sm .middle } **Configuration**

    ***

    The `[network]` settings shared by `jaff.toml` and `jaffgen.toml` (one schema): reactions, temperature cutoffs, radiation, dust, thermodynamics, and how the layers stack.

    [:octicons-arrow-right-24: Configuration](../configuration/index.md)

- :phosphor-dna:{ .sm .middle } **Species**

    ***

    Work with individual `Specie` objects and the `Species` catalogue: lookup, filtering, mass, charge, and elemental composition.

    [:octicons-arrow-right-24: Species](species.md)

- :phosphor-atom:{ .sm .middle } **Elements**

    ***

    Query the `Elements` catalogue derived from all species: density matrices, truth matrices, atomic properties.

    [:octicons-arrow-right-24: Elements](elements.md)

- :phosphor-lightning:{ .sm .middle } **Reactions**

    ***

    Inspect individual `Reaction` objects and the `Reactions` catalogue: filtering, conservation checks, rate expressions, and code generation.

    [:octicons-arrow-right-24: Reactions](reactions.md)

- :phosphor-terminal-window:{ .sm .middle } **jaffx CLI**

    ***

    Quick network inspection and rate-coefficient export from the command line — without running the full code-generation pipeline.

    [:octicons-arrow-right-24: jaffx](jaffx.md)

</div>
