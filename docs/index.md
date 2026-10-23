---
tags:
    - Introduction
icon: phosphor/house
---

# JAFF

<div align="left" markdown>
<div style="display: flex; align-items: center;">
    <img class="jaff-logo" src="assets/brand/logo.png" alt="JAFF Logo" width="60" style="margin-right: 20px;">
    <div>
    <strong>Just Another Fancy Format</strong>
    <br>
    <em>An astrochemical network parser and multi-language code generator</em>
    </div>
</div>

---

<!--## Welcome to JAFF-->

JAFF is a comprehensive library for working with astrochemical reaction networks. It provides an unified interface for parsing multiple network formats, analyzing chemical species and reactions, and generating optimized code for numerical simulations in multiple programming languages: `C`, `C++`, `Fortran`, `Python`, `Rust`, `Julia`, and `R` and is the first automated code generator capable of treating photo-chemistry explicitly.

<div class="grid cards" markdown>

- :phosphor-atom:{ .sm .middle } **Multi-format support**

    ***

    Parses `KIDA`, `UDFA`, `PRIZMO`, `KROME` and `UCLCHEM` networks with automatic format detection.

- :phosphor-chart-line:{ .sm .middle } **Validation and Analysis**

    ***

    Automatically validates mass and charge conservation, and extracts elemental composition.

- :phosphor-code:{ .sm .middle } **Code Generation**

    ***

    Generates optimized code for rates, ODEs, Jacobians, and fluxes `C`, `C++`, `Fortran`, `Python`, `Rust`, `Julia`, and `R` with optional CSE.

- :phosphor-chart-bar:{ .sm .middle } **Template System**

    ***

    Powerful template language (JAFF directives) for customizing code generation to match any simulation framework.

</div>

---

## Quick Example

```python
from jaff import Network

# Load a chemical network
net = Network("networks/h_photoionization/h_photo.jet")

# Access species information
print(f"Network contains {net.species.count} species")
print(f"First species: {net.species[0].name}, mass: {net.species[0].mass} gm")

# Access reactions
rea = net.reactions["H + _PHOTON -> H+ + e-"]
print(f"Reactants: {', '.join(rea.reactants.names())}")
print(f"Products: {', '.join(rea.products.names())}")
```

**Output**

```text
Network contains 3 species
First species: H, mass: 1.673773e-24 gm
Reactants: H, _PHOTON
Products: H+, e-
```

---

## Supported Network Formats

| Format      | Description                                          | Reference                                                            |
| ----------- | ---------------------------------------------------- | -------------------------------------------------------------------- |
| **KIDA**    | Kinetic Database for Astrochemistry                  | [A&A, 689, A63 (2024)](https://doi.org/10.1051/0004-6361/202450606)  |
| **UDFA**    | UMIST Database for Astrochemistry                    | [A&A, 682, A109 (2024)](https://doi.org/10.1051/0004-6361/202346908) |
| **PRIZMO**  | Protoplanetary disk chemistry evolution code         | [MNRAS 494, 4471 (2020)](https://doi.org/10.1093/mnras/staa971)      |
| **KROME**   | Library for astrophysical chemistry and microphysics | [MNRAS 439, 2386 (2014)](https://doi.org/10.1093/mnras/stu114)       |
| **UCLCHEM** | Gas-grain astrochemical code for Python modelling    | [AJ 154 38 (2017)](https://doi.org/10.3847/1538-3881/aa773f)         |

---

## What's Next?

<div class="grid cards" markdown>

- :phosphor-rocket-launch:{ .sm .middle } **Getting Started**

    ***

    New to JAFF? The installation documentation provides a comprehensive guide to help you get started with JAFF quickly

    [:octicons-arrow-right-24: Installation Guide](getting-started/installation.md)

- :phosphor-book-open:{ .sm .middle } **User Guide**

    ***

    The user guide provides in-depth information on the key concepts of Networks, Reactions and Species and how to utilize them

    [:octicons-arrow-right-24: User Guide](user-guide/designing-networks/index.md)

- :phosphor-brackets-curly:{ .sm .middle } **Code Generation**

    ***

    Use JAFF's templated code generation capabilites to generate code any any of the major programming languages to simulate chemical reactions

    [:octicons-arrow-right-24: Code Generation Guide](user-guide/code-generation/index.md)

- :phosphor-kanban:{ .sm .middle } **API Reference**

    ***

    The reference guide contains a detailed description of the functions, modules, and objects used by JAFF, assuming that you have an understanding of basic concepts

    [:octicons-arrow-right-24: API Docs](reference/index.md)

</div>

---

## Community & Support

- **GitHub Issues**: Report bugs and request features at [github.com/jaff-chemistry/jaff/issues](https://github.com/jaff-chemistry/jaff/issues)
- **Contributing**: See our [Contributing Guide](development/contributing.md)

---

## Citation

If you use JAFF in your research, please cite:

```bibtex
@software{jaff2024,
  title = {JAFF: Just Another Fancy Format},
  author = {JAFF Team},
  year = {2024},
  url = {https://github.com/jaff-chemistry/jaff},
  version = {0.1.5}
}
```

---

## License

JAFF is released under the [MIT License](about/license.md).

---

<!-- prettier-ignore -->
!!! tip "New to astrochemistry?"
    Check out our [Basic Concepts](getting-started/concepts.md) page to learn about chemical networks, reaction rates, and how JAFF can help your research.

<!-- prettier-ignore -->
!!! example "Ready to dive in?"
    Jump straight to the [Quick Start Guide](getting-started/quickstart.md) to start using JAFF in minutes!
