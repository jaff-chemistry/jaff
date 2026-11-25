<div align="center">
  <img src="assets/logo.png" alt="JAFF logo" width="120" />

  <h1>JAFF</h1>

  <p><em>Just Another Fancy Format</em></p>

  <p>A fast, multi-format astrochemical network parser with analysis, code generation, and explicit photochemistry.</p>

  <p>
    <a href="https://github.com/jaff-chemistry/jaff/blob/main/LICENSE"><img alt="License: MIT" src="https://img.shields.io/badge/License-MIT-8b6cff?style=for-the-badge&labelColor=241b2f&logo=opensourceinitiative&logoColor=EDE4FF"></a>
    <a href="https://www.python.org/downloads/"><img alt="Python 3.11+" src="https://img.shields.io/badge/Python-3.11%2B-e05fb0?style=for-the-badge&labelColor=241b2f&logo=python&logoColor=EDE4FF"></a>
    <img alt="Version" src="https://img.shields.io/badge/Version-0.1.5-ff6a5a?style=for-the-badge&labelColor=241b2f&logo=semanticrelease&logoColor=EDE4FF">
    <img alt="Status" src="https://img.shields.io/badge/Status-alpha-ffc24b?style=for-the-badge&labelColor=241b2f&logo=rocket&logoColor=EDE4FF">
    <a href="https://jaff-chemistry.github.io/jaff/"><img alt="Docs" src="https://img.shields.io/badge/Docs-online-8b6cff?style=for-the-badge&labelColor=241b2f&logo=readthedocs&logoColor=EDE4FF"></a>
  </p>

  <p>
    <a href="#installation">Installation</a> ·
    <a href="#quick-start">Quick Start</a> ·
    <a href="#supported-network-formats">Formats</a> ·
    <a href="https://jaff-chemistry.github.io/jaff/">Documentation</a> ·
    <a href="#contributing">Contributing</a>
  </p>
</div>

---

## Overview

**JAFF** is an astrochemical network parser that reads multiple network formats, validates and analyses them, generates simulation code, and handles explicit photochemistry — all from a single tool.

- **Multi-format** — parse KIDA, UDFA, PRIZMO, KROME, and UCLCHEM networks
- **Automatic validation** — catch malformed reactions and inconsistent species
- **Analysis** — inspect species, reactions, masses, charges, and rates
- **Code generation** — emit ODE solvers in C, C++, Python, Fortran, Rust, Julia, and R
- **Explicit photochemistry** — first-class treatment of photoreactions

> Full guides and API reference live in the [Documentation](https://jaff-chemistry.github.io/jaff/).

## Installation

Clone the repository and install with `pip`:

```bash
git clone https://github.com/jaff-chemistry/jaff.git
cd jaff
pip install .
```

Requires **Python 3.11+**.

## Quick Start

### Command Line

Installation provides two commands: `jaffx` for quick network inspection and
export, and `jaffgen` for code generation.

```bash
# Count species and reactions in a network
jaffx get num-species  --network networks/COthin/react_COthin.jet
jaffx get num-reactions --network networks/COthin/react_COthin.jet

# Export a table of rate coefficients over a temperature range to a text file
jaffx export txt --network networks/GOW/GOW.jet --tmin 10 --tmax 1e4 --file rates.txt
```

### Network Parsing

The format is auto-detected — the same `Network` call reads KIDA, KROME,
PRIZMO, UDFA, and UCLCHEM files, plus JAFF's own `.jaff` format:

```python
from jaff import Network

network = Network("networks/COthin/react_COthin.jet")
network = Network("networks/kida_uva_2024/gas_reactions_kida.uva.2024.jet")
```

### Analysis

```python
# Access species
for species in network.species:
    print(f"{species.name}: mass={species.mass}, charge={species.charge}")

# Access reactions
for reaction in network.reactions:
    print(f"{reaction.reactants}")
```

### Code Generation

```bash
jaffgen --template microphysics --network networks/GOW/GOW.jet
```

**Supported languages:** C · C++ · Python · Fortran · Rust · Julia · R

## Supported Network Formats

| Format      | Reference                                                                    |
| ----------- | ---------------------------------------------------------------------------- |
| **KIDA**    | [A&A, 689, A63 (2024)](https://doi.org/10.1051/0004-6361/202450606)          |
| **UDFA**    | [A&A, 682, A109 (2024)](https://doi.org/10.1051/0004-6361/202346908)         |
| **PRIZMO**  | [MNRAS 494, 4471–4491 (2020)](https://doi.org/10.1093/mnras/staa971)         |
| **KROME**   | [MNRAS 439, 2386–2419 (2014)](https://doi.org/10.1093/mnras/stu114)          |
| **UCLCHEM** | [J. Holdship et al 2017 AJ 154 38](https://doi.org/10.3847/1538-3881/aa773f) |

## Examples

Example network files can be found in the [`networks/`](networks/) directory.

## Contributing

Contributions are welcome! To contribute or modify JAFF, please refer to our
[Contributing Guide](CONTRIBUTING.md).

## License

Distributed under the **MIT License**. See [`LICENSE`](LICENSE) for details.

---

<div align="center">
  <img src="./assets/xkcd.png" alt="xkcd 927: Standards" width="500" /><br>
  <sub><a href="https://xkcd.com/927/">xkcd 927</a></sub>
</div>
