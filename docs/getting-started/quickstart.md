---
tags:
    - Introduction
icon: phosphor/rocket-launch
---

# Quick Start Guide

JAFF does two things: it **loads a chemical reaction network** and lets you
inspect it, and it **generates solver code** from that network in the language
of your choice. This guide is a five-minute tour of both — just enough to see
the shape of the tool. Every step links to the [User Guide](../user-guide/working-with-networks/index.md)
where the same ground is covered in detail.

## Load and explore a network

A [`Network`](../user-guide/working-with-networks/network.md) is the backbone of
JAFF — it holds every species and reaction, with all their derived properties.
Load one from a file:

```python
from jaff import Network

net = Network("networks/h_photoionization/h_photo.jet")

print(f"Label:     {net.label}")
print(f"Species:   {net.species.count}")
print(f"Reactions: {net.reactions.count}")
```

```text
Label:     h_photo
Species:   3
Reactions: 2
```

JAFF reads all the common community formats — KIDA, UDFA, PRIZMO, KROME, and
UCLCHEM — and detects which one a file uses automatically, so you never pass a
format flag. See [Network Formats](../user-guide/designing-networks/network-formats.md).

### Species

`net.species` is a catalogue of every [`Specie`](../user-guide/working-with-networks/species.md);
iterate it to see each one's properties:

```python
for s in net.species:
    print(f"{s.index}: {s.name:<3} mass={s.mass:.5e} g  charge={s.charge:+d}")
```

```text
0: H   mass=1.67377e-24 g  charge=+0
1: H+  mass=1.67377e-24 g  charge=+1
2: e-  mass=9.10938e-28 g  charge=-1
```

### Reactions

`net.reactions` holds every [`Reaction`](../user-guide/working-with-networks/reactions.md);
`verbatim` is the human-readable form:

```python
for r in net.reactions:
    print(r.verbatim)
```

```text
H -> H+ + e-
H+ + e- -> H
```

<!-- prettier-ignore -->
!!! tip "Look up by key instead of looping"
    Indexing a catalogue returns the object directly — no search needed:
    `#!python net.species["H"]` gives the `H` specie, and `#!python net.reactions[0]`
    the first reaction.

The [Working with Networks](../user-guide/working-with-networks/index.md) guide
covers lookup, filtering, elemental composition, conservation checks, and export.

## Generate code

Code generation is primarily **template-driven**: JAFF expands ordinary source files that
contain `$$JAFF` directives into network-specific code. The fastest way to see it
work is a **built-in template** — a ready-made collection you don't have to
write — run through the [`jaffgen`](../user-guide/code-generation/jaffgen.md) CLI:

```bash
jaffgen \
    --network  networks/h_photoionization/h_photo.jet \
    --template fortran_dlsodes \
    --outdir   generated/
```

```text
INFO     Network loaded successfully!
INFO     Successfully generated files
```

`generated/` now holds a complete, buildable Fortran solver:

```text
commons.f90  fluxes.f90  ode.f90  reactions.f90  main.f90  Makefile
opkda1.f     opkda2.f    opkdmain.f
```

Each template names its own language with a `$$JAFF LANG` line or by its file
extension; see [jaffgen](../user-guide/code-generation/jaffgen.md) for the other
built-in templates and input options.

From here the [Code Generation](../user-guide/code-generation/index.md) guide
shows how to write your own templates with the
[directive language](../user-guide/code-generation/template-syntax.md). The
[`Builder`](../user-guide/advanced-code-generation/builder.md) API does the same
generation but is intended for manual code generation when the template syntax is
insufficient.

## Next steps

1. **Concepts** — [Basic Concepts](concepts.md) explains chemical networks and the JAFF model.
2. **Working with Networks** — [inspect, filter, and export](../user-guide/working-with-networks/index.md) loaded networks.
3. **Code Generation** — the [template directive language](../user-guide/code-generation/template-syntax.md) and the [`jaffgen`](../user-guide/code-generation/jaffgen.md) CLI.
4. **Examples** — step-by-step [example notebooks](examples.md) from loading a network to running generated solvers.
5. **API Reference** — the complete [API documentation](../api/index.md).

## Example networks

JAFF ships with several networks to experiment with:

- `networks/demos` — small test networks
- `networks/h_photoionization` — hydrogen photo-ionization (used above)
- `networks/COthin` — CO chemistry
- `networks/GOW` — the Gong–Ostriker–Wolfire network
