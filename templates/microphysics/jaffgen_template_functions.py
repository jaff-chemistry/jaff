# ABOUTME: $$JAFF CALL functions of the microphysics template: the species table of
# ABOUTME: pynucastro.net (name, enum-safe short name, A, Z), one line per species

from __future__ import annotations

from jaff.codegen import jaffgen_function
from jaff.physics.constants import u


#: Short names that clash with Microphysics identifiers, and their replacement
#: (``C`` is the constants namespace; Microphysics' metal_chem also uses Carbon)
RESERVED = {"C": "Carbon"}


@jaffgen_function
def species_table(net) -> str:
    """One ``name short_name A Z`` line per species, in network order.

    The short name is the j/k form (``+`` -> ``j``, ``-`` -> ``k``) capitalised as
    Microphysics' enum does, with :data:`RESERVED` names replaced; A is the mass
    in atomic mass units, Z the charge.
    """
    amu = u.cgs.value
    shorts = [
        RESERVED.get(short, short)
        for short in (name.capitalize() for name in net.species.normalized_names())
    ]
    clashes = sorted({s for s in shorts if shorts.count(s) > 1})
    if clashes:
        raise ValueError(f"species share a Microphysics short name: {clashes}")

    width = max(len(s.name) for s in net.species) + 2
    return "\n".join(
        f"{s.name:<{width}} {short:<{width}} {s.mass / amu:.10g}  {s.charge}"
        for s, short in zip(net.species, shorts)
    )
