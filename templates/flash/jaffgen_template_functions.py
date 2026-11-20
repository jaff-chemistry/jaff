# ABOUTME: $$JAFF CALL functions of the flash template: the SPECIES lines of the unit
# ABOUTME: Config, refusing species whose FLASH names (lower-case j/k form) coincide

from __future__ import annotations

from jaff.codegen import jaffgen_function


@jaffgen_function
def flash_species(net) -> str:
    """One ``SPECIES <name>`` line per species, in network order.

    FLASH names are the lower-case j/k form (``+`` -> ``j``, ``-`` -> ``k``);
    FLASH ignores case, so species that differ only in case (``CO``, ``Co``)
    would become one variable: they are an error.
    """
    names = list(net.species.normalized_names())
    clashes = sorted({n for n in names if names.count(n) > 1})
    if clashes:
        raise ValueError(f"species share a FLASH name: {clashes}")

    return "\n".join(f"SPECIES {name}" for name in names)
