# ABOUTME: GET's properties: name -> (entity kind, getter of the value for one entity)
# ABOUTME: Add a property here; the GET command validates and renders from this table

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable, Dict, Tuple

if TYPE_CHECKING:
    from ...context import RenderContext

# Property name -> (entity kind, getter taking the context and the entity)
_Getter = Callable[["RenderContext", str], Any]
PROPERTIES: Dict[str, Tuple[str, _Getter]] = {
    "element_idx": ("element", lambda c, e: c.net.elements.symbols().index(e)),
    "specie_idx": ("specie", lambda c, s: c.net.species[s].index),
    "reaction_idx": ("reaction", lambda c, r: c.net.reactions[r].catalogue_index),
    "specie_mass": ("specie", lambda c, s: c.net.species[s].mass),
    "specie_charge": ("specie", lambda c, s: c.net.species[s].charge),
    "specie_latex": ("specie", lambda c, s: c.net.species[s].latex()),
    "reaction_tmin": ("reaction", lambda c, r: c.net.reactions[r].tmin),
    "reaction_tmax": ("reaction", lambda c, r: c.net.reactions[r].tmax),
    "reaction_verbatim": ("reaction", lambda c, r: c.net.reactions[r].verbatim),
    "aux_func": (
        "aux",
        lambda c, f: c.cg.lang.code_gen(
            c.net.spec.aux_funcs[f].expr, strict=False, allow_unknown_functions=True
        ),
    ),
}
PROPS: Tuple[str, ...] = tuple(PROPERTIES)
