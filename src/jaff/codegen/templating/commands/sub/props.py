# ABOUTME: SUB's tokens: each name templates may substitute and how its value is computed
# ABOUTME: Add a token as one entry; the SUB command validates and renders from this table

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable, Dict

if TYPE_CHECKING:
    from ...context import RenderContext

# Token name -> getter of its value from the render context
Getter = Callable[["RenderContext"], Any]
TOKENS: Dict[str, Getter] = {
    "nspec": lambda c: c.net.species.count,
    "nelem": lambda c: c.net.elements.count,
    "nreact": lambda c: c.net.reactions.count,
    "nbands": lambda c: c.net.radiation.nbands if c.net.radiation else 0,
    "label": lambda c: c.net.label,
    "filename": lambda c: c.file.name,
    "filepath": lambda c: c.file,
    "dedt": lambda c: c.cg.get_dedt(),
    "dtdt": lambda c: c.cg.get_dtdt(),
    "e_idx": lambda c: c.net.species.e_idx(),
}
