# ABOUTME: HAS's entity kinds: each kind templates may test and the collection it checks
# ABOUTME: GET uses the same table to check its entity; a collection must support `in`

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable, Dict

if TYPE_CHECKING:
    from ...context import RenderContext

# Entity kind -> getter of the network collection holding that kind
Getter = Callable[["RenderContext"], Any]
KINDS: Dict[str, Getter] = {
    "specie": lambda c: c.net.species,
    "reaction": lambda c: c.net.reactions,
    "element": lambda c: c.net.elements,
}
