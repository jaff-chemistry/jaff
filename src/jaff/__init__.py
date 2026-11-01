# ABOUTME: JAFF package root: astrochemical network parser and multi-language code
# ABOUTME: generator; re-exports the public API (Network, Codegen, ...).
"""JAFF - Just Another Fancy Format

Astrochemical reaction network parser supporting KROME, PRIZMO, UDFA, KIDA, and UCLChem formats.
"""

__version__ = "0.1.5"

from .codegen.codegen import Codegen
from .core.elements import Element, Elements
from .core.network import Network
from .core.reaction import Reaction, Reactions
from .core.species import Specie, Species

__all__ = [
    "Element",
    "Elements",
    "Network",
    "Reaction",
    "Reactions",
    "Specie",
    "Species",
    "Codegen",
]
