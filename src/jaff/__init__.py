# ABOUTME: JAFF package root: astrochemical network parser and multi-language code
# ABOUTME: generator; re-exports the public API (Network, Codegen, Builder, ...).
"""JAFF - Just Another Fancy Format

Astrochemical reaction network parser supporting KROME, PRIZMO, UDFA, KIDA, and UCLChem formats.
"""

__version__ = "0.1.5"

from .codegen.builder import Builder
from .codegen.codegen import Codegen
from .codegen.preprocessor import Preprocessor
from .core.elements import Element, Elements
from .core.network import Network, NetworkSpec
from .core.reaction import Reaction, Reactions
from .core.species import Specie, Species

__all__ = [
    "Element",
    "Elements",
    "Network",
    "NetworkSpec",
    "Reaction",
    "Reactions",
    "Specie",
    "Species",
    "Builder",
    "Codegen",
    "Preprocessor",
]
