# ABOUTME: Core data model of JAFF: networks, species, reactions and elements.
# ABOUTME: Re-exports the public classes and typed dicts of the jaff.core subpackages.
"""Primary data-model classes for chemical reaction networks.

``jaff.core`` contains the classes for loading, parsing and representing
astrochemical reaction networks:

- :class:`Network` loads and manages a reaction network from file.
- :class:`Species` / :class:`Specie` hold the network's species, each with
  mass, charge and other attributes and helper methods.
- :class:`Reactions` / :class:`Reaction` hold the network's reactions, each
  with a symbolic rate expression and other attributes and helper methods.
- :class:`Elements` / :class:`Element` extract the chemical elements of a set
  of species and their collective atomic properties.

Examples
--------
Load a bundled network by name and inspect it:

>>> from jaff import Network
>>> net = Network("COthin")
>>> net.species.count, net.reactions.count
(37, 238)
>>> net.species["CO"].charge
0
>>> net.elements.symbols()
['C', 'H', 'He', 'O', 'Si']

Elements can also be built directly from a list of species names:

>>> from jaff.core import Elements
>>> elements = Elements(["CO", "H2O"])
>>> elements.symbols()
['C', 'H', 'O']
>>> elements.density_matrix()
[[1, 0], [0, 2], [1, 1]]
"""

from jaff.core._typing._reaction import ReactionProps

from . import elements, network, reaction, species
from ._typing import ElementProps
from .elements import Element, Elements
from .network import Network, NetworkSpec, NetworkSymbols
from .reaction import (
    PhotoReaction,
    Reaction,
    Reactions,
)
from .species import Specie, Species

__all__ = [
    "elements",
    "network",
    "reaction",
    "species",
    "Element",
    "Elements",
    "Network",
    "NetworkSpec",
    "Reaction",
    "Reactions",
    "Specie",
    "Species",
    "ElementProps",
    "ReactionProps",
    "NetworkSymbols",
    "PhotoReaction",
]
