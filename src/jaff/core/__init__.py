# ABOUTME: Core data model of JAFF: networks, species, reactions and elements.
# ABOUTME: Re-exports the public classes of the jaff.core subpackages.
"""Primary data-model classes for chemical reaction networks.

``jaff.core`` contains the classes for loading, parsing and representing
astrochemical reaction networks:

- :class:`Network` loads and manages a reaction network from file;
  :class:`Kinetics` (``net.kinetics``) holds its chemical rates of change (species
  ODEs, weighted density rates); :class:`OdeSystem` (``net.ode_system(...)``) is the
  state vector and rhs a solver integrates.
- :class:`Species` / :class:`Specie` hold the network's species, each with
  mass, charge and other attributes and helper methods.
- :class:`Reactions` / :class:`Reaction` hold the network's reactions, each
  with a symbolic rate expression and other attributes and helper methods.
- :class:`Elements` / :class:`Element` extract the chemical elements of a set
  of species and their collective atomic properties; each :class:`Element` is
  a record of :func:`~jaff.chemistry.periodic_table`.

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

from . import analysis, elements, network, reaction, species
from .analysis import (
    Analysis,
    AnalysisSpec,
    AnalysisTable,
    Densities,
    Entry,
    Functions,
    Parameters,
)
from .elements import Element, Elements
from .network import Kinetics, Network, NetworkSymbols, OdeSystem, StateVariable
from .reaction import (
    PhotoReaction,
    Rate,
    ReactionEnergy,
    ReactionRadiation,
    Reaction,
    Reactions,
)
from .species import Specie, Species

__all__ = [
    "analysis",
    "elements",
    "network",
    "reaction",
    "species",
    "Element",
    "Elements",
    "Network",
    "Reaction",
    "Reactions",
    "Specie",
    "Species",
    "NetworkSymbols",
    "Kinetics",
    "OdeSystem",
    "StateVariable",
    "PhotoReaction",
    "ReactionRadiation",
    "ReactionEnergy",
    "Rate",
    "Analysis",
    "AnalysisSpec",
    "AnalysisTable",
    "Entry",
    "Parameters",
    "Densities",
    "Functions",
]
