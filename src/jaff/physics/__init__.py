# ABOUTME: Physics package: symbolic ODE/flux generators, radiation, dust and EOS
# ABOUTME: models re-exported for Network and codegen
"""Physics models behind a network's ODEs: constants, radiation, dust and EOS.

A :class:`~jaff.core.network.Network` builds these objects for itself and
exposes them as attributes; they are rarely constructed directly:

* ``net.thermodynamics`` -- a :class:`Thermodynamics`: the equation of state
  (built by :class:`EosFactory` from ``[network.thermodynamics.eos]``), the
  heating/cooling rates ``dE/dt`` and the temperature rates ``dT/dt``.
* ``net.radiation`` -- a :class:`Radiation` of :class:`RadiationGroup` bands,
  present when a radiation field is configured (``[network.radiation]``).
* ``net.dust`` -- a :class:`Dust`, present when dust is configured
  (``[network.dust]``).
* :class:`Photochemistry` resolves photo-reaction cross sections
  (:class:`XsecProps`) and shielding factors from the bundled databases.
* :func:`get_sfluxes`, :func:`get_sodes` and :func:`get_sradodes` assemble the
  symbolic reaction fluxes, species ODEs and radiation-moment ODEs.

Notes
-----
:mod:`jaff.physics.constants` exposes physical constants as module-level
:class:`astropy.units.Quantity` objects.  There are no separate ``cgs`` /
``si`` / ``gaussian`` / ``natural`` tables: select the unit system per
quantity at the call site with ``.cgs``, ``.si`` or ``.to(...)``.

Examples
--------
>>> from jaff import Network
>>> from jaff.physics import constants as c
>>> c.k_B.cgs  # Boltzmann constant [erg/K]
<Quantity 1.380649e-16 erg / K>
>>> c.k_B.to("eV / K")
<Quantity 8.61733326e-05 eV / K>
>>> float(c.m_e.cgs.value)  # bare float [g] when a number is required
9.1093837139e-28
>>> net = Network("COthin", thermodynamics={"eos": {"type": "ideal", "gamma": 1.4}})
>>> net.thermodynamics.eos.per_particle.simplify()  # k_B T / (gamma - 1) [erg]
3.4516225e-16*tgas
"""

from . import constants
from ._equations import get_sfluxes, get_sodes, get_sradodes
from .dust import Dust
from .photo_reactions._photo_band import PhotoBand
from .photo_reactions._photochemistry import Photochemistry, XsecProps
from .photo_reactions._radiation import Radiation, RadiationGroup
from .thermodynamics import Thermodynamics
from .thermodynamics.eos import EosFactory

__all__ = [
    "XsecProps",
    "constants",
    "Photochemistry",
    "get_sfluxes",
    "get_sodes",
    "get_sradodes",
    "Radiation",
    "RadiationGroup",
    "PhotoBand",
    "Dust",
    "EosFactory",
    "Thermodynamics",
]
