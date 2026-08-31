"""Dust physics for a JAFF network.

:class:`Dust` is the container for the network's dust-related physics.  It is
attached to a :class:`~jaff.core.network.Network` when the network is built
with ``dust=True`` and groups the individual dust processes as sub-objects
(currently photoelectric emission; more may be added later).

The symbols exposed by these sub-objects (e.g. the scaled photoelectric-band
radiation field ``chi_pe`` supplied by :attr:`Dust.pe`) are substituted into
reaction-rate expressions during network standardisation.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from .photoelectric_emission import PhotoelectricEmission

if TYPE_CHECKING:
    from ...core.network import Network


class Dust:
    """Dust-physics container for a network.

    Parameters
    ----------
    network : Network
        The parent network this dust model belongs to.

    Attributes
    ----------
    network : Network
        Back-reference to the parent network.
    pe : PhotoelectricEmission
        Photoelectric-emission model, source of the ``chi_pe`` symbol
        substitution (the Draine field scaled to the photoelectric band).
    """

    def __init__(self, network: Network):
        """Build the dust model and its process sub-objects.

        Parameters
        ----------
        network : Network
            The parent network.
        """
        self.network: Network = network
        self.pe: PhotoelectricEmission = PhotoelectricEmission(network)
