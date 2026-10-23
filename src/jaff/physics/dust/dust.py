"""Dust physics for a JAFF network.

:class:`Dust` is the container for the network's dust-related physics.  It is
attached to a :class:`~jaff.core.network.Network` when the network is built
from its ``[network.dust]`` config (:class:`~jaff.physics.dust.config.DustConfig`)
and groups the individual
dust processes as sub-objects (currently photoelectric emission and tabular
dust cross sections; more may be added later).

The symbols exposed by these sub-objects (e.g. the scaled photoelectric-band
radiation field ``chi_pe`` supplied by :attr:`Dust.pe`) are substituted into
reaction-rate expressions during network standardisation.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ...drivers.pooch import download_dust
from .photoelectric_emission import PhotoelectricEmission
from .tabular import Tabular

if TYPE_CHECKING:
    from ...core.network import Network
    from .config import DustConfig


class Dust:
    """Dust-physics container for a network.

    Parameters
    ----------
    network : Network
        The parent network this dust model belongs to.
    config : DustConfig
        The ``[network.dust]`` config: dust model selection, radiation
        reduction settings and photoelectric band.

    Attributes
    ----------
    net : Network
        Back-reference to the parent network.
    pe : PhotoelectricEmission
        Photoelectric-emission model, source of the ``chi_pe`` symbol
        substitution (the Draine field scaled to the photoelectric band).
    rv : float
        Total-to-selective extinction ratio R_V (from ``props``) selecting
        the WD01/Draine dust model.
    u_reduction : str
        Dust cross-section kind used to attenuate the radiation *density*
        (0th) moment (from ``props``).
    f_reduction : str
        Dust cross-section kind used to attenuate the radiation *flux*
        (1st) moment (from ``props``).
    tabular : Tabular
        Tabulated per-H-nucleus dust cross sections for the selected R_V.
    """

    def __init__(
        self,
        network: Network,
        config: DustConfig,
    ):
        """Build the dust model and its process sub-objects.

        Parameters
        ----------
        network : Network
            The parent network.
        config : DustConfig
            The ``[network.dust]`` config: the R_V dust-model selection
            (``rv``), the radiation reduction kinds (``u_reduction``,
            ``f_reduction``) and the photoelectric band (``pe``).
        """
        download_dust()

        self.net: Network = network
        self.pe: PhotoelectricEmission = PhotoelectricEmission(network, config.pe)
        self.rv: float = config.rv
        self.u_reduction: str = config.u_reduction
        self.f_reduction: str = config.f_reduction
        self.tabular: Tabular = Tabular(self)
