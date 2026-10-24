# ABOUTME: Dust: container of a network's dust physics (tabular grain cross
# ABOUTME: sections and photoelectric emission)
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

from sympy import Expr, Float

from ...drivers.pooch import download_dust
from .photoelectric_emission import PhotoelectricEmission
from .tabular import Tabular

if TYPE_CHECKING:
    from ...core.network import Network
    from ..photo_reactions import RadiationGroup
    from .config import DustConfig


class Dust:
    """Dust-physics container for a network.

    :class:`~jaff.core.Network` builds one from its dust config, given as the
    ``dust`` argument (a dict or a
    :class:`~jaff.physics.dust.config.DustConfig`) or as the ``[network.dust]``
    table of a config file, and stores it as ``net.dust`` (``None`` without
    either).  It groups the dust processes as sub-objects: photoelectric
    emission (:attr:`pe`) and the tabulated dust cross sections
    (:attr:`tabular`).  With a radiation field, the radiation equations
    subtract the dust attenuation (:meth:`attenuation`) from each band's
    moment equations:

    ```
    d(den[i])/dt   -= Zd * c * den[i]   * n_hnuc * <σ_u>_i
    d(rflux[i])/dt -= Zd * c * rflux[i] * n_hnuc * <σ_f>_i
    ```

    with ``<σ_u>_i`` and ``<σ_f>_i`` the per-H-nucleus cross sections of kinds
    ``config.u_reduction`` and ``config.f_reduction`` averaged over band *i*
    (``tabular.avg_cross_section_per_hnuc``); a kind of ``"none"`` drops the
    term.

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
    config : DustConfig
        The ``[network.dust]`` settings (read-only): ``rv`` (R_V: ``3.1``,
        ``4.0`` or ``5.5``) selecting the dust model, and the cross-section
        kinds attenuating the radiation density (``u_reduction``) and flux
        (``f_reduction``): ``"extinction"``, ``"absorption"``,
        ``"scattering"``, ``"transport"`` or ``"none"``.
    pe : PhotoelectricEmission
        Photoelectric-emission model, source of the ``chi_pe`` symbol
        substitution (the radiation in the photoelectric band scaled to the
        network's reference background field, Draine by default).
    tabular : Tabular
        Tabulated per-H-nucleus dust cross sections for the selected R_V.

    Examples
    --------
    >>> from jaff import Network
    >>> net = Network(
    ...     "h_photoionization",
    ...     radiation={"bands": [1, 13.6, 100, "inf"]},
    ...     dust={"rv": 3.1},
    ... )
    >>> net.dust.config.rv, net.dust.config.u_reduction, net.dust.config.f_reduction
    (3.1, 'absorption', 'transport')
    >>> float(net.dust.tabular.avg_cross_section_per_hnuc("absorption", (13.6, 100.0)))
    1.3480313134398857e-21
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
        self.config: DustConfig = config
        self.pe: PhotoelectricEmission = PhotoelectricEmission(network, config.pe)
        self.tabular: Tabular = Tabular(self)

    def attenuation(self, group: RadiationGroup) -> tuple[Expr, Expr]:
        """The radiation the dust removes from *group*, per unit time.

        ``Zd · c · den[i] · n_hnuc · <σ_u>_i`` from the band's density and
        ``Zd · c · rflux[i] · n_hnuc · <σ_f>_i`` from its flux, with ``<σ_u>_i``
        and ``<σ_f>_i`` the per-H-nucleus cross sections of the
        ``u_reduction`` and ``f_reduction`` kinds averaged over the band
        (``tabular.avg_cross_section_per_hnuc``); a kind of ``"none"`` gives
        ``0``.

        Parameters
        ----------
        group : RadiationGroup
            The radiation band (``net.radiation`` must be set).

        Returns
        -------
        tuple[sympy.Expr, sympy.Expr]
            ``(density_loss, flux_loss)``, to subtract from the band's density
            and flux equations.
        """
        radiation = self.net.radiation
        assert radiation is not None

        symbols = self.net.symbols
        base = symbols.zd * radiation.c * symbols.n_hnuc
        bounds = (group.lower, group.upper)

        def loss(kind: str, moment: Expr) -> Expr:
            if kind == "none":
                return Float(0.0)

            return base * moment * self.tabular.avg_cross_section_per_hnuc(kind, bounds)

        return (
            loss(self.config.u_reduction, radiation.den[group.index]),
            loss(self.config.f_reduction, radiation.flux[group.index]),
        )
