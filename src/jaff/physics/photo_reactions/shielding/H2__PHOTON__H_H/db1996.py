# ABOUTME: Draine & Bertoldi (1996) H2 self-shielding (shielding.type = "db1996"),
# ABOUTME: bound to the H2 photodissociation reaction H2._PHOTON__H.H
"""
H2 shielding by Draine & Bertoldi 1996
DOI:https://ui.adsabs.harvard.edu/link_gateway/1996ApJ...468..269D/doi:10.1086/177689
"""

from __future__ import annotations

from typing import TYPE_CHECKING, cast

from sympy import Expr

from jaff.physics.photo_reactions.shielding.config import ShieldingParams

from .. import ShieldingFunction, _register
from ._utils import shielding

if TYPE_CHECKING:
    from jaff.core.network import Network
    from jaff.core.reaction import Reaction


class DB1996Params(ShieldingParams):
    """Settings of the Draine & Bertoldi (1996) H2 self-shielding.

    Attributes
    ----------
    min_ncol : float
        Floor on the dimensionless column ``x = N(H2) / 5e14 cm^-2`` of the fit
        (a floor on ``x``, not on the column density in cm^-2).
    min_vdisp : float
        Floor on the dimensionless Doppler parameter
        ``b5 = sqrt(2) * vdisp / 1e5 cm/s`` of the fit (a floor on ``b5``, not
        on the velocity dispersion in cm/s); keeps ``x / b5`` finite when
        ``vdisp`` is 0.
    """

    min_ncol: float = 1.0e-50
    min_vdisp: float = 1.0e-50


@_register
class DB1996(ShieldingFunction):
    """Draine & Bertoldi (1996) H2 self-shielding (``shielding.type = "db1996"``)."""

    name = "db1996"
    reaction = "H2._PHOTON__H.H"

    @staticmethod
    def config_model() -> type[DB1996Params]:
        """Return the settings model: the ``min_ncol`` / ``min_vdisp`` floors."""
        return DB1996Params

    def get_shielding(self, reaction: "Reaction", network: "Network") -> Expr:
        """Return the Draine & Bertoldi (1996) H2 self-shielding factor.

        Thin wrapper over :func:`shielding` with ``alpha = 2.0``, using the
        optional floors ``shielding.min_ncol`` (on the dimensionless column
        ``x``) and ``shielding.min_vdisp`` (on the dimensionless Doppler
        parameter ``b5``).

        Raises
        ------
        ParserError
            If a shielding setting is unknown or has an invalid value.
        """
        params = cast(DB1996Params, self.params(reaction))

        return shielding(
            alpha=2.0,
            min_ncol=params.min_ncol,
            min_vdisp=params.min_vdisp,
        )
