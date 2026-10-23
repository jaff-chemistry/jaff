# ABOUTME: Hartwig et al. (2015) H2 self-shielding (shielding.type = "hg2015"),
# ABOUTME: bound to the H2 photodissociation reaction H2._PHOTON__H.H
"""
H2 shielding by Hartwig et.al. 2015
DOI:https://doi.org/10.1093/mnras/stv1368
"""

from __future__ import annotations

from typing import TYPE_CHECKING, cast

from sympy import Expr

from jaff.physics.photo_reactions.shielding.config import ShieldingParams

from jaff.physics.photo_reactions.shielding import _register
from jaff.physics.photo_reactions.shielding._base import ShieldingFunction

from ._utils import shielding

if TYPE_CHECKING:
    from jaff.core.network import Network
    from jaff.core.reaction import Reaction


class HG2015Params(ShieldingParams):
    """Settings of the Hartwig et al. (2015) H2 self-shielding.

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
class HG2015(ShieldingFunction):
    """Hartwig et al. (2015) H2 self-shielding (``shielding.type = "hg2015"``)."""

    name = "hg2015"
    reaction = "H2._PHOTON__H.H"

    @staticmethod
    def config_model() -> type[HG2015Params]:
        """Return the settings model: the ``min_ncol`` / ``min_vdisp`` floors."""
        return HG2015Params

    def get_shielding(self, reaction: "Reaction", network: "Network") -> Expr:
        """Return the Hartwig et al. (2015) H2 self-shielding factor.

        Thin wrapper over :func:`shielding` with ``alpha = 1.1``, using the
        optional floors ``shielding.min_ncol`` (on the dimensionless column
        ``x``) and ``shielding.min_vdisp`` (on the dimensionless Doppler
        parameter ``b5``).

        Raises
        ------
        ParserError
            If a shielding setting is unknown or has an invalid value.
        """
        params = cast(HG2015Params, self.params(reaction))

        return shielding(
            alpha=1.1,
            min_ncol=params.min_ncol,
            min_vdisp=params.min_vdisp,
        )
