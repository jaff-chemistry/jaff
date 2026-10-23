"""
H2 shielding by Draine & Bertoldi 1996
DOI:https://ui.adsabs.harvard.edu/link_gateway/1996ApJ...468..269D/doi:10.1086/177689
"""

from __future__ import annotations

from typing import TYPE_CHECKING, cast

from sympy import Expr

from jaff.config.schema import ShieldingParams

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
        Lower floor on the H2 column density (cm^-2).
    min_vdisp : float
        Lower floor on the velocity dispersion (cm/s).
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
        optional floors ``shielding.min_ncol`` (cm^-2) and
        ``shielding.min_vdisp`` (cm/s).

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
