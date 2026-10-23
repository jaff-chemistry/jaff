"""Leiden tabulated line-shielding (global shielding function).

Builds the shielding factor for a photo-reaction from the collapsed Leiden
shielding tables (``data/shielding/leiden.hdf5``, one group per reaction).  For
each requested shielding species the relevant column-density grid (``N``) and
shielding-factor column are extracted into a per-reaction
``shielding_<reaction>.hdf5`` next to the generated code, and the total factor
is the product of one interpolation call per shielding species.

A reaction selects this function via ``shielding.type = "leiden"`` and must
declare ``shielding.shielded_by`` (a subset of :data:`LeidenSpecies`); the
optional ``shielding.radiation`` field picks the radiation-field subgroup
(default ``"ISRF"``).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, List, Literal, cast

from sympy import Expr, Float, Function

from jaff.config import SHIELDING_DATA_DIR
from jaff.config.base import LowerStr
from jaff.physics.photo_reactions.shielding.config import ShieldingParams
from jaff.drivers import HDF5

from .. import ShieldingFunction, _register

if TYPE_CHECKING:
    from jaff.core.network import Network
    from jaff.core.reaction import Reaction

#: Species with a shielding column in the Leiden tables; ``"self"`` stands for
#: the reaction's own reactant.
LeidenSpecies = Literal["H2", "H", "C", "N2", "CO", "self"]


class LeidenParams(ShieldingParams):
    """Settings of the Leiden tabulated shielding.

    Attributes
    ----------
    shielded_by : list of str
        Shielding species, ``"self"`` being the reaction's own reactant; their
        factors multiply.
    radiation : str
        Radiation-field subgroup of the Leiden table (case-insensitive).
    """

    shielded_by: List[LeidenSpecies]
    radiation: LowerStr = "isrf"


@_register
class Leiden(ShieldingFunction):
    """Leiden tabulated line-shielding (``shielding.type = "leiden"``)."""

    name = "leiden"
    reaction = None  # global: applies to any reaction

    @staticmethod
    def config_model() -> type[LeidenParams]:
        """Return the settings model: ``shielded_by`` and ``radiation``."""
        return LeidenParams

    def get_shielding(self, reaction: "Reaction", network: "Network") -> Expr:
        """Return the Leiden tabulated shielding factor for *reaction*.

        Extracts the column-density grid and one shielding-factor column per
        species in the ``shielded_by`` setting from
        ``data/shielding/leiden.hdf5``, keeps them on the reaction as its
        ``shielding_<reaction>.hdf5`` table (``photo.shielding.tables``), and
        returns the product of one
        ``interp_<index>_shielding_<species>(ncol_<species>)`` interpolation
        call per shielding species.

        Recognised ``shielded_by`` entries are :data:`LeidenSpecies`; ``"self"``
        resolves to the reaction's reactant (its own column density).  The
        optional ``radiation`` setting selects the radiation-field subgroup and
        defaults to ``"ISRF"``.

        Raises
        ------
        ParserError
            If ``shielded_by`` is missing, names a species not in
            :data:`LeidenSpecies`, or another setting is unknown or invalid.
        """
        params = cast(LeidenParams, self.params(reaction))

        # The N grid is shared at the reaction-group level; the shielding-factor
        # columns live under the radiation-field subgroup.
        h5group_core = f"{SHIELDING_DATA_DIR}/leiden.hdf5::{reaction.serialized}"
        h5group_rad = f"{h5group_core}/{params.radiation}"
        h5obj = HDF5()

        ncol = h5obj.to_dict(h5group_core, include="N")
        shielding_table = h5obj.to_dict(
            h5group_rad, include=[str(sp) for sp in params.shielded_by]
        )
        h5dict = {**ncol, **shielding_table}

        # Keep the per-reaction table the generated code interpolates over.
        assert reaction.photo is not None and reaction.photo.shielding is not None
        tables = reaction.photo.shielding.tables
        tables[f"shielding_{reaction.serialized}.hdf5"] = h5dict

        # Total factor = product of one interpolation call per shielding species.
        from jaff.core import NetworkSymbols

        shielding: Expr = Float(1.0)
        for specie in params.shielded_by:
            # "self" shields by the reaction's own reactant column density.
            column = f"{reaction.reactants[0]}" if specie == "self" else specie
            interp = Function(f"interp_{reaction.index}_shielding_{specie}")
            shielding *= interp(NetworkSymbols.ncol(column))

        return shielding
