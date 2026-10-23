# ABOUTME: Leiden tabulated line shielding (shielding.type = "leiden"), a global
# ABOUTME: shielding function interpolating the shipped leiden.hdf5 tables
"""Leiden tabulated line-shielding (global shielding function).

Builds the shielding factor for a photo-reaction from the collapsed Leiden
shielding tables (``data/shielding/leiden.hdf5``, one group per reaction).  For
each requested shielding species the column-density grid (``N``) and the
shielding-factor column are extracted into a per-reaction table, returned with
the factor (kept as ``photo.shielding_tables``) and written out as
``shielding_<reaction>.hdf5`` by
:meth:`~jaff.core.network.Network.write_shielding_tables` (``jaffgen`` writes it
to the output directory).  The total factor is the product of one
interpolation call per shielding species.

A reaction selects this function via ``shielding.type = "leiden"`` and must
declare ``shielding.shielded_by`` (a subset of :data:`LeidenSpecies`); the
optional ``shielding.radiation`` field picks the radiation-field subgroup
(default ``"isrf"``).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, List, Literal, cast

from sympy import Expr, Float, Function

from jaff.config import SHIELDING_DATA_DIR
from jaff.config.base import LowerStr
from jaff.physics.photo_reactions.shielding.config import ShieldingParams
from jaff.drivers import HDF5

from .. import ShieldingFunction, ShieldingResult, register

if TYPE_CHECKING:
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
        Radiation field the shielding factors were computed for: the name of a
        subgroup of the reaction's group in the shipped ``leiden.hdf5``
        (case-insensitive).  Every reaction has ``"isrf"`` (the default),
        ``"bb-4000"``, ``"bb-10000"``, ``"bb-20000"``, ``"ly-alpha"``,
        ``"solar"`` and ``"tw-hya"``; all but ``"Zn._PHOTON__Zn+.e-"`` have
        ``"mathis-1983"``; only ``"HF._PHOTON__F.H"`` and
        ``"HF._PHOTON__HF+.e-"`` have ``"gondhalekar-1980"`` and
        ``"habing-1968"``.  Some subgroups lack some species (``"ly-alpha"``
        is empty for 84 of the 149 reactions).  The value is not checked when
        the config is read: a missing subgroup fails when the shielding table
        is read from ``leiden.hdf5``.
    """

    shielded_by: List[LeidenSpecies]
    radiation: LowerStr = "isrf"


@register
class Leiden(ShieldingFunction):
    """Leiden tabulated line-shielding (``shielding.type = "leiden"``)."""

    name = "leiden"
    reaction = None  # global: applies to any reaction

    @staticmethod
    def config_model() -> type[LeidenParams]:
        """Return the settings model: ``shielded_by`` and ``radiation``."""
        return LeidenParams

    def get_shielding(self, reaction: "Reaction") -> ShieldingResult:
        """Return the Leiden tabulated shielding factor for *reaction* and its table.

        Extracts the column-density grid and one shielding-factor column per
        species in the ``shielded_by`` setting from
        ``data/shielding/leiden.hdf5`` as the reaction's
        ``shielding_<reaction>.hdf5`` table, and returns it with the product of
        one
        ``interp_<index>_shielding_<species>(ncol_<species>)`` interpolation
        call per shielding species.

        Recognised ``shielded_by`` entries are :data:`LeidenSpecies`; ``"self"``
        resolves to the reaction's reactant (its own column density).  The
        optional ``radiation`` setting selects the radiation-field subgroup and
        defaults to ``"isrf"``.

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

        # The per-reaction table the generated code interpolates over.
        tables = {f"shielding_{reaction.serialized}.hdf5": h5dict}

        # Total factor = product of one interpolation call per shielding species.
        from jaff.core import NetworkSymbols

        shielding: Expr = Float(1.0)
        for specie in params.shielded_by:
            # "self" shields by the reaction's own reactant column density.
            column = f"{reaction.reactants[0]}" if specie == "self" else specie
            interp = Function(f"interp_{reaction.index}_shielding_{specie}")
            shielding *= interp(NetworkSymbols.ncol(column))

        return ShieldingResult(shielding, tables)
