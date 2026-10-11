"""Leiden tabulated line-shielding (global shielding function).

Builds the shielding factor for a photo-reaction from the collapsed Leiden
shielding tables (``data/shielding/leiden.hdf5``, one group per reaction).  For
each requested shielding species the relevant column-density grid (``N``) and
shielding-factor column are extracted into a per-reaction
``shielding_<reaction>.hdf5`` next to the generated code, and the total factor
is the product of one interpolation call per shielding species.

A reaction selects this function via ``shielding.type = "leiden"`` and must
declare ``shielding.shielded_by`` (a subset of :data:`LEIDEN_SPECIES_MAP`); the
optional ``shielding.radiation`` field picks the radiation-field subgroup
(default ``"ISRF"``).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from sympy import Expr, Float, Function

from jaff.config import SHIELDING_DATA_DIR
from jaff.drivers import HDF5
from jaff.errors import ParserError

from .. import ShieldingFunction, _register

if TYPE_CHECKING:
    from jaff.core.network import Network
    from jaff.core.reaction import Reaction

LEIDEN_SPECIES_MAP: dict[str, str] = {
    "H2": "H2",
    "H": "H",
    "C": "C",
    "N2": "N2",
    "CO": "CO",
    "self": "",
}


@_register
class Leiden(ShieldingFunction):
    """Leiden tabulated line-shielding (``shielding.type = "leiden"``)."""

    name = "leiden"
    reaction = None  # global: applies to any reaction

    def get_shielding(self, reaction: "Reaction", network: "Network") -> Expr:
        """Return the Leiden tabulated shielding factor for *reaction*.

        Extracts the column-density grid and one shielding-factor column per
        species in ``reaction._metadata["shielding"]["shielded_by"]`` from
        ``data/shielding/leiden.hdf5``, writes them to a per-reaction
        ``shielding_<reaction>.hdf5`` in the generator output directory, and
        returns the product of one
        ``interp_<index>_shielding_<species>(ncol_<species>)`` interpolation
        call per shielding species.

        Recognised ``shielded_by`` entries are the keys of
        :data:`LEIDEN_SPECIES_MAP`; ``"self"`` resolves to the reaction's
        reactant (its own column density).  The optional ``shielding.radiation``
        field selects the radiation-field subgroup and defaults to ``"ISRF"``.

        Raises
        ------
        ParserError
            If ``shielded_by`` is missing, or names a species not in
            :data:`LEIDEN_SPECIES_MAP`.
        """
        sprops = reaction._metadata["shielding"]
        if "shielded_by" not in sprops:
            raise ParserError(
                "'shielded_by' must be specified for Leiden shielding tables"
            )

        if "radiation" not in sprops:
            sprops["radiation"] = "ISRF"

        # "self" shields by the reaction's own reactant column density.
        species_map = {**LEIDEN_SPECIES_MAP, "self": f"{reaction.reactants[0]}"}
        if any(sp not in species_map for sp in sprops["shielded_by"]):
            raise ParserError(
                f"Invalid shielding specie detected for reaction: {reaction}"
            )

        sprops["radiation"] = sprops["radiation"].lower()
        # The N grid is shared at the reaction-group level; the shielding-factor
        # columns live under the radiation-field subgroup.
        h5group_core = f"{SHIELDING_DATA_DIR}/leiden.hdf5::{reaction.serialized}"
        h5group_rad = f"{h5group_core}/{sprops['radiation']}"
        h5obj = HDF5()

        ncol = h5obj.to_dict(h5group_core, include="N")
        shielding_table = h5obj.to_dict(h5group_rad, include=sprops["shielded_by"])
        h5dict = {**ncol, **shielding_table}

        # Emit a per-reaction shielding table the generated code interpolates over.
        h5obj.from_dict(
            f"{reaction._metadata['jaffgen']['jaffgen_object'].jaffgen_config['output_dir']}/shielding_{reaction.serialized}.hdf5",
            h5dict,
            mode="w",
        )

        # Total factor = product of one interpolation call per shielding species.
        from jaff.core import NetworkSymbols

        shielding: Expr = Float(1.0)
        for specie in sprops["shielded_by"]:
            interp = Function(f"interp_{reaction.index}_shielding_{specie}")
            shielding *= interp(NetworkSymbols.ncol(species_map[specie]))

        return shielding
