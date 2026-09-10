"""
Look-ups for photochemical cross sections.

Two backends are exposed:

- :func:`get_verner_xsec` -- analytic Verner (1996) photoionisation fits
  stored as SymPy strings in the ``verner_cross_sections`` SQLite table.
- :func:`get_xsec` -- tabulated Leiden / NORAD cross sections, indexed by the
  ``photo_reaction_cross_sections`` table and read from the corresponding
  HDF5 group as numpy arrays.

Both are keyed by ``reaction.serialized``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from sympy import Basic, Expr, sympify

from ...config import JAFF_DIR
from ...drivers import HDF5, JaffDb
from ...drivers.pooch import (
    download_background_radiation,
    download_shielding,
    download_xsecs,
)
from ._typing import XsecsProps
from .shielding import _get_shielding_function

if TYPE_CHECKING:
    from ...core import Reaction
    from ...core.network import Network


class Photochemistry:
    """Cross-section provider for photo-reactions.

    Resolves per-reaction photoionisation / photodissociation / photoabsorption
    cross sections from the JAFF database and the tabulated Leiden, NORAD, and
    Verner data files.
    """

    def __init__(self):
        """Ensure the cross-section and shielding data files are available locally.

        Constructing a :class:`Photochemistry` triggers
        :func:`~jaff.drivers.pooch.download_xsecs` (Leiden / NORAD / Verner
        cross sections) and :func:`~jaff.drivers.pooch.download_shielding` (the
        Leiden line-shielding tables), downloading both on first use (a network
        fetch unless already cached). Instantiate once and reuse rather than per
        reaction.
        """
        download_xsecs()
        download_shielding()

    def get_verner_xsec(self, reaction: Reaction) -> Basic | None:
        """
        Query the JAFF database for the Verner photoionisation cross section.

        Verner cross sections are analytic fits to photoionisation cross
        sections from Verner et al. (1996) stored as SymPy-parseable strings
        in the ``verner_cross_sections`` SQLite table.

        Parameters
        ----------
        reaction : Reaction
            Reaction whose serialised key is used as the database look-up.

        Returns
        -------
        sympy.Basic or None
            The SymPy expression for σ(E) if the reaction is found, or
            ``None`` if no entry exists (e.g. for non-photoionisation
            reactions).

        Notes
        -----
        The expression uses the symbol ``E`` (photon energy in erg) as the
        independent variable and returns cross sections in cm².

        References
        ----------
        Verner, D. A. et al. 1996, ApJ, 465, 487
        """
        with JaffDb() as jdb:
            table = jdb.table("verner_cross_sections")
            rows: list = table.rows(conditions=f"reaction = '{reaction.serialized}'")

        if not rows:
            return None

        # Convert the stored string representation back to a SymPy expression.
        return sympify(rows[0]["xsecs"])

    def get_xsec(self, reaction: Reaction) -> XsecsProps | None:
        """
        Load the tabulated Leiden / NORAD cross sections for a reaction.

        The ``photo_reaction_cross_sections`` table maps the reaction key to an
        HDF5 group (Leiden preferred, NORAD as fallback) and records which
        processes are present.  The group is read into numpy arrays.

        Parameters
        ----------
        reaction : Reaction
            Reaction whose serialised key is used as the database look-up.

        Returns
        -------
        XsecsProps or None
            Dict with ``units`` (photon energy in eV, cross section in cm²),
            ``_equations`` (``pa`` photo-absorption flag and ``decay_type``,
            either ``"dissociation"`` or ``"ionization"``), the ``photon_energy``
            grid, the optional ``photo_absorption`` array, and the reaction's
            single ``photodecay`` cross-section array (``None`` when absent).
            Returns ``None`` if the reaction has no cross-section entry.
        """
        with JaffDb() as jdb:
            table = jdb.table("photo_reaction_cross_sections")
            rows: list = table.rows(conditions=f"reaction = '{reaction.serialized}'")

        if not rows:
            return None

        row = rows[0]
        loc: str = row["leiden"] if row["leiden"] else row["norad"]
        h5group = str(JAFF_DIR.resolve() / loc)
        pr_xsec = HDF5().to_dict(h5group)

        xsecs: XsecsProps = {
            "units": {
                "photon_energy": "eV",
                "cross_section": "cm^2",
            },
            "_equations": {
                "pa": bool(row["photo_absorption"]),
                "decay_type": row["decay_type"],
            },
            "photon_energy": pr_xsec.get("photon_energy", {}).get("_data", None),
            "photo_absorption": pr_xsec.get("photoabsorption", {}).get("_data", None),
            "photodecay": pr_xsec.get("photodecay", {}).get("_data", None),
        }

        return xsecs

    @staticmethod
    def shielding(reaction: Reaction, network: Network) -> Expr:
        """Build the symbolic shielding factor for a photo-reaction.

        The shielding function named by
        ``reaction._metadata["shielding"]["type"]`` is resolved from the
        shielding registry via
        :func:`~jaff.physics.photo_reactions.shielding._get_shielding_function`.
        Lookup is keyed by ``(type, reaction.serialized)`` and prefers a
        reaction-specific (local) function, falling back to a global one
        registered with ``reaction=None``.  The resolved
        :class:`~jaff.physics.photo_reactions.shielding._base.ShieldingFunction`
        instance's ``get_shielding(reaction, network)`` produces the factor.

        The result is cached on ``reaction._metadata["shielding"]["value"]`` so
        repeated calls (e.g. once per radiation band) reuse it.

        Parameters
        ----------
        reaction : Reaction
            Reaction to shield; its ``metadata["shielding"]["type"]`` selects
            the shielding function.
        network : Network
            Network the reaction belongs to, forwarded to the shielding
            function for species/column-density look-ups.

        Returns
        -------
        sympy.Expr
            Dimensionless shielding factor multiplying the photo-rate.

        Raises
        ------
        ParserError
            If no shielding function is registered for ``type`` either locally
            (for this reaction) or as a global fallback.
        """
        sprops = reaction._metadata["shielding"]

        shielding_fn = _get_shielding_function(sprops["type"], reaction.serialized)
        shielding_expr = shielding_fn.get_shielding(reaction, network)
        reaction._metadata["shielding"]["value"] = shielding_expr

        return shielding_expr
