"""
Look-ups for photochemical cross sections.

Three databases are used:

- **NORAD** and **Leiden** -- tabulated cross sections, indexed by the
  ``photo_reaction_cross_sections`` table and read from the corresponding HDF5
  group as numpy arrays.
- **Verner** (1996) -- analytic photoionisation fits stored as SymPy strings
  in the ``verner_cross_sections`` table, returned as an expression in ``ph_E``.

Photodissociation always comes from Leiden.  Photoionisation follows the
reaction's ``pi_database`` override, then the global
``network.reactions.photo.pi_database``, then the remaining databases (see
:meth:`Photochemistry.get_xsec`).  All look-ups are keyed by
``reaction.serialized`` (or the normalised proxy string).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
from sympy import Expr, sympify

from ...config import JAFF_DIR
from ...drivers import HDF5, JaffDb
from ...drivers.pooch import (
    download_shielding,
    download_xsecs,
)
from ...errors import ParserError
from ...io import JaffLogger
from ._radiation_props import PI_DATABASES
from .shielding import _get_shielding_function

if TYPE_CHECKING:
    from ...core import Reaction
    from ...core.network import Network


@dataclass
class XsecUnits:
    photon_energy: str
    cross_section: str


@dataclass
class XsecProcess:
    absorption: bool
    decay_type: Literal["dissociation", "ionization"]


@dataclass(eq=False)
class XsecProps:
    units: XsecUnits
    process: XsecProcess
    database: Literal["norad", "verner", "leiden", "custom"]
    photon_energy: np.ndarray | None
    photo_absorption: np.ndarray | None
    photodecay: np.ndarray | None
    photodecay_expr: Expr | None
    # Cross section averaged over the whole radiation spectrum (cm²); 0.0 until
    # the rate is computed when radiation is enabled, None without radiation.
    # An expression when a photoCrossSection uses other symbols (e.g. tgas).
    total: float | Expr | None = None


class Photochemistry:
    """Cross-section provider for photo-reactions.

    Resolves per-reaction photoionisation / photodissociation / photoabsorption
    cross sections from the JAFF database and the tabulated Leiden, NORAD, and
    Verner data files.
    """

    def __init__(self, network: Network):
        """Ensure the cross-section and shielding data files are available locally.

        Constructing a :class:`Photochemistry` triggers
        :func:`~jaff.drivers.pooch.download_xsecs` (Leiden / NORAD / Verner
        cross sections) and :func:`~jaff.drivers.pooch.download_shielding` (the
        Leiden line-shielding tables), downloading both on first use (a network
        fetch unless already cached). Instantiate once and reuse rather than per
        reaction.
        """
        self.net = network
        self.logger: logging.Logger = JaffLogger().get_logger()
        # (key, requested, used) for photoionizations resolved from a fallback
        # database; reported once by log_fallback_summary().
        self.fallbacks: list[tuple[str, str, str]] = []
        download_xsecs()
        download_shielding()

    def _lookup_key(self, reaction: Reaction) -> str:
        """Database key: the proxy string when the reaction uses proxy lookups."""
        photo = self.net.config.reactions.resolve(reaction.serialized).photo
        if photo is not None and photo.use_proxy_photoreaction:
            return reaction.normalized_proxy_reaction_str()

        return reaction.serialized

    def get_verner_xsec(self, reaction: Reaction) -> Expr | None:
        """
        Query the JAFF database for the Verner photoionisation cross section.

        Verner cross sections are analytic fits to photoionisation cross
        sections from Verner et al. (1996) stored as SymPy-parseable strings
        in the ``verner_cross_sections`` SQLite table.

        Parameters
        ----------
        reaction : Reaction
            Reaction whose look-up key (see :meth:`_lookup_key`) is queried.

        Returns
        -------
        sympy.Expr or None
            The SymPy expression for σ(E) if the reaction is found, or
            ``None`` if no entry exists.

        Notes
        -----
        The expression uses the symbol ``ph_E`` (photon energy in eV) as the
        independent variable and returns cross sections in cm²; it is zero
        outside ``[E_th, E_max]``.

        References
        ----------
        Verner, D. A. et al. 1996, ApJ, 465, 487
        """
        return self._verner_expr(self._lookup_key(reaction))

    @staticmethod
    def _verner_expr(key: str) -> Expr | None:
        with JaffDb() as jdb:
            table = jdb.table("verner_cross_sections")
            rows: list = table.rows(conditions=f"reaction = '{key}'")

        return sympify(rows[0]["xsecs"]) if rows else None

    @staticmethod
    def _tabulated_row(key: str) -> dict[str, Any] | None:
        with JaffDb() as jdb:
            table = jdb.table("photo_reaction_cross_sections")
            rows: list = table.rows(conditions=f"reaction = '{key}'")

        return rows[0] if rows else None

    @staticmethod
    def _is_ionization(key: str, row: dict[str, Any] | None) -> bool:
        """Decay type from the table row, else from an ``e-`` product."""
        if row is not None:
            return row["decay_type"] == "ionization"

        return "e-" in key.partition("__")[2].split(".")

    @staticmethod
    def _pi_override(reaction: Reaction) -> str | None:
        """The reaction's own ``pi_database``, or ``None`` when it has none."""

        return reaction.photo.pi_database if reaction.photo is not None else None

    def _pi_order(self, reaction: Reaction) -> list[str]:
        """Databases to try: the chosen one, then :data:`PI_DATABASES` order."""
        chosen = (
            self._pi_override(reaction)
            or self.net.config.reactions.photo.pi_database
            or "norad"
        )

        return list(dict.fromkeys([chosen, *PI_DATABASES]))

    def _default_total(self) -> float | None:
        """Initial :attr:`XsecProps.total`: ``0.0`` with radiation, else ``None``."""
        return 0.0 if self.net.radiation is not None else None

    @staticmethod
    def _tabulated_xsec(
        db: str, loc: str, row: dict[str, Any], total: float | None
    ) -> XsecProps:
        """Read one tabulated (Leiden / NORAD) HDF5 group into an XsecProps."""
        pr_xsec = HDF5().to_dict(str(JAFF_DIR.resolve() / loc))
        photo_absorption = pr_xsec.get("photoabsorption", {}).get("_data", None)

        return XsecProps(
            units=XsecUnits(photon_energy="eV", cross_section="cm^2"),
            process=XsecProcess(
                absorption=bool(row["photo_absorption"]) and photo_absorption is not None,
                decay_type=row["decay_type"],
            ),
            database=db,  # type: ignore
            photon_energy=pr_xsec.get("photon_energy", {}).get("_data", None),
            photo_absorption=photo_absorption,
            photodecay=pr_xsec.get("photodecay", {}).get("_data", None),
            photodecay_expr=None,
            total=total,
        )

    def _load(self, db: str, key: str, row: dict[str, Any] | None) -> XsecProps | None:
        """Cross sections for *key* from database *db*, or ``None`` if absent."""
        if db == "verner":
            expr = self._verner_expr(key)

            if expr is None:
                return None

            return XsecProps(
                units=XsecUnits(photon_energy="eV", cross_section="cm^2"),
                process=XsecProcess(absorption=False, decay_type="ionization"),
                database="verner",
                photon_energy=None,
                photo_absorption=None,
                photodecay=None,
                photodecay_expr=expr,
                total=self._default_total(),
            )

        if row is None or not row[db]:
            return None

        return self._tabulated_xsec(db, row[db], row, self._default_total())

    def get_custom_xsec(self, reaction: Reaction, expr: Expr) -> XsecProps:
        """Cross sections from a ``.jfunc`` ``photoCrossSection`` function.

        Parameters
        ----------
        reaction : Reaction
            Reaction the cross section belongs to; its decay type is resolved
            as in :meth:`get_xsec`.
        expr : Expr
            Photodecay cross section (cm²) as an expression in ``ph_E`` (eV).

        Returns
        -------
        XsecProps
            Symbolic cross sections with ``database="custom"`` and no
            photo-absorption.
        """
        key = self._lookup_key(reaction)
        ionization = self._is_ionization(key, self._tabulated_row(key))

        return XsecProps(
            units=XsecUnits(photon_energy="eV", cross_section="cm^2"),
            process=XsecProcess(
                absorption=False,
                decay_type="ionization" if ionization else "dissociation",
            ),
            database="custom",
            photon_energy=None,
            photo_absorption=None,
            photodecay=None,
            photodecay_expr=expr,
            total=self._default_total(),
        )

    def get_xsec(self, reaction: Reaction, required: bool = False) -> XsecProps | None:
        """
        Resolve the cross sections for a photo-reaction.

        Photodissociation comes from Leiden.  Photoionisation tries, in order,
        ``reaction.photo.pi_database`` (or the global
        ``network.reactions.photo.pi_database``, ``"norad"`` by default),
        then ``norad``, ``verner``, ``leiden`` (duplicates removed).  Fallbacks to
        a later database are recorded and reported once by
        :meth:`log_fallback_summary`.

        The look-up key is ``reaction.normalized_proxy_reaction_str()`` when
        the reaction's ``use_proxy_photoreaction`` setting is on, otherwise
        ``reaction.serialized``.

        Parameters
        ----------
        reaction : Reaction
            Reaction to resolve.
        required : bool, optional
            Raise when no database has a photoionisation cross section
            (default ``False``: return ``None``).  The network sets this when
            radiation is enabled and the reaction has no custom rate.

        Returns
        -------
        XsecProps or None
            Object with ``units``, ``process`` (``absorption`` flag and
            ``decay_type``), the ``database`` used, and either tabulated arrays
            (``photon_energy``, optional ``photo_absorption``, ``photodecay``)
            or a SymPy ``photodecay_expr`` in ``ph_E`` (Verner).  ``None`` if no
            entry exists and *required* is false.

        Raises
        ------
        ParserError
            If ``reaction.photo.pi_database`` is set on a non-photoionisation
            reaction, or *required* and no database has the reaction.
        """
        key = self._lookup_key(reaction)
        row = self._tabulated_row(key)

        if not self._is_ionization(key, row):
            if self._pi_override(reaction) is not None:
                raise ParserError(
                    f"pi_database is only valid for photoionization reactions: {key}"
                )

            return self._load("leiden", key, row)

        order = self._pi_order(reaction)
        for db in order:
            xsecs = self._load(db, key, row)
            if xsecs is None:
                continue

            if db != order[0]:
                self.fallbacks.append((key, order[0], db))

            return xsecs

        if required:
            raise ParserError(
                f"No photoionization cross section for {key} in any database "
                f"(tried {', '.join(order)})"
            )

        return None

    def log_fallback_summary(self) -> None:
        """Log one warning listing every photoionization database fallback.

        Groups the recorded fallbacks as ``requested -> used: key, key`` and
        clears them.  Does nothing when no fallback occurred.
        """
        if not self.fallbacks:
            return

        groups: dict[tuple[str, str], list[str]] = {}
        for key, requested, used in self.fallbacks:
            groups.setdefault((requested, used), []).append(key)

        lines = [
            f"{req} -> {used}: {', '.join(keys)}" for (req, used), keys in groups.items()
        ]
        self.logger.warning(
            "Photoionization cross sections taken from a fallback database:\n"
            + "\n".join(lines)
        )

        self.fallbacks.clear()

    @staticmethod
    def shielding(reaction: Reaction, network: Network) -> Expr:
        """Build the symbolic shielding factor for a photo-reaction.

        The shielding function named by ``reaction.photo.shielding.type`` is
        resolved from the
        shielding registry via
        :func:`~jaff.physics.photo_reactions.shielding._get_shielding_function`.
        Lookup is keyed by ``(type, reaction.serialized)`` and prefers a
        reaction-specific (local) function, falling back to a global one
        registered with ``reaction=None``.  The resolved
        :class:`~jaff.physics.photo_reactions.shielding._base.ShieldingFunction`
        instance's ``get_shielding(reaction, network)`` produces the factor.

        The result is not cached: ``photo.shielding.props`` holds only the
        user's settings, so callers evaluate it once per reaction.

        Parameters
        ----------
        reaction : Reaction
            Reaction to shield; its ``photo.shielding.type`` selects
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
        shielding = reaction.photo.shielding if reaction.photo is not None else None
        assert shielding is not None, f"{reaction} is not shielded"

        shielding_fn = _get_shielding_function(shielding.type, reaction.serialized)

        return shielding_fn.get_shielding(reaction, network)
