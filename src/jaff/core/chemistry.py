# ABOUTME: Chemistry reference data as module functions (not objects): the element
# ABOUTME: table used by species and elements, read once per process

"""Chemistry reference data shared by species, elements and networks.

Module-level functions, like :mod:`jaff.physics.constants`: nothing here is an
object to build or configure.  Each table is read from the bundled JAFF
database on first use and cached for the rest of the process, so a
:class:`~jaff.core.Specie` or :class:`~jaff.core.Element` works the same with or
without a :class:`~jaff.core.Network`.

Examples
--------
>>> from jaff.core import chemistry
>>> chemistry.element_table()["C"]["name"]
'Carbon'
"""

from __future__ import annotations

from functools import cache
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._typing import ElementProps

__all__ = ["element_table"]


@cache
def element_table() -> dict[str, ElementProps]:
    """The element data, keyed by element symbol; read once per process.

    From the bundled ``atomic_masses`` table (:class:`~jaff.drivers.sqlite.JaffDb`).
    Every caller gets the same dict: do not modify it.

    Returns
    -------
    dict[str, ElementProps]
        Mapping from element symbol (e.g. ``"H"``) to an
        :class:`~jaff.core.ElementProps` dict with keys ``name``, ``mass`` (g),
        ``atomic_mass``, ``protons``, ``neutrons`` and ``electrons``.

    Examples
    --------
    >>> from jaff.core import chemistry
    >>> chemistry.element_table()["H"]["mass"]
    1.673773e-24
    """
    from ..drivers.sqlite import JaffDb

    with JaffDb() as jdb:
        rows = jdb.table("atomic_masses").all_rows()

    return {
        row["element"]: {
            "name": row["name"],
            "mass": row["mass"],
            "atomic_mass": row["atomic_mass"],
            "protons": row["protons"],
            "neutrons": row["neutrons"],
            "electrons": row["electrons"],
        }
        for row in rows
    }
