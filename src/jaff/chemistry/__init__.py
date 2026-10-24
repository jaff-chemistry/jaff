# ABOUTME: jaff.chemistry: chemistry reference data shared by species, elements and
# ABOUTME: networks: the periodic table of Element records, read once per process
"""Chemistry reference data shared by species, elements and networks.

Like :mod:`jaff.physics.constants`, this is reference data, not network state:
nothing here is built or configured per network.  The bundled element table is
read on first use and cached for the rest of the process, so a
:class:`~jaff.core.Specie` or :class:`~jaff.core.Element` works the same with or
without a :class:`~jaff.core.Network`.

* :func:`periodic_table` -- the whole table, or one
  :class:`~jaff.core.Element` by symbol.

Examples
--------
>>> from jaff.chemistry import periodic_table
>>> periodic_table("C").name
'Carbon'
>>> len(periodic_table()) > 100
True
"""

from ._periodic_table import periodic_table

__all__ = ["periodic_table"]
