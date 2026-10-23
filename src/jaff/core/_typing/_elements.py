# ABOUTME: ElementProps: typed dict for one entry of the JAFF mass dictionary
# ABOUTME: Keys become the attributes of the same name on Element

from typing import TypedDict


class ElementProps(TypedDict):
    """Properties of one entry of the JAFF mass dictionary.

    The mass dictionary maps an element symbol (e.g. ``"H"``) to an
    ``ElementProps`` dict.  It is loaded from the bundled ``atomic_masses``
    table by default and can be replaced with ``Elements.configure``
    (or ``Element.configure``).  Each key becomes the attribute of the
    same name on :class:`Element`.

    Attributes
    ----------
    name : str
        Full element name (e.g. ``"Hydrogen"``).
    mass : float
        Mass of the element in grams (CGS).
    atomic_mass : float
        Standard atomic weight in atomic mass units.
    protons : int
        Number of protons (atomic number).
    neutrons : int
        Number of neutrons in the most common isotope.
    electrons : int
        Number of electrons in the neutral atom.

    Notes
    -----
    In the bundled table, pseudo-elements such as ``"e-"`` and ``"GRAIN"``
    have ``None`` for ``atomic_mass``, ``protons``, ``neutrons`` and
    ``electrons``.

    Examples
    --------
    >>> from jaff.core import ElementProps
    >>> hydrogen: ElementProps = {
    ...     "name": "Hydrogen",
    ...     "mass": 1.673773e-24,
    ...     "atomic_mass": 1.008,
    ...     "protons": 1,
    ...     "neutrons": 0,
    ...     "electrons": 1,
    ... }
    """

    name: str
    mass: float
    atomic_mass: float
    protons: int
    neutrons: int
    electrons: int
