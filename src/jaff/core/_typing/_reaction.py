from typing import TYPE_CHECKING, Any, NotRequired, TypedDict

from sympy import Expr

if TYPE_CHECKING:
    from ...physics import XsecProps
    from .. import Network, Reaction, Specie, Species
else:
    Specie = "Specie"
    Species = "Species"
    Reaction = "Reaction"
    Network = "Network"


RateSegmentProps = TypedDict(
    "RateSegmentProps",
    {
        "rate": Expr,
        "tmin": float | None,
        "tmax": float | None,
    },
)


ReactionProps = TypedDict(
    "ReactionProps",
    {
        "reactants": list["Specie"],
        "products": list["Specie"],
        "rate": Expr,
        "dE": Expr,
        "custom_rad_rate": bool,
        "reaction_type": str,
        "format": NotRequired[str | None],
        "tmin": float | None,
        "tmax": float | None,
        "t_cutoff": NotRequired[str],
        "rate_segments": NotRequired[list["RateSegmentProps"]],
        "original_string": str,
        "index": NotRequired[int],
        "pi_database": NotRequired[str | None],
        "shielding": NotRequired[dict[str, Any] | None],
        "xsec_props": NotRequired["XsecProps | None"],
        "spectral_yield": NotRequired[Expr | None],
    },
)
"""Typed dictionary of one reaction as stored in a ``.jaff`` file.

:func:`~jaff.io._io.from_jaff_file` decodes each stored reaction into one
of these (the ``reactions`` list of ``JaffProps``), and
:class:`~jaff.core.Network` rebuilds a :class:`~jaff.core.Reaction` from
it when loading a ``.jaff`` file.

Attributes
----------
reactants : list[Specie]
    Reactant species.
products : list[Specie]
    Product species.
rate : Expr
    Rate coefficient (``Reaction.rate_coefficient``).
dE : Expr
    Energy released per reaction event (erg).
custom_rad_rate : bool
    ``True`` when a photo reaction's rate came from a ``.jfunc``
    ``chemRate<N>`` (``Reaction.photo.custom_rad_rate``).
reaction_type : str
    Reaction type (``Reaction.type``), e.g. ``"photo"``.
format : str or None
    Optional key.  Name of the parser that read the reaction; absent or
    ``None`` in files written before it was stored.
tmin : float or None
    Lower temperature bound in Kelvin, or ``None``.
tmax : float or None
    Upper temperature bound in Kelvin, or ``None``.
t_cutoff : str
    Optional key.  Out-of-range behaviour of the rate segments (``"clip"``
    when absent).
rate_segments : list[RateSegmentProps]
    Optional key.  Per-temperature-range pieces, each a
    ``{"rate", "tmin", "tmax"}`` dict.
original_string : str
    The raw network-file line of the reaction.
index : int
    Optional key.  File-side source number (``Reaction.index``); the list
    position when absent.
pi_database : str or None
    Optional key.  Per-reaction photoionization database override of a
    photo reaction.
shielding : dict or None
    Optional key.  Shielding of a photo reaction as
    ``{"type": ..., "props": {...}}`` (the arguments of
    :class:`~jaff.core.reaction.ShieldingProps`).
xsec_props : XsecProps or None
    Optional key.  Cross sections of a photo reaction
    (``Reaction.photo.xsec``).
spectral_yield : Expr or None
    Optional key.  Emitted radiation per event
    (``Reaction.rad.spectral_yield``).
"""
