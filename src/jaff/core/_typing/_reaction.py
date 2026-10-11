from typing import TYPE_CHECKING, NotRequired, TypedDict

from sympy import Basic, Expr

if TYPE_CHECKING:
    from ...physics.photo_reactions._photochemistry import XsecsProps
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
        "dE": Basic,
        "dRad": Basic,
        "custom_rad_rate": bool,
        "reaction_type": str,
        "tmin": float | None,
        "tmax": float | None,
        "t_cutoff": NotRequired[str],
        "rate_segments": NotRequired[list["RateSegmentProps"]],
        "original_string": str,
        "pi_database": NotRequired[str | None],
        "xsecs_dict": "XsecsProps",
    },
)
