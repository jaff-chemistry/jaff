from typing import Literal, TypedDict

import numpy as np
from sympy import Basic

units = TypedDict("units", {"photon_energy": str, "cross_section": str})
_equations = TypedDict(
    "_equations",
    {
        "pa": bool,
        "decay_type": Literal["dissociation", "ionization"],
    },
)

XsecsProps = TypedDict(
    "XsecsProps",
    {
        "units": units,
        "_equations": _equations,
        "database": Literal["norad", "verner", "leiden"],
        "photon_energy": np.ndarray | None,
        "photo_absorption": np.ndarray | None,
        "photodecay": np.ndarray | None,
        "photodecay_expr": Basic | None,
    },
)
