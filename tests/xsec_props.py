# ABOUTME: Shared helper building XsecProps cross-section records for tests, with
# ABOUTME: neutral defaults so each test spells out only the fields it exercises

from typing import Literal, Optional

import numpy as np
from sympy import Expr

from jaff.core.reaction import PhotoReaction
from jaff.core.reaction.config import PhotoConfig, ReactionConfig, ReactionsConfig
from jaff.physics.photo_reactions._photochemistry import (
    XsecProcess,
    XsecProps,
    XsecUnits,
)


def make_xsec_props(
    *,
    photon_energy: Optional[np.ndarray] = None,
    photo_absorption: Optional[np.ndarray] = None,
    photodecay: Optional[np.ndarray] = None,
    photodecay_expr: Optional[Expr] = None,
    absorption: bool = False,
    decay_type: Literal["dissociation", "ionization"] = "ionization",
    database: Literal["norad", "verner", "leiden"] = "norad",
) -> XsecProps:
    """Build an :class:`XsecProps` in eV / cm² with the given data.

    Args:
        photon_energy: Tabulated photon energy grid (eV)
        photo_absorption: Tabulated photo-absorption cross section (cm²)
        photodecay: Tabulated photodecay cross section (cm²)
        photodecay_expr: Symbolic photodecay cross section in ``ph_E`` (eV)
        absorption: Weight the rate by the photo-absorption cross section
        decay_type: Decay channel of the reaction
        database: Source database label

    Returns:
        Cross-section record carrying the given fields
    """
    return XsecProps(
        units=XsecUnits(photon_energy="eV", cross_section="cm^2"),
        process=XsecProcess(absorption=absorption, decay_type=decay_type),
        database=database,
        photon_energy=photon_energy,
        photo_absorption=photo_absorption,
        photodecay=photodecay,
        photodecay_expr=photodecay_expr,
    )


def default_reaction_config() -> ReactionConfig:
    """Resolved default reaction settings (unshielded), for reaction stand-ins."""
    return ReactionsConfig().resolve("stand-in")


def stub_photo(
    xsec: XsecProps | None = None, config: PhotoConfig | None = None
) -> PhotoReaction:
    """A PhotoReaction with *xsec* and default (or given) photo settings, no reaction."""
    photo = config if config is not None else default_reaction_config().photo
    assert photo is not None

    return PhotoReaction(None, photo, xsec=xsec)  # type: ignore[arg-type]
