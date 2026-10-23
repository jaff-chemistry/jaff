# ABOUTME: Per-reaction data split out of Reaction: per-band photo data and the
# ABOUTME: reaction's auxiliary functions

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from sympy import Expr

from ...core.parsers import AuxiliaryFunctionProps
from ...errors import ParserError

if TYPE_CHECKING:
    pass


class GroupedPhotoReactionProps:
    """One photo reaction's data for one radiation band.

    Built by :meth:`Radiation.get_reaction_rate_props` for every band, in
    ``net.radiation.groups`` order, and stored in ``reaction.photo.rgrouped``.
    A photo reaction whose rate comes from a ``.jfunc`` ``chemRate<N>`` is not
    coupled to the bands and gets none.

    Attributes
    ----------
    k : Expr
        Band loss coefficient: ``k`` times the reactant densities is the
        radiation the reaction removes from the band per unit time.  It is
        ``c · den[i] · <σ>_i`` (with any shielding), where ``<σ>_i`` is the
        absorption cross section weighted by photon number in ``"nph"`` mode
        and by photon energy in ``"u"`` mode (each absorbed photon takes its
        own energy).  The reaction rate itself always counts photons; see
        :meth:`Radiation.get_reaction_rate_props`.
    xsec : float or sympy.Expr
        Photon-number-weighted band-average cross section ``<σ>_i`` (cm²): of
        the photo-absorption table when the cross section has one
        (``xsec.process.absorption``), else of the photodecay cross section.
        An expression when a ``photoCrossSection`` uses other symbols (e.g.
        ``tgas``).
    xsec_frac : float or sympy.Expr
        ``xsec`` divided by the photodecay cross section averaged over the
        whole spectrum (``reaction.photo.xsec.total``); ``0.0`` when that
        average is zero.  Not a fraction of the rate: the values need not sum
        to 1 over the bands and a band where σ is above its spectrum average
        exceeds 1.

    Examples
    --------
    >>> from jaff import Network
    >>> net = Network(
    ...     "h_photoionization",
    ...     radiation={"bands": [1, 13.6, 100, "inf"]},
    ... )
    >>> props = net.reactions.photo_reactions()[0].photo.rgrouped[1]
    >>> props.k
    5.85154082319921e-8*photden[1]
    >>> float(props.xsec), float(props.xsec_frac)
    (1.9518639202054948e-18, 15.735906475383057)
    """

    def __init__(self, k: Expr, xsec: float | Expr | None, xsec_frac: float):
        """Store the band data; see the class attributes."""
        self.k: Expr = k
        self.xsec: float | Expr | None = xsec
        self.xsec_frac: float = xsec_frac


@dataclass
class AuxFunctions:
    chemrate: AuxiliaryFunctionProps | None = None
    deltae: AuxiliaryFunctionProps | None = None
    photocrosssection: AuxiliaryFunctionProps | None = None
    spectralyield: AuxiliaryFunctionProps | None = None

    def __post_init__(self) -> None:
        """Reject a reaction given both a custom rate and a custom cross section."""
        if self.chemrate is not None and self.photocrosssection is not None:
            raise ParserError(
                "A reaction cannot have both a chemRate and a photoCrossSection "
                "auxiliary function"
            )
