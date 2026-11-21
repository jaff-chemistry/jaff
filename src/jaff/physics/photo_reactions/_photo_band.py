# ABOUTME: PhotoBand: one photo reaction's data for one radiation band, holding
# ABOUTME: the RadiationGroup it belongs to (frozen; held in reaction.photo.bands)

"""The :class:`PhotoBand` per-band data of a photo reaction."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from sympy import Expr

if TYPE_CHECKING:
    from ._radiation import RadiationGroup


@dataclass(frozen=True)
class PhotoBand:
    """One photo reaction's data for one radiation band.

    Built for every band, in ``net.radiation.groups`` order, and held in
    ``reaction.photo.bands``.  A photo reaction whose rate comes from a
    ``.jfunc`` ``chemRateCoeff<N>`` is not coupled to the bands and gets none.
    Immutable: :meth:`~jaff.core.reaction.PhotoReaction.substitute` replaces
    the bands.

    Attributes
    ----------
    group : RadiationGroup
        The radiation band (``net.radiation.groups`` entry) the data is for.
    k : Expr
        Band loss coefficient: ``k`` times the reactant densities is the
        radiation the reaction removes from the band per unit time.  It is
        ``c · den[i] · <σ>_i`` (with any shielding), where ``<σ>_i`` is the
        absorption cross section weighted by photon number in ``"nph"`` mode
        and by photon energy in ``"u"`` mode (each absorbed photon takes its
        own energy).  The reaction rate itself always counts photons; see
        :meth:`~jaff.core.reaction.PhotoReaction.band_rates`.
    xsec : float or sympy.Expr
        Photon-number-weighted band-average cross section ``<σ>_i`` (cm²): of
        the photo-absorption table when the cross section has one
        (``xsec.process.absorption``), else of the photodecay cross section.
        An expression when a ``photoCrossSection`` uses other symbols (e.g.
        ``tgas``).
    xsec_frac : float or sympy.Expr
        ``xsec`` divided by the same table's cross section averaged over the
        whole spectrum (photo-absorption when the cross section has one, else
        photodecay, i.e. ``reaction.photo.xsec_total``); ``0.0`` when that
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
    >>> band = net.reactions.photo_reactions()[0].photo.bands[1]
    >>> band.group.index, band.k
    (1, 5.85154082319921e-8*photden[1])
    >>> float(band.xsec), float(band.xsec_frac)
    (1.9518639202054948e-18, 15.735906475383057)
    """

    group: RadiationGroup
    k: Expr
    xsec: float | Expr
    xsec_frac: float | Expr
