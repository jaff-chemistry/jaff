# ABOUTME: Per-reaction property groups split out of Reaction: photo-reaction state
# ABOUTME: (cross sections, bands) and radiation coupling (energy added to the field)

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Dict, List

from sympy import Expr

from ...core.parsers import AuxiliaryFunctionProps
from ...errors import ParserError

if TYPE_CHECKING:
    from ...physics import XsecProps


class PhotoReactionProps:
    """Photo-reaction state, held as ``Reaction.photo`` (``None`` otherwise).

    Attributes
    ----------
    xsec : XsecProps or None
        Cross sections of the reaction's single decay channel, or ``None``
        when no database has an entry.  ``xsec.total`` holds the cross
        section averaged over the whole radiation spectrum (cm²).
    pi_database : str or None
        Per-reaction photoionization database override (``"norad"``,
        ``"verner"`` or ``"leiden"``); ``None`` uses the global
        :attr:`RadiationProps.pi_database`.
    override_custom_rate : bool or None
        Per-reaction override of :attr:`RadiationProps.override_custom_rate`;
        ``None`` uses the global flag.
    override_custom_cross_section : bool or None
        Per-reaction override of
        :attr:`RadiationProps.override_custom_cross_section`; ``None`` uses the
        global flag (``False`` without radiation).
    custom_rad_rate : bool
        ``True`` when the rate comes from a ``.jfunc`` aux function rather
        than from the cross sections.
    rgrouped : list[GroupedPhotoReactionProps] or None
        One entry per radiation band, in ``net.radiation.groups`` order, set
        when a radiation field is configured; ``None`` otherwise.  See
        :attr:`Reaction.band_xsecs`.
    shielding : ShieldingProps or None
        Shielding applied to the photo rate, or ``None`` when the reaction is
        not shielded.
    """

    def __init__(
        self,
        xsec: XsecProps | None = None,
        pi_database: str | None = None,
        override_custom_rate: bool | None = None,
        override_custom_cross_section: bool | None = None,
        custom_rad_rate: bool = False,
        rgrouped: List[GroupedPhotoReactionProps] | None = None,
        shielding: ShieldingProps | None = None,
    ):
        """Store the photo-reaction state; see the class attributes."""
        self.xsec: XsecProps | None = xsec
        self.pi_database: str | None = pi_database
        self.override_custom_rate: bool | None = override_custom_rate
        self.override_custom_cross_section: bool | None = override_custom_cross_section
        self.custom_rad_rate: bool = custom_rad_rate
        self.rgrouped: List[GroupedPhotoReactionProps] | None = rgrouped
        self.shielding: ShieldingProps | None = shielding


class GroupedPhotoReactionProps:
    """One photo reaction's data for one radiation band.

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
    xsec : float, sympy.Expr or None
        Photon-number-weighted band-average cross section (cm²): an expression
        when a ``photoCrossSection`` uses other symbols (e.g. ``tgas``), or
        ``None`` for a custom-rate reaction.
    xsec_frac : float
        Fraction of the reaction's rate attributed to this band.
    """

    def __init__(self, k: Expr, xsec: float | Expr | None, xsec_frac: float):
        """Store the band data; see the class attributes."""
        self.k: Expr = k
        self.xsec: float | Expr | None = xsec
        self.xsec_frac: float = xsec_frac


class ReactionRadiationProps:
    """Radiation emitted by any reaction, held as ``Reaction.rad``.

    Attributes
    ----------
    spectral_yield : Expr or None
        Radiation energy the reaction adds to the field per event per unit
        photon energy (erg/eV), integrated over solid angle, as an expression in
        ``ph_E`` (eV); from a ``.jfunc`` ``spectralYield`` function.  ``None``
        when the reaction emits nothing.
    group_yields : list[Expr] or None
        ``spectral_yield`` integrated over each radiation group, in
        ``net.radiation.groups`` order, in the units of the group density per
        event (erg in ``"u"`` mode, photons in ``"nph"`` mode); see
        :meth:`Radiation.integrate_over_groups`.  ``None`` without radiation or
        without a yield.
    """

    def __init__(
        self,
        spectral_yield: Expr | None = None,
        group_yields: List[Expr] | None = None,
    ):
        """Store the emission; see the class attributes."""
        self.spectral_yield: Expr | None = spectral_yield
        self.group_yields: List[Expr] | None = group_yields


class ShieldingProps:
    """Shielding of a photo rate, held as ``Reaction.photo.shielding``.

    Attributes
    ----------
    type : str
        Shielding model (e.g. ``"leiden"``, ``"db1996"``).
    props : dict[str, Any]
        Model-specific settings (e.g. ``shielded_by``, ``min_ncol``) as given
        by the user.
    tables : dict[str, Any]
        Data tables the model built for the generated code to interpolate
        over, keyed by file name (e.g. ``shielding_<reaction>.hdf5``); written
        out by :meth:`~jaff.core.network.Network.write_shielding_tables`.
    """

    def __init__(self, type: str, props: Dict[str, Any] | None = None):
        """Store the model; ``props=None`` (default) starts empty."""
        self.type: str = type
        self.props: Dict[str, Any] = props if props is not None else {}
        self.tables: Dict[str, Any] = {}


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
