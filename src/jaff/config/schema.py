"""Pydantic models describing the ``jaff.toml`` / ``jaffgen.toml`` config schema.

The models mirror the class structure they configure: ``[network.reactions]``
holds the settings applied to every reaction, and
``[network.reactions.reaction."<serialized>"]`` overrides them for one
reaction.  Every model rejects unknown keys, is strict about types (a TOML
string is never coerced into a number or a bool) and is frozen once validated.

String choices (``t_cutoff``, ``pi_database``, ...) are case-insensitive and
stored lower-cased.

A ``shielding`` table is checked against the shielding plugin its ``type``
names; a per-reaction table may name a plugin bound to that reaction.
"""

from __future__ import annotations

from typing import Annotated, Any, ClassVar, Dict, List, Literal, Self, TypeVar

import astropy.constants as _ac
from pydantic import (
    BaseModel,
    BeforeValidator,
    ConfigDict,
    ValidationError,
    field_validator,
    model_validator,
)

from ..errors import ParserError


def _lower(value: Any) -> Any:
    """Lower-case *value* when it is a string, so choices are case-insensitive."""
    return value.lower() if isinstance(value, str) else value


LowerStr = Annotated[str, BeforeValidator(_lower)]
TCutoff = Annotated[Literal["clip", "extrapolate"], BeforeValidator(_lower)]
DuplicatePolicy = Annotated[
    Literal["preserve-first", "preserve-last", "error"], BeforeValidator(_lower)
]
PiDatabase = Annotated[Literal["norad", "verner", "leiden"], BeforeValidator(_lower)]
RadiationMode = Annotated[Literal["nph", "u"], BeforeValidator(_lower)]
BackgroundField = Annotated[
    Literal[
        "bb_4000",
        "bb_10000",
        "bb_20000",
        "draine",
        "habing",
        "mathis",
        "solar",
        "tw_hydra",
    ],
    BeforeValidator(_lower),
]

_ModelT = TypeVar("_ModelT", bound=BaseModel)


class _Config(BaseModel):
    """Base of every config model: unknown keys are errors, types are strict and
    instances are frozen."""

    model_config = ConfigDict(extra="forbid", strict=True, frozen=True)


def _overlay(base: _ModelT, over: BaseModel | None) -> _ModelT:
    """Return *base* with every key set (not ``None``) in *over* applied.

    Nested models are overlaid field by field rather than replaced, so an
    override that sets one nested key keeps the others from *base* -- except
    an atomic model (``_atomic = True``), which replaces the base one whole.
    """
    if over is None:
        return base

    update: dict[str, Any] = {}
    for name in type(over).model_fields:
        value = getattr(over, name)
        if value is None:
            continue

        current = getattr(base, name)
        atomic = getattr(type(value), "_atomic", False)
        if isinstance(value, BaseModel) and current is not None and not atomic:
            value = _overlay(current, value)

        update[name] = value

    return base.model_copy(update=update)


class ShieldingParams(_Config):
    """Base of the settings of one shielding type (every key but ``type``).

    Each shielding plugin exposes its subclass through
    :meth:`~jaff.physics.photo_reactions.shielding.ShieldingFunction.config_model`;
    this base has no keys, so a plugin without settings accepts none.
    """


class ShieldingConfig(_Config):
    """Shielding of a photo reaction (a ``shielding`` table).

    Written flat in TOML (``{ type = "leiden", shielded_by = ["self"] }``):
    ``type`` names the shielding plugin and every other key is one of its
    settings, checked by the plugin's ``config_model()`` once the reaction it
    applies to is known (see :class:`ReactionsConfig`).

    The table is atomic: a per-reaction one replaces the global one whole.

    Attributes
    ----------
    type : str
        Shielding plugin name, lower-cased.
    params : dict
        The plugin's settings, as written.
    """

    _atomic: ClassVar[bool] = True

    type: LowerStr = "leiden"
    params: Dict[str, Any] = {}

    @model_validator(mode="before")
    @classmethod
    def _split_type(cls, data: Any) -> Any:
        """Split a flat table into ``type`` and the plugin ``params``."""
        if not isinstance(data, dict):
            return data

        params = {k: v for k, v in data.items() if k != "type"}

        return {"type": data.get("type", "leiden"), "params": params}


def _check_shielding(shielding: ShieldingConfig | None, reaction: str | None) -> None:
    """Check *shielding* against its plugin, for *reaction* (``None``: global).

    Raises
    ------
    ValueError
        If no plugin of that type applies to *reaction*, or the plugin rejects
        a setting.
    """
    if shielding is None:
        return

    # Imported here: the plugins import this module.
    from ..physics.photo_reactions.shielding import _get_shielding_function

    where = (
        "reactions.photo.shielding"
        if reaction is None
        else f'reactions.reaction."{reaction}".photo.shielding'
    )
    try:
        plugin = _get_shielding_function(shielding.type, reaction)
        plugin.config_model().model_validate(shielding.params)
    except (ParserError, ValidationError) as e:
        raise ValueError(f"{where}: {e}") from e


class PhotoConfig(_Config):
    """Photo-reaction settings (``[network.reactions.photo]``).

    Attributes
    ----------
    pi_database : str or None
        Photoionisation cross-section database: ``"norad"``, ``"verner"`` or
        ``"leiden"``.
    override_custom_rate : bool or None
        Recompute a ``.jfunc`` ``chemRate<N>`` rate from the cross sections.
    override_custom_cross_section : bool or None
        Ignore ``.jfunc`` ``photoCrossSection<N>`` and use the database.
    use_proxy_photoreaction : bool or None
        Look cross sections up under the proxy reaction (``CHx`` → ``CH``).
    shielding : ShieldingConfig or None
        Line shielding of the rate; ``None`` globally means no shielding.

    ``None`` means unset: the value comes from the global settings.
    """

    pi_database: PiDatabase | None = None
    override_custom_rate: bool | None = None
    override_custom_cross_section: bool | None = None
    use_proxy_photoreaction: bool | None = None
    shielding: ShieldingConfig | None = None


_DEFAULT_PHOTO = PhotoConfig(
    pi_database="norad",
    override_custom_rate=False,
    override_custom_cross_section=False,
    use_proxy_photoreaction=True,
)


class ReactionConfig(_Config):
    """Settings of one reaction (``[network.reactions.reaction."<serialized>"]``).

    Attributes
    ----------
    duplicate_policy : str or None
        How two rate coefficients over the same temperature range are
        resolved: ``"preserve-first"``, ``"preserve-last"`` or ``"error"``.
    t_cutoff : str or None
        Rate behaviour outside the temperature range: ``"clip"`` or
        ``"extrapolate"``.
    photo : PhotoConfig or None
        Photo-reaction settings.

    ``None`` means unset: the value comes from the global settings.
    """

    duplicate_policy: DuplicatePolicy | None = None
    t_cutoff: TCutoff | None = None
    photo: PhotoConfig | None = None


class ReactionsConfig(ReactionConfig):
    """Settings of every reaction (``[network.reactions]``).

    Holds the same keys as :class:`ReactionConfig`, with concrete defaults,
    plus the per-reaction entries.

    Attributes
    ----------
    reaction : dict[str, ReactionConfig]
        Per-reaction settings keyed by serialized reaction; a key set here
        wins over the global one (see :meth:`resolve`).
    """

    duplicate_policy: DuplicatePolicy = "preserve-first"
    t_cutoff: TCutoff = "clip"
    photo: PhotoConfig = _DEFAULT_PHOTO
    reaction: Dict[str, ReactionConfig] = {}

    @field_validator("photo", mode="after")
    @classmethod
    def _fill_photo_defaults(cls, photo: PhotoConfig) -> PhotoConfig:
        """Complete a partially given global ``photo`` table with the defaults."""
        return _overlay(_DEFAULT_PHOTO, photo)

    @model_validator(mode="after")
    def _check_shielding_tables(self) -> Self:
        """Check every shielding table against its plugin and reaction."""
        _check_shielding(self.photo.shielding, None)
        for serialized, entry in self.reaction.items():
            if entry.photo is not None:
                _check_shielding(entry.photo.shielding, serialized)

        return self

    def resolve(self, serialized: str) -> ReactionConfig:
        """Return the settings of one reaction, every key set.

        The global settings, with the keys set in the reaction's own entry (if
        any) applied over them field by field.

        Parameters
        ----------
        serialized : str
            Serialized reaction (e.g. ``"H+.e-__H"``).

        Returns
        -------
        ReactionConfig
            The reaction's resolved settings.
        """
        base = ReactionConfig(
            **{name: getattr(self, name) for name in ReactionConfig.model_fields}
        )

        return _overlay(base, self.reaction.get(serialized))


class RadiationConfig(_Config):
    """Radiation field (``[network.radiation]``); the table enables radiation.

    Only the shape of the settings is checked here; whether the band-averaged
    spectrum converges is checked by the radiation field itself.

    Attributes
    ----------
    bands : list of (float or "inf")
        Ordered photon-energy band edges in eV, at least two; ``"inf"`` is
        allowed as the last edge only.
    profile_index : float or list of float
        Spectral index *α* of the photon-number spectrum ``n(E) ∝ E^(α-2)``:
        one for every band, or a list with one per band.
    mode : str
        Tracked quantity: ``"nph"`` (photon number density) or ``"u"``
        (energy density).
    c : float or str
        Speed of light in cm/s, or a symbol name (e.g. ``"c_hat"``) for a
        reduced speed of light.
    background_field : str
        Reference radiation field used to scale ``chi_pe``.
    absorption_attenuates_flux : bool
        Whether the ``globalAbsorptionCoefficient`` also attenuates the flux
        equations, not only the density equations.
    """

    bands: List[float | Literal["inf"]]
    profile_index: float | List[float] = 0.0
    mode: RadiationMode = "nph"
    c: float | str = _ac.c.cgs.value  # ty: ignore[unresolved-attribute]
    background_field: BackgroundField = "draine"
    absorption_attenuates_flux: bool = True

    @field_validator("bands", mode="after")
    @classmethod
    def _check_bands(cls, bands: List[float | str]) -> List[float | str]:
        """At least two edges, with ``"inf"`` only as the last one."""
        if len(bands) < 2:
            raise ValueError(f"need at least two edges (one band), got {bands}")

        if "inf" in bands[:-1]:
            raise ValueError(f"'inf' is allowed as the last edge only, got {bands}")

        return bands

    @model_validator(mode="after")
    def _check_profile_index_per_band(self) -> Self:
        """A list ``profile_index`` has one entry per band."""
        nbands = len(self.bands) - 1
        if isinstance(self.profile_index, list) and len(self.profile_index) != nbands:
            raise ValueError(
                f"profile_index has {len(self.profile_index)} entries but there are "
                f"{nbands} bands; give one per band or a single value"
            )

        return self
