"""Pydantic models describing the ``jaff.toml`` / ``jaffgen.toml`` config schema.

The models mirror the class structure they configure: ``[network.reactions]``
holds the settings applied to every reaction, and
``[network.reactions.reaction."<serialized>"]`` overrides them for one
reaction.  Every model rejects unknown keys and is frozen once validated.

String choices (``t_cutoff``, ``pi_database``, ...) are case-insensitive and
stored lower-cased.
"""

from __future__ import annotations

from typing import Annotated, Any, Dict, Literal, TypeVar

from pydantic import BaseModel, BeforeValidator, ConfigDict, field_validator


def _lower(value: Any) -> Any:
    """Lower-case *value* when it is a string, so choices are case-insensitive."""
    return value.lower() if isinstance(value, str) else value


TCutoff = Annotated[Literal["clip", "extrapolate"], BeforeValidator(_lower)]
DuplicatePolicy = Annotated[
    Literal["preserve-first", "preserve-last", "error"], BeforeValidator(_lower)
]
PiDatabase = Annotated[Literal["norad", "verner", "leiden"], BeforeValidator(_lower)]

_ModelT = TypeVar("_ModelT", bound=BaseModel)


class _Config(BaseModel):
    """Base of every config model: unknown keys are errors, instances are frozen."""

    model_config = ConfigDict(extra="forbid", frozen=True)


def _overlay(base: _ModelT, over: BaseModel | None) -> _ModelT:
    """Return *base* with every key set (not ``None``) in *over* applied.

    Nested models are overlaid field by field rather than replaced, so an
    override that sets one nested key keeps the others from *base*.
    """
    if over is None:
        return base

    update: dict[str, Any] = {}
    for name in type(over).model_fields:
        value = getattr(over, name)
        if value is None:
            continue

        current = getattr(base, name)
        if isinstance(value, BaseModel) and current is not None:
            value = _overlay(current, value)

        update[name] = value

    return base.model_copy(update=update)


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

    ``None`` means unset: the value comes from the global settings.
    """

    pi_database: PiDatabase | None = None
    override_custom_rate: bool | None = None
    override_custom_cross_section: bool | None = None
    use_proxy_photoreaction: bool | None = None


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
        base = ReactionConfig.model_validate(self.model_dump(exclude={"reaction"}))

        return _overlay(base, self.reaction.get(serialized))
