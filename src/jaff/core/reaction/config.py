# ABOUTME: Config models of the reactions: [network.reactions] global settings,
# ABOUTME: per-reaction overrides and the photo-reaction (photo) table
"""Config of the reactions (``[network.reactions]``).

``[network.reactions]`` holds the settings applied to every reaction, and
``[network.reactions.reaction."<serialized>"]`` overrides them for one
reaction (see :meth:`ReactionsConfig.resolve`).  String choices are
case-insensitive and stored lower-cased.
"""

from __future__ import annotations

from typing import Annotated, Dict, Literal, Self

from pydantic import BeforeValidator, field_validator, model_validator

from ...config.base import ConfigModel, lower, overlay
from ...physics.photo_reactions.shielding.config import ShieldingConfig, check_shielding

__all__ = [
    "DuplicatePolicy",
    "PhotoConfig",
    "PiDatabase",
    "ReactionConfig",
    "ReactionsConfig",
    "TCutoff",
]


TCutoff = Annotated[Literal["clip", "extrapolate"], BeforeValidator(lower)]
DuplicatePolicy = Annotated[
    Literal["preserve-first", "preserve-last", "error"], BeforeValidator(lower)
]
PiDatabase = Annotated[Literal["norad", "verner", "leiden"], BeforeValidator(lower)]


class PhotoConfig(ConfigModel):
    """Photo-reaction settings (``[network.reactions.photo]``).

    The same keys appear in the global ``[network.reactions.photo]`` table and
    in a reaction's own ``[network.reactions.reaction."<serialized>".photo]``
    table.  A key left out of a reaction's table takes the global value; a key
    left out of the global table takes its default: ``pi_database = "norad"``,
    ``override_custom_rate = false``, ``override_custom_cross_section = false``,
    ``use_proxy_photoreaction = true`` and no ``shielding``.

    Attributes
    ----------
    pi_database : str or None
        Photoionisation cross-section database: ``"norad"``, ``"verner"`` or
        ``"leiden"``.
    override_custom_rate : bool or None
        Replace a reaction's ``.jfunc`` ``chemRate<N>`` rate with one computed
        from its cross section over the radiation bands.  Has no effect without
        a ``[network.radiation]`` table: the ``chemRate`` is then always kept.
    override_custom_cross_section : bool or None
        Ignore a reaction's ``.jfunc`` ``photoCrossSection<N>`` and take its
        cross section from the database.  The global value only applies with a
        ``[network.radiation]`` table; a value in a reaction's own table
        applies with or without one.
    use_proxy_photoreaction : bool or None
        Look cross sections up under the proxy reaction (``CHx`` → ``CH``).
    shielding : ShieldingConfig or None
        Line shielding of the rate; without it, the rate is not shielded.
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


class ReactionConfig(ConfigModel):
    """Settings of one reaction (``[network.reactions.reaction."<serialized>"]``).

    ``None`` means unset: the value comes from the global settings.

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
    """

    duplicate_policy: DuplicatePolicy | None = None
    t_cutoff: TCutoff | None = None
    photo: PhotoConfig | None = None


class ReactionsConfig(ReactionConfig):
    """Settings of the reactions (``[network.reactions]``).

    Each key applies to every reaction, unless a reaction's own table
    (``reaction."<serialized>"``) sets it.

    Attributes
    ----------
    reaction : dict[str, ReactionConfig]
        Per-reaction settings keyed by serialized reaction; a key set there
        wins over the same key here.
    """

    duplicate_policy: DuplicatePolicy = "preserve-first"
    t_cutoff: TCutoff = "clip"
    photo: PhotoConfig = _DEFAULT_PHOTO
    reaction: Dict[str, ReactionConfig] = {}

    @field_validator("photo", mode="after")
    @classmethod
    def _fill_photo_defaults(cls, photo: PhotoConfig) -> PhotoConfig:
        """Complete a partially given global ``photo`` table with the defaults."""
        return overlay(_DEFAULT_PHOTO, photo)

    @model_validator(mode="after")
    def _check_shielding_tables(self) -> Self:
        """Check every shielding table against its plugin and reaction."""
        check_shielding(self.photo.shielding, None)
        for serialized, entry in self.reaction.items():
            if entry.photo is not None:
                check_shielding(entry.photo.shielding, serialized)

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

        return overlay(base, self.reaction.get(serialized))
