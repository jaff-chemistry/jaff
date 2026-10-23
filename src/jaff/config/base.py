"""Shared building blocks of the module-owned config models.

Each module that is configured from ``jaff.toml`` / ``jaffgen.toml`` defines
its model in its own ``config.py``, deriving from :class:`ConfigModel`; a
parent model lists its children's models as fields.  Every model rejects
unknown keys, is strict about types (a TOML string is never coerced into a
number or a bool) and is frozen once validated.
"""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Any, TypeVar

from pydantic import BaseModel, BeforeValidator, ConfigDict


def lower(value: Any) -> Any:
    """Lower-case *value* when it is a string, so choices are case-insensitive."""
    return value.lower() if isinstance(value, str) else value


LowerStr = Annotated[str, BeforeValidator(lower)]


def _to_path(value: Any) -> Any:
    """Turn a string into a :class:`~pathlib.Path`; leave anything else as is."""
    return Path(value) if isinstance(value, str) else value


#: A path written as a string; relative paths are resolved by the loader.
PathLike = Annotated[Path, BeforeValidator(_to_path)]

_ModelT = TypeVar("_ModelT", bound=BaseModel)


class ConfigModel(BaseModel):
    """Base of every config model: unknown keys are errors, types are strict and
    instances are frozen."""

    model_config = ConfigDict(extra="forbid", strict=True, frozen=True)


def overlay(base: _ModelT, over: BaseModel | None) -> _ModelT:
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
            value = overlay(current, value)

        update[name] = value

    return base.model_copy(update=update)
