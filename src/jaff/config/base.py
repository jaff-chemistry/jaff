"""Shared building blocks of the module-owned config models.

Each module that is configured from ``jaff.toml`` / ``jaffgen.toml`` defines
its model in its own ``config.py``, deriving from :class:`ConfigModel`; a
parent model lists its children's models as fields.  Every model rejects
unknown keys, is strict about types (a TOML string is never coerced into a
number or a bool) and is frozen once validated.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Annotated, Any, ClassVar, List, Tuple, TypeVar

from pydantic import BaseModel, BeforeValidator, ConfigDict, ValidationError


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
    instances are frozen.

    A model whose keys depend on its ``type`` sets ``_atomic = True``: a table of
    it is replaced whole, never merged with another (see :func:`overlay` and
    :func:`~jaff.config.loader.atomic_tables`).
    """

    _atomic: ClassVar[bool] = False

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
        atomic = isinstance(value, ConfigModel) and value._atomic
        if isinstance(value, BaseModel) and current is not None and not atomic:
            value = overlay(current, value)

        update[name] = value

    return base.model_copy(update=update)


#: A key path into a config table, e.g. ``("network", "reactions", "photo")``.
KeyPath = Tuple[str, ...]

_BARE_KEY = re.compile(r"[A-Za-z0-9_-]+")

#: Plain wording for the pydantic error types a config user meets most.
_MESSAGES = {"extra_forbidden": "unknown key", "missing": "required key missing"}


def dotted(path: Tuple[str | int, ...]) -> str:
    """Write a key path as in TOML: ``network.reactions.reaction."H+.e-__H"``.

    Keys that are not bare TOML keys are quoted; list positions are ``[i]``.
    """
    text = ""
    for part in path:
        if isinstance(part, int):
            text += f"[{part}]"
            continue

        key = part if _BARE_KEY.fullmatch(part) else f'"{part}"'
        text += f".{key}" if text else key

    return text


def format_validation_errors(
    error: ValidationError, prefix: Tuple[str | int, ...] = ()
) -> List[str]:
    """Return one ``"<key path>: <message>"`` line per problem in *error*.

    Parameters
    ----------
    error : ValidationError
        Error raised validating a config model.
    prefix : tuple, optional
        Key path of that model within the config, prepended to every path.

    Returns
    -------
    list of str
        The problems; a validator's message spanning several lines gives one
        line each, under the same path.
    """
    lines: List[str] = []
    for problem in error.errors(include_url=False):
        path = dotted((*prefix, *problem["loc"]))
        message = _MESSAGES.get(problem["type"], problem["msg"])
        message = message.removeprefix("Value error, ")
        lines.extend(f"{path}: {line}" if path else line for line in message.splitlines())

    return lines
