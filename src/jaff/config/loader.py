"""Loading of the ``jaff.toml`` / ``jaffgen.toml`` config.

Every config source (a file, the command line, a Python ``dict``) is a *layer*:
a nested ``dict`` in the shape of the config file.  The layers are merged with
:func:`deep_merge`, the later one winning, and the result is validated against
the module-owned config models.
"""

from __future__ import annotations

import tomllib
from pathlib import Path
from types import UnionType
from typing import (
    Annotated,
    Any,
    Callable,
    Dict,
    List,
    Mapping,
    Sequence,
    TypeVar,
    Union,
    get_args,
    get_origin,
)

import tomli_w
from pydantic import BaseModel, ValidationError

from ..errors import ParserError
from .base import KeyPath, format_validation_errors

__all__ = [
    "atomic_tables",
    "deep_merge",
    "dump_config",
    "load_config",
    "load_layer",
    "resolve_paths",
]

_ModelT = TypeVar("_ModelT", bound=BaseModel)


def dump_config(config: BaseModel, schema: Path | None = None) -> str:
    """Write *config* as TOML text: every key, defaults included.

    Paths are written as strings and shielding tables flat, as in a config
    file; keys whose value is ``None`` (unset optionals) and empty tables are
    omitted, as TOML has no null; a ``[jaffgen]`` section comes first.  The
    text loads back into an equal config.

    Parameters
    ----------
    config : BaseModel
        A validated config (e.g. a ``JaffgenTomlConfig``).
    schema : Path or None, optional
        JSON Schema to name in a leading ``#:schema`` line, for editors.

    Returns
    -------
    str
        The TOML text.
    """
    data = _drop_none(config.model_dump(mode="json"))
    if "jaffgen" in data:
        data = {"jaffgen": data.pop("jaffgen"), **data}

    header = f"#:schema {schema}\n\n" if schema is not None else ""

    return header + tomli_w.dumps(data)


def _drop_none(value: Any) -> Any:
    """Copy *value*, leaving out ``None`` values and empty tables."""
    if isinstance(value, dict):
        items = {key: _drop_none(item) for key, item in value.items() if item is not None}

        return {key: item for key, item in items.items() if item != {}}

    if isinstance(value, list):
        return [_drop_none(item) for item in value]

    return value


def load_config(layers: Sequence[Mapping[str, Any]], model: type[_ModelT]) -> _ModelT:
    """Merge config *layers* and validate the result against *model*.

    Parameters
    ----------
    layers : Sequence[Mapping[str, Any]]
        Layers, lowest precedence first (paths already resolved).
    model : type[BaseModel]
        Model of the whole config (e.g. ``JaffgenTomlConfig``); it also says
        which tables are replaced whole between layers (:func:`atomic_tables`).

    Returns
    -------
    BaseModel
        The validated config.

    Raises
    ------
    ParserError
        Listing every problem as ``<key path>: <message>``.
    """
    merged = deep_merge(*layers, atomic=atomic_tables(model))
    try:
        return model.model_validate(merged)
    except ValidationError as e:
        problems = "\n  ".join(format_validation_errors(e))
        raise ParserError(f"Invalid config:\n  {problems}") from e


def load_layer(path: Path, model: type[BaseModel]) -> Dict[str, Any]:
    """Read the TOML file *path* as a layer, its paths resolved against its directory.

    Only its top-level sections are checked here, so an error names the file;
    the rest is validated once all layers are merged.

    Parameters
    ----------
    path : Path
        TOML config file.
    model : type[BaseModel]
        Model of the whole file (e.g. ``JaffTomlConfig``), locating its paths.

    Returns
    -------
    dict
        The layer.

    Raises
    ------
    FileNotFoundError
        If *path* does not exist.
    ParserError
        If *path* is not valid TOML or has a top-level section *model* does not
        hold (e.g. ``[jaffgen]`` in a ``jaff.toml``).
    """
    with open(path, "rb") as f:
        try:
            layer = tomllib.load(f)
        except tomllib.TOMLDecodeError as e:
            raise ParserError(f"Invalid TOML in {path}: {e}") from e

    unknown = [key for key in layer if key not in model.model_fields]
    if unknown:
        raise ParserError(
            f"{path}: unknown section(s) {', '.join(unknown)}; "
            f"this file holds only {', '.join(model.model_fields)}"
        )

    return resolve_paths(layer, model, path.parent)


def resolve_paths(
    layer: Mapping[str, Any], model: type[BaseModel], base_dir: Path
) -> Dict[str, Any]:
    """Return a copy of *layer* with its relative paths resolved against *base_dir*.

    A path is a string (or :class:`~pathlib.Path`) under a key the *model*
    types as a path (a ``PathLike`` field, or a list of them).  Anything else -- including unknown keys and
    values of the wrong type, left for validation to report -- is copied as is.

    Parameters
    ----------
    layer : Mapping[str, Any]
        Layer to resolve.
    model : type[BaseModel]
        Model of the whole layer (e.g. ``JaffgenTomlConfig``).
    base_dir : Path
        Directory relative paths are resolved against: the file's directory,
        or the working directory for the command line.

    Returns
    -------
    dict
        The resolved layer.
    """
    return _resolve(layer, [model], base_dir)


def _resolve(value: Any, types: List[Any], base_dir: Path) -> Any:
    """Copy *value*, typed by any of *types*, resolving the paths in it."""
    if isinstance(value, Mapping):
        return {
            key: _resolve(item, [c for t in types for c in _children(t, key)], base_dir)
            for key, item in value.items()
        }

    if isinstance(value, list):
        items = [
            m for t in types if get_origin(t) is list for m in _members(get_args(t)[0])
        ]

        return [_resolve(item, items, base_dir) for item in value]

    if isinstance(value, (str, Path)) and Path in types:
        path = Path(value)

        return path if path.is_absolute() else (base_dir / path).resolve()

    return value


def atomic_tables(model: type[BaseModel]) -> Callable[[KeyPath], bool]:
    """Return the :func:`deep_merge` ``atomic`` check derived from *model*.

    A key path is atomic when it leads, through *model*'s fields, per-key
    ``dict`` entries and union members, to a config model flagged
    ``_atomic = True`` (e.g. a ``shielding`` table or the ``eos`` table, whose
    keys depend on their ``type``).

    Parameters
    ----------
    model : type[BaseModel]
        Model of a whole config file (e.g. ``JaffgenTomlConfig``).

    Returns
    -------
    Callable[[KeyPath], bool]
        ``True`` for the key path of an atomic table.
    """

    def atomic(path: KeyPath) -> bool:
        nodes: List[Any] = [model]
        for key in path:
            nodes = [child for node in nodes for child in _children(node, key)]

        return any(getattr(node, "_atomic", False) for node in nodes)

    return atomic


def _children(node: Any, key: str) -> List[Any]:
    """Return the types a table of type *node* can hold under *key*."""
    if isinstance(node, type) and issubclass(node, BaseModel):
        field = node.model_fields.get(key)

        return [] if field is None else _members(field.annotation)

    if get_origin(node) is dict:
        return _members(get_args(node)[1])

    return []


def _members(annotation: Any) -> List[Any]:
    """Return the types *annotation* stands for, unwrapping ``Annotated`` and unions."""
    origin = get_origin(annotation)
    if origin is Annotated:
        return _members(get_args(annotation)[0])

    if origin is Union or origin is UnionType:
        return [member for arg in get_args(annotation) for member in _members(arg)]

    return [annotation]


def _never_atomic(path: KeyPath) -> bool:
    """Default for :func:`deep_merge`: every table merges."""
    return False


def deep_merge(
    *layers: Mapping[str, Any],
    atomic: Callable[[KeyPath], bool] = _never_atomic,
) -> Dict[str, Any]:
    """Merge config *layers*, a later layer winning over an earlier one.

    Tables present in both layers merge key by key; any other value (a list, a
    scalar, or a table meeting a non-table) is replaced by the later layer.  A
    table for which ``atomic(path)`` is true is replaced whole instead of
    merged.  The inputs are not modified.

    Parameters
    ----------
    *layers : Mapping[str, Any]
        Layers, lowest precedence first.
    atomic : Callable[[KeyPath], bool], optional
        Called with the key path of each table present in two layers; ``True``
        replaces the earlier table with the later one.  By default every table
        merges.

    Returns
    -------
    dict
        The merged layer.
    """
    merged: Dict[str, Any] = {}
    for layer in layers:
        merged = _merge(merged, layer, (), atomic)

    return merged


def _merge(
    base: Mapping[str, Any],
    over: Mapping[str, Any],
    path: KeyPath,
    atomic: Callable[[KeyPath], bool],
) -> Dict[str, Any]:
    """Return a fresh merge of *over* into *base*; see :func:`deep_merge`."""
    merged = {key: _copy(value) for key, value in base.items()}

    for key, value in over.items():
        current = merged.get(key)
        key_path = (*path, key)
        if (
            isinstance(current, Mapping)
            and isinstance(value, Mapping)
            and not atomic(key_path)
        ):
            merged[key] = _merge(current, value, key_path, atomic)
        else:
            merged[key] = _copy(value)

    return merged


def _copy(value: Any) -> Any:
    """Copy nested tables and lists, so the merge never shares them with a layer."""
    if isinstance(value, Mapping):
        return {key: _copy(item) for key, item in value.items()}

    if isinstance(value, list):
        return [_copy(item) for item in value]

    return value
