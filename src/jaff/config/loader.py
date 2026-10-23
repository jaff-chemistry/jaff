"""Loading of the ``jaff.toml`` / ``jaffgen.toml`` config.

Every config source (a file, the command line, a Python ``dict``) is a *layer*:
a nested ``dict`` in the shape of the config file.  The layers are merged with
:func:`deep_merge`, the later one winning, and the result is validated against
the module-owned config models.
"""

from __future__ import annotations

from types import UnionType
from typing import (
    Annotated,
    Any,
    Callable,
    Dict,
    List,
    Mapping,
    Tuple,
    Union,
    get_args,
    get_origin,
)

from pydantic import BaseModel

#: A key path into a layer, e.g. ``("network", "reactions", "photo")``.
KeyPath = Tuple[str, ...]


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
