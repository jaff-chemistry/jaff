# ABOUTME: Shared helpers for the jaffgen/jaffx CLIs: option→config-key mapping,
# ABOUTME: --set parsing, --funcfile parsing and error reporting for commands.
"""Shared argument helpers for the JAFF command-line entry points."""

from __future__ import annotations

import functools
import json
import tomllib
from enum import Enum
from typing import Any, Callable, Dict, List, Tuple, TypeVar, cast, get_args

import typer

from ..config.loader import deep_merge
from ..core.reaction import config as reaction_config
from ..errors import ParserError

#: Command-line option (its ``args`` attribute) → config key path it sets.
CLI_KEYS: Dict[str, Tuple[str, ...]] = {
    "template": ("jaffgen", "template"),
    "indir": ("jaffgen", "input_dir"),
    "files": ("jaffgen", "input_files"),
    "outdir": ("jaffgen", "output_dir"),
    "network": ("jaffgen", "network_file"),
    "network_config": ("jaffgen", "network_config"),
    "echo": ("jaffgen", "echo"),
    "label": ("network", "label"),
    "funcfile": ("network", "funcfile"),
    "expand_nuclei": ("network", "expand_nuclei"),
    "strict": ("network", "strict"),
    "duplicate_policy": ("network", "reactions", "duplicate_policy"),
}

#: The :data:`CLI_KEYS` setting ``[network]`` keys (the options a network tool may
#: take; :func:`cli_layer` skips any option a tool does not define).
NETWORK_KEYS: Dict[str, Tuple[str, ...]] = {
    option: path for option, path in CLI_KEYS.items() if path[0] == "network"
}


_CommandT = TypeVar("_CommandT", bound=Callable[..., Any])


def reports_errors(command: _CommandT) -> _CommandT:
    """Make a CLI command print input errors instead of a traceback.

    A :class:`~jaff.errors.ParserError` (a bad config, network or option) or a
    missing file is printed to stderr and the command exits with status 1;
    any other exception is a bug and propagates with its traceback.

    Parameters
    ----------
    command : callable
        The typer command function; its signature is kept.

    Returns
    -------
    callable
        The wrapped command.
    """

    @functools.wraps(command)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        try:
            return command(*args, **kwargs)
        except (ParserError, FileNotFoundError) as error:
            typer.echo(f"Error: {error}", err=True)
            raise typer.Exit(1) from None

    return cast(_CommandT, wrapper)


#: ``--duplicate-policy`` choices of the JAFF CLI tools: built from the config's
#: allowed values (:data:`~jaff.core.reaction.config.DuplicatePolicy`), so the two
#: cannot drift apart.
DuplicatePolicyChoice = Enum(  # type: ignore[misc]
    "DuplicatePolicyChoice",
    {
        value.replace("-", "_"): value
        for value in get_args(get_args(reaction_config.DuplicatePolicy)[0])
    },
    type=str,
)


def funcfile_arg(value: str) -> bool | str:
    """
    Parse the ``--funcfile`` flag, mapping ``true``/``false`` onto booleans.

    A command line cannot express a boolean literal the way a ``jaffgen.toml``
    can, so ``--funcfile false`` (case-insensitive) is the CLI spelling for
    ``funcfile=False`` — skip auxiliary-function loading entirely — and
    ``--funcfile true`` for ``funcfile=True``, which scans the network
    directory.  ``--funcfile true`` (or ``--set network.funcfile=true``)
    overrides a ``funcfile`` path set in a config file back to a directory
    scan.

    Parameters
    ----------
    value : str
        Raw value supplied on the command line.

    Returns
    -------
    bool | str
        ``False`` to skip, ``True`` to scan, otherwise *value* unchanged.
    """
    lowered = value.lower()
    if lowered == "false":
        return False
    if lowered == "true":
        return True

    return value


def cli_layer(args: Any, keys: Dict[str, Tuple[str, ...]] = CLI_KEYS) -> Dict[str, Any]:
    """Return the config layer of the command-line options in *args*.

    Each option given (an attribute of *args* that is not ``None``) sets the
    config key path :data:`CLI_KEYS` maps it to; options not given, or absent
    from *args*, set nothing.

    Parameters
    ----------
    args : Any
        Parsed command-line options, one attribute per option.
    keys : dict, optional
        Option name → config key path (default :data:`CLI_KEYS`).

    Returns
    -------
    dict
        The layer, shaped like a ``jaffgen.toml``.
    """
    layer: Dict[str, Any] = {}
    for option, path in keys.items():
        value = getattr(args, option, None)
        if value is None:
            continue

        table = layer
        for key in path[:-1]:
            table = table.setdefault(key, {})

        table[path[-1]] = value

    return layer


def parse_set(items: List[str]) -> Dict[str, Any]:
    """Return the config layer of ``--set key.path=value`` items, later ones winning.

    The key path and the value are written as in TOML (a reaction key is
    quoted: ``network.reactions.reaction."H2._PHOTON__H.H".t_cutoff=clip``); a
    value that is not a TOML value is taken as a string, so
    ``network.label=gow`` needs no quotes.

    Parameters
    ----------
    items : list of str
        The ``--set`` values, in command-line order.

    Returns
    -------
    dict
        The layer, shaped like a config file.

    Raises
    ------
    ParserError
        If an item has no ``=`` or its key path is not a TOML dotted key.
    """
    return deep_merge(*(_set_item(item) for item in items))


def _set_item(item: str) -> Dict[str, Any]:
    """Return the config layer of one ``--set`` item; see :func:`parse_set`."""
    key, sep, value = item.partition("=")
    if not sep or not key.strip():
        raise ParserError(f"Invalid --set {item!r}: expected key.path=value")

    # The value as TOML, else as a plain string.
    for text in (f"{key} = {value}", f"{key} = {json.dumps(value)}"):
        try:
            return tomllib.loads(text)
        except tomllib.TOMLDecodeError:
            continue

    raise ParserError(f"Invalid --set {item!r}: {key!r} is not a valid key path")


def command_line_layer(
    args: Any, keys: Dict[str, Tuple[str, ...]] = CLI_KEYS
) -> Dict[str, Any]:
    """Return the config layer of a command line: its options, then ``--set``.

    Parameters
    ----------
    args : Any
        Parsed command-line options; ``args.set`` (if any) holds the
        ``--set`` items, which win over the named options.
    keys : dict, optional
        Option name → config key path (default :data:`CLI_KEYS`).

    Returns
    -------
    dict
        The layer, shaped like a config file.
    """
    return deep_merge(cli_layer(args, keys), parse_set(getattr(args, "set", None) or []))
