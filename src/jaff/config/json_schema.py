"""JSON Schemas of ``jaff.toml`` and ``jaffgen.toml``, generated from the config models.

Editors with a TOML language server (Taplo, Tombi) use them to check a config
file while it is written.  The schemas are committed under :data:`SCHEMA_DIR`;
regenerate them after changing a config model with::

    python -m jaff.config.json_schema

``python -m jaff.config.json_schema --path jaff`` (or ``jaffgen``) prints the
installed schema's path, for a ``#:schema`` line at the top of a config file.
TOML has no null, so an optional key is simply absent: the schemas never offer
``null`` as a value.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List

from pydantic import BaseModel

#: Directory of the committed schema files (``<name>.schema.json``).
SCHEMA_DIR: Path = Path(__file__).resolve().parent / "schemas"


def json_schema(model: type[BaseModel]) -> Dict[str, Any]:
    """Return the JSON Schema (draft 2020-12) of the config file *model*."""
    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        **_without_null(model.model_json_schema()),
    }


def _without_null(node: Any) -> Any:
    """Copy a JSON schema without ``null``: TOML has no null value.

    An ``anyOf`` alternative ``{"type": "null"}`` is dropped (an ``anyOf`` left
    with one alternative is replaced by it) and a ``null`` default removed, so
    an optional key reads as its type and editors report errors inside it.
    """
    if isinstance(node, list):
        return [_without_null(item) for item in node]

    if not isinstance(node, dict):
        return node

    node = {
        key: _without_null(value)
        for key, value in node.items()
        if not (key == "default" and value is None)
    }
    if "anyOf" in node:
        options: List[Any] = [o for o in node["anyOf"] if o != {"type": "null"}]
        if len(options) == 1:
            rest = {key: value for key, value in node.items() if key != "anyOf"}
            node = {**options[0], **rest}
        else:
            node["anyOf"] = options

    return node


def _file_models() -> Dict[str, type[BaseModel]]:
    """Return the model of each config file, keyed by file name stem."""
    # Imported here: the models live above this package in the import graph.
    from ..cli.jaffgen.config import JaffgenTomlConfig
    from ..core.network.config import JaffTomlConfig

    return {"jaff": JaffTomlConfig, "jaffgen": JaffgenTomlConfig}


def schema_files() -> Dict[str, Dict[str, Any]]:
    """Return the schema of each config file, keyed by file name stem."""
    return {name: json_schema(model) for name, model in _file_models().items()}


def write_schemas(directory: Path = SCHEMA_DIR) -> None:
    """Write ``<name>.schema.json`` for each config file into *directory*."""
    directory.mkdir(parents=True, exist_ok=True)
    for name, schema in schema_files().items():
        text = json.dumps(schema, indent=2, ensure_ascii=False) + "\n"
        (directory / f"{name}.schema.json").write_text(text, encoding="utf-8")


def main(argv: List[str] | None = None) -> None:
    """Write the schemas, or with ``--path NAME`` print one's installed path."""
    parser = argparse.ArgumentParser(prog="python -m jaff.config.json_schema")
    parser.add_argument(
        "--path",
        choices=["jaff", "jaffgen"],
        help="print the path of the schema of jaff.toml or jaffgen.toml",
    )
    args = parser.parse_args(argv)
    if args.path is not None:
        print(SCHEMA_DIR / f"{args.path}.schema.json")
        return

    write_schemas()


if __name__ == "__main__":
    main()
