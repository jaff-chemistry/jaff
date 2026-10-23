"""JSON Schemas of ``jaff.toml`` and ``jaffgen.toml``, generated from the config models.

Editors with a TOML language server (Taplo, Tombi) use them to check a config
file while it is written.  The schemas are committed under :data:`SCHEMA_DIR`;
regenerate them after changing a config model with::

    python -m jaff.config.json_schema
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict

from pydantic import BaseModel

#: Directory of the committed schema files (``<name>.schema.json``).
SCHEMA_DIR: Path = Path(__file__).resolve().parent / "schemas"


def json_schema(model: type[BaseModel]) -> Dict[str, Any]:
    """Return the JSON Schema (draft 2020-12) of the config file *model*."""
    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        **model.model_json_schema(),
    }


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


if __name__ == "__main__":
    write_schemas()
