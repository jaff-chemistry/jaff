# ABOUTME: ConfigTable: turns one [[jaffgen.table]] config entry into an HDF5Dict or
# ABOUTME: DataFrame (path remapping, composite datasets, computed attributes).
"""
Configuration table engine for the JAFF CLI.

This module implements :class:`ConfigTable`, which reads a single
``[[jaffgen.table]]`` block from a ``jaffgen.toml`` configuration file and converts the
described data table into either an :class:`~jaff.types.HDF5Dict` or a
:class:`pandas.DataFrame`, ready to be written to disk by the caller.

A ``[[jaffgen.table]]`` block has two sub-sections:

``source`` (required)
    Describes where to read the input data from.  Supported formats are HDF5
    and CSV.  The special path value ``"default"`` expands to
    ``<network_dir>/<network_stem>.hdf5``.

``target`` (optional)
    Describes the output format and path.  Defaults to the same format and
    filename as the source.  HDF5 targets additionally support per-dataset
    path remapping and attribute injection via a dot-notation syntax
    (``<source_path>.<property>``).

Attribute notation
------------------
Within a ``[[jaffgen.table]]`` block's ``target`` HDF5 dataset entry the
``attrs`` key accepts mappings of the form:

```toml
[jaffgen.table.target."/some/dataset"]
attrs = { "my_attr" = "/other/dataset.max" }
```

This reads the ``max`` property of ``/other/dataset`` and stores it as the
HDF5 attribute ``my_attr`` on ``/some/dataset``.  Recognised property names
are: ``max``, ``min``, ``mean``, ``median``, ``length`` (``median`` currently
raises :class:`AttributeError`; see :attr:`ConfigTable.attr_dict`).

A value may also be a **list**, resolved element-wise, so computed references
and literals can be mixed:

```toml
[jaffgen.table.target."/co".attrs]
Ndim = 2                                  # literal
Nx   = ["/co/TCO.length", "/co/x.length"] # list of computed references
spacing = ["log", "log"]                  # list of literals
```
"""

import copy
from functools import cached_property
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ...common._helper import CSV_EXTENSIONS, HDF_EXTENSIONS
from ...types import HDF5Dict


class ConfigTable:
    """
    Parse a single ``[[jaffgen.table]]`` TOML entry and produce structured output.

    :class:`ConfigTable` bridges the TOML configuration layer and the JAFF
    data drivers.  It validates source/target paths, loads the source data
    (HDF5 or CSV), applies any requested column/path remapping, injects HDF5
    attributes, and returns a ready-to-write object.

    Parameters
    ----------
    table_dict : dict[str, Any]
        The parsed ``[[jaffgen.table]]`` dictionary from the TOML config.  Must
        contain at least a ``"source"`` key.
    file : Path
        Path to the ``jaffgen.toml`` file; relative source paths resolve
        against its directory, and it is named in error messages.
    network_file : Path
        Path to the network file.  Its parent directory and stem are
        used to resolve the ``"default"`` source path.

    Raises
    ------
    RuntimeError
        If the network directory does not contain a default data file when
        ``source.path = "default"`` is used.
    KeyError
        If the ``"source"`` key is missing from *table_dict*.
    ValueError
        If the source or target file extension is not a recognised HDF5 or
        CSV extension.
    """

    def __init__(self, table_dict: dict[str, Any], file: Path, network_file: Path):
        """Initialise and immediately validate source/target configuration.

        Parameters
        ----------
        table_dict : dict[str, Any]
            Parsed ``[[jaffgen.table]]`` entry from the TOML config.
        file : Path
            Path to the ``jaffgen.toml`` file (relative source paths resolve
            against its directory; also used in error messages).
        network_file : Path
            Path to the network file; its parent and stem resolve
            the ``"default"`` source path alias.

        Raises
        ------
        RuntimeError
            If the network directory does not contain a default data file.
        KeyError
            If the ``"source"`` key is missing from *table_dict*.
        ValueError
            If the source or target file extension is unsupported.
        """
        self.config: dict[str, Any] = table_dict
        # Explicit relative source paths resolve against the config file's dir.
        self.config_dir = file.parent
        self.network_dir = network_file.parent
        self.default_data: Path = network_file.with_suffix(".hdf5")
        if (
            table_dict.get("source", {}).get("path") == "default"
            and not self.default_data.exists()
        ):
            raise RuntimeError(
                f"{self.network_dir} doesn't contain a default data file "
                f"({self.default_data.name})"
            )

        if "source" not in table_dict:
            raise KeyError(f"source must be specified in {file}")
        self.source_config: dict[str, Any] = self.config["source"]
        self.target_config: dict[str, Any] = self.config.get("target", {})

        self.source_props: dict[str, Any] = {}
        self.target_props: dict[str, Any] = {}
        self.__set_source_props()
        self.__set_target_props()

        # Load the source data tree immediately so it is available during parse().
        self.source_tree: dict[str, Any] | HDF5Dict = self.__get_source_tree()

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def parse(self) -> HDF5Dict | pd.DataFrame | None:
        """
        Convert the source data to the target format.

        Dispatches to :meth:`__parse_to_hdf5` or :meth:`__parse_to_csv`
        based on ``target_props["type"]``.

        Returns
        -------
        HDF5Dict or pandas.DataFrame or None
            Structured data ready to be written to disk, or ``None`` if the
            target type is unsupported.
        """
        if self.target_props["type"] == "hdf5":
            return self.__parse_to_hdf5()
        elif self.target_props["type"] == "csv":
            return self.__parse_to_csv()

    # ------------------------------------------------------------------
    # HDF5 output
    # ------------------------------------------------------------------

    def __parse_to_hdf5(self) -> HDF5Dict:
        """
        Build an :class:`~jaff.types.HDF5Dict` for the target HDF5 file.

        Handles two source types:

        * **HDF5 → HDF5**: Copies the source tree, applying any path
          remappings declared under HDF5 path keys (keys starting with
          ``"/"``).
        * **CSV → HDF5**: Wraps each requested column as a ``"linear"``
          HDF5 dataset under the declared path.

        After path remapping, any ``attrs`` entries are resolved via
        :meth:`__resolve_attr_value` and injected onto their target path.

        Returns
        -------
        HDF5Dict
            Nested dictionary representing the target HDF5 structure.

        Raises
        ------
        ValueError
            If an attribute reference uses an invalid dot-notation path.
        """
        # Collect all HDF5 path-level overrides from the target config
        # (keys that start with "/" are interpreted as dataset paths).
        target_hdf_tree = {
            k: v for k, v in self.target_config.items() if k.startswith("/")
        }
        if self.source_props["type"] == "hdf5":
            # ``flatten()`` yields a plain dict keyed by absolute HDF5 path.
            assert isinstance(self.source_tree, dict)
            if not target_hdf_tree:
                # No remapping requested — return the source tree as-is.
                return HDF5Dict(self.source_tree)

            target_tree = copy.deepcopy(self.source_tree)
            # Remap dataset paths: move datasets from their source h5path to
            # the declared target path, removing the old key if different.
            for path, items in target_hdf_tree.items():
                if "h5path" not in items:
                    continue

                h5path = items["h5path"]
                # A list of source paths declares a *composite* dataset: several
                # source datasets folded into a single compound (named-column) or
                # N-D array dataset. A bare string keeps the single-move behaviour.
                if isinstance(h5path, list):
                    target_tree[path] = self.__build_composite(h5path, items)
                    # The folded-in source datasets are consumed; drop them so
                    # they don't also pass through at their old locations.
                    for src in map(self.__norm_path, h5path):
                        target_tree.pop(src, None)
                    continue

                target_tree[path] = self.source_tree[h5path]
                if h5path != path and h5path in target_tree:
                    target_tree.pop(h5path)

        elif self.source_props["type"] == "csv":
            assert isinstance(self.source_tree, dict)
            target_tree = HDF5Dict({})

            # Wrap each requested CSV column as a linear HDF5 dataset.
            for path, items in target_hdf_tree.items():
                if "col" not in items:
                    continue

                target_tree[path] = {
                    "_kind": "linear",  # scalar/1-D array dataset
                    "_name": items["col"],
                    "_data": self.source_tree[items["col"]],
                    "_dtype": HDF5Dict._from_np()[self.source_tree[items["col"]].dtype],
                }

        # Inject HDF5 attributes. Each value is resolved to either a computed
        # statistic (``"/path.property"`` reference), a plain literal, or a list
        # of those resolved element-wise.
        for path, items in target_hdf_tree.items():
            if "attrs" not in items:
                continue

            for var, value in items["attrs"].items():
                resolved = self.__resolve_attr_value(target_tree, value)
                target_tree[path] = {
                    **target_tree.get(path, {}),
                    "_attrs": {
                        **target_tree.get(path, {}).get("_attrs", {}),
                        var: resolved,
                    },
                }

        return HDF5Dict(target_tree)

    # ------------------------------------------------------------------
    # Composite datasets
    # ------------------------------------------------------------------

    @staticmethod
    def __norm_path(path: str) -> str:
        """Normalise a source path to an absolute HDF5 path (leading ``/``)."""
        return path if path.startswith("/") else "/" + path

    def __build_composite(self, h5path: list[str], items: dict[str, Any]) -> dict:
        """
        Fold several source datasets into one composite HDF5 dataset.

        Triggered when a ``[jaffgen.table.target."/path"]`` heading declares its
        ``h5path`` as a **list** of source paths.  The ``type`` key selects the
        storage layout (default ``"compound"``):

        * ``"compound"`` — a NumPy structured (record) array, one named column
          per source dataset.  Column names come from the ``names`` list, or
          from each source path's final component when ``names`` is absent.
        * ``"ndarray"`` — a single N-D array of shape ``(Ncols, *xlens)``, where
          ``Ncols`` is the number of source datasets and ``xlens`` are the
          lengths of the axis datasets listed under ``regrid.x``.  Each source
          column is reshaped to the axis grid and stacked along a new leading
          axis.

        Parameters
        ----------
        h5path : list of str
            Source dataset paths to combine.  Missing leading slashes are
            normalised.
        items : dict
            The target heading's sub-config (may hold ``type``, ``names``,
            ``regrid``).

        Returns
        -------
        dict
            A leaf dataset dict (``_kind``/``_data``/``_dtype``/``_attrs``).

        Raises
        ------
        ValueError
            If the source columns differ in length, if ``names`` length does
            not match, if ``type`` is unrecognised, or (``ndarray``) if
            ``regrid.x`` is missing or the reshape does not fit.
        """
        paths = [self.__norm_path(p) for p in h5path]
        leaves = [self.source_tree[p] for p in paths]

        # Every column must share a length to form a single dataset.
        lengths = {len(leaf["_data"]) for leaf in leaves}
        if len(lengths) != 1:
            raise ValueError(
                f"Composite columns must share a length; got {lengths} for {paths}"
            )

        kind = items.get("type", "compound")
        if kind == "compound":
            return self.__build_compound(paths, leaves, items)
        if kind == "ndarray":
            return self.__build_ndarray(paths, leaves, items)
        raise ValueError(
            f"Unknown composite type {kind!r} (expected 'compound' or 'ndarray')"
        )

    def __build_compound(
        self, paths: list[str], leaves: list[dict], items: dict[str, Any]
    ) -> dict:
        """Build a structured-array (compound) leaf from *leaves*.

        Column names are taken from ``items['names']`` when present, otherwise
        from each source path's final component.
        """
        names = items.get("names") or [p.rsplit("/", 1)[-1] for p in paths]
        if len(names) != len(paths):
            raise ValueError(f"names length {len(names)} != h5path length {len(paths)}")

        np_fields = [
            (name, np.asarray(leaf["_data"]).dtype) for name, leaf in zip(names, leaves)
        ]
        data = np.empty(len(leaves[0]["_data"]), dtype=np_fields)
        for name, leaf in zip(names, leaves):
            data[name] = leaf["_data"]

        return {
            "_kind": "compound",
            "_data": data,
            # Reuse each source leaf's own JAFF dtype token for the field schema.
            "_dtype": {name: leaf["_dtype"] for name, leaf in zip(names, leaves)},
            "_attrs": {},
        }

    def __build_ndarray(
        self, paths: list[str], leaves: list[dict], items: dict[str, Any]
    ) -> dict:
        """Build an ``(Ncols, *xlens)`` N-D array leaf from *leaves*.

        Axis lengths are read from the ``regrid.x`` dataset list; each source
        column is reshaped to that grid and stacked along a new leading axis.
        """
        xpaths = items.get("regrid", {}).get("x")
        if xpaths is None:
            raise ValueError(
                "ndarray composite requires a 'regrid.x' axis list to define shape"
            )
        if isinstance(xpaths, str):
            xpaths = [xpaths]

        xlens = [len(self.source_tree[self.__norm_path(x)]["_data"]) for x in xpaths]
        expected = int(np.prod(xlens))

        cols = []
        for path, leaf in zip(paths, leaves):
            arr = np.asarray(leaf["_data"])
            if arr.size != expected:
                raise ValueError(
                    f"{path}: size {arr.size} does not match axis grid {xlens} "
                    f"(product {expected})"
                )
            cols.append(arr.reshape(xlens))

        data = np.stack(cols, axis=0)  # (Ncols, *xlens)
        return {
            "_kind": "linear",
            "_data": data,
            "_dtype": HDF5Dict._from_np()[data.dtype],
            "_attrs": {},
        }

    # ------------------------------------------------------------------
    # Attribute helpers
    # ------------------------------------------------------------------

    def __resolve_attr_value(self, target_tree: dict, value: Any) -> Any:
        """
        Resolve one attribute value from the config to its stored form.

        The resolution rules are:

        * **list** — resolved element-wise (a list of resolved values).
        * **string of the form** ``"/path.property"`` where *property* is a key
          of :attr:`attr_dict` — a computed statistic read from the dataset at
          ``/path`` in *target_tree* (e.g. ``"/co/TCO.max"``).
        * **anything else** — a literal, stored verbatim (numbers, booleans, and
          plain strings such as ``"log"`` or ``"K"``).

        Parameters
        ----------
        target_tree : dict
            The in-progress target tree the references are computed against.
        value : Any
            Raw attribute value from the TOML config.

        Returns
        -------
        Any
            The resolved attribute value (scalar, literal, or list thereof).

        Raises
        ------
        ValueError
            If a value parses as a ``"/path.property"`` reference but ``/path``
            is not present in *target_tree*.
        """
        if isinstance(value, list):
            return [self.__resolve_attr_value(target_tree, v) for v in value]

        if isinstance(value, str):
            # A computed reference is "<path>.<property>" where <property> is a
            # known statistic; the "." also distinguishes it from a plain literal
            # such as "log". Anything else is stored verbatim.
            prop_path, sep, prop = value.rpartition(".")
            if sep and prop in self.attr_dict:
                if prop_path not in target_tree:
                    raise ValueError(
                        f"Attribute reference path not found in target tree: {value}"
                    )
                return self.get_attr(target_tree, prop_path, prop)

        return value

    def get_attr(self, target_tree: dict, prop_path: str, prop: str) -> Any:
        """
        Compute a statistical attribute from a dataset in *target_tree*.

        Parameters
        ----------
        target_tree : dict
            The in-progress target tree that contains *prop_path*.
        prop_path : str
            HDF5 path of the dataset whose ``_data`` is used.
        prop : str
            Name of the statistic to compute.  Must be a key of
            :attr:`attr_dict` (``"max"``, ``"min"``, ``"mean"``,
            ``"median"``, ``"length"``).

        Returns
        -------
        Any
            The computed scalar attribute value.

        Raises
        ------
        ValueError
            If the dataset at *prop_path* has no ``_data`` field.
        """
        if target_tree[prop_path].get("_data") is None:
            raise ValueError(f"Path doesn't contain any data: {prop_path}")

        return self.attr_dict[prop](target_tree[prop_path]["_data"])

    @cached_property
    def attr_dict(self) -> dict:
        """
        Mapping of attribute property names to aggregation callables.

        Each callable receives a dataset's ``_data`` (a numpy array) and
        returns a scalar.  ``"median"`` calls ``data.median()``, which numpy
        arrays do not have, so it currently raises :class:`AttributeError`.

        Returns
        -------
        dict
            Keys are ``"max"``, ``"min"``, ``"mean"``, ``"median"``,
            ``"length"``; values are single-argument lambdas.
        """
        return {
            "max": lambda data: data.max(),
            "min": lambda data: data.min(),
            "mean": lambda data: data.mean(),
            "median": lambda data: data.median(),
            "length": lambda data: data.size,
        }

    # ------------------------------------------------------------------
    # CSV output
    # ------------------------------------------------------------------

    def __parse_to_csv(self) -> pd.DataFrame | None:
        """
        Build a :class:`pandas.DataFrame` from a CSV source.

        Currently only supports a CSV source type.  Returns ``None`` for
        other source types (e.g. HDF5 → CSV is not yet implemented).

        Returns
        -------
        pandas.DataFrame or None
            DataFrame of the source CSV data, or ``None`` if the source is
            not CSV.
        """
        if self.source_props["type"] == "csv":
            return pd.DataFrame(self.source_tree)

    # ------------------------------------------------------------------
    # Property setters
    # ------------------------------------------------------------------

    def __set_source_props(self) -> None:
        """
        Resolve and validate source file properties from the config.

        Expands the ``"default"`` path alias to the network's default HDF5
        file, resolves any other relative path against the config file's
        directory, and infers the format from the file extension.  For CSV
        sources, also reads delimiter and comment-character settings.

        Stores the result in ``self.source_props``.

        Raises
        ------
        ValueError
            If the source file extension is not a supported HDF5 or CSV
            extension.
        """
        props: dict[str, Any] = {
            # "default" expands to <network_dir>/<network_stem>.hdf5
            "path": self.default_data
            if self.source_config["path"] == "default"
            else self.config_dir / self.source_config["path"]
        }
        if props["path"].suffix.lower() not in HDF_EXTENSIONS + CSV_EXTENSIONS:
            raise ValueError(f"Unsupported source format found: {props['path']}")

        props["type"] = (
            "csv" if props["path"].suffix.lower() in CSV_EXTENSIONS else "hdf5"
        )
        if props["type"] == "csv":
            # Read CSV parsing options; JAFF defaults are a space delimiter and
            # "#" comments (not pandas read_csv's "," and no comment character).
            props["delimiter"] = self.source_config.get("delimiter", " ")
            props["comment"] = self.source_config.get("comment", "#")

        self.source_props = props

    def __set_target_props(self) -> None:
        """
        Resolve and validate target file properties from the config.

        Defaults to the source file's name (and therefore format) when no
        explicit target path is given.  Reads format-specific options: the CSV
        delimiter/comment, and the HDF5 ``default_group`` (stored in
        ``target_props`` but not currently used when writing).

        Stores the result in ``self.target_props``.

        Raises
        ------
        ValueError
            If the target file extension is not a supported HDF5 or CSV
            extension.
        """
        props: dict[str, Any] = {
            # Default to the same filename as the source if no explicit path given.
            "path": Path(self.target_config.get("path", self.source_props["path"].name))
        }
        if props["path"].suffix.lower() not in HDF_EXTENSIONS + CSV_EXTENSIONS:
            raise ValueError(f"Unsupported target format found: {props['path']}")

        props["type"] = (
            "csv" if props["path"].suffix.lower() in CSV_EXTENSIONS else "hdf5"
        )
        if props["type"] == "hdf5":
            props["default_group"] = self.target_config.get("default_group", "/")
        elif props["type"] == "csv":
            props["delimiter"] = self.target_config.get("delimiter", " ")
            props["comment"] = self.target_config.get("comment", "#")

        self.target_props = props

    # ------------------------------------------------------------------
    # Source data loading
    # ------------------------------------------------------------------

    def __get_source_tree(self) -> dict[str, Any]:
        """
        Load the source data into a flat dictionary.

        For HDF5 sources, returns the flattened :class:`~jaff.types.HDF5Dict`
        where every key is an absolute HDF5 path.  For CSV sources, returns a
        plain ``{column_name: numpy_array}`` dictionary.

        Returns
        -------
        dict[str, Any]
            Loaded source data keyed by path (HDF5) or column name (CSV).
        """
        if self.source_props["type"] == "hdf5":
            # Flatten the nested HDF5 structure into absolute path → dataset dict.
            return HDF5Dict(self.source_props["path"]).flatten()

        # CSV: read selected columns and convert to numpy arrays.
        df = pd.read_csv(
            self.source_props["path"],
            sep=self.source_props["delimiter"],
            comment=self.source_props["comment"],
            usecols=self.source_config.get("cols", None),
        )
        return {col: df[col].to_numpy() for col in df.columns}
