# ABOUTME: HDF5Dict: dict mirror of an HDF5 file's group/dataset hierarchy with
# ABOUTME: leaf-dataset schema (_kind/_data/_dtype/_attrs), filters and flatten/nest.
"""
HDF5 in-memory dictionary representation.

This module provides :class:`HDF5Dict`, a :class:`dict` subclass that can be
constructed directly from an HDF5 file (or an open :mod:`h5py` file handle)
and mirrors its group/dataset hierarchy as nested Python dictionaries.

Schema
------
Each HDF5 dataset is represented by a **leaf dictionary** with the following
keys:

``_kind``
    ``"linear"`` for a plain N-D array; ``"compound"`` for a NumPy
    structured-dtype (record) array.
``_data``
    The actual data as a NumPy array (result of ``dataset[()]``).  Text
    (variable- or fixed-length strings, and the text fields of a compound
    dataset) is decoded to ``str`` in an ``object`` array, not ``bytes``.
``_dtype``
    For ``"linear"``: a JAFF dtype token string (e.g. ``"f64"``).
    For ``"compound"``: a ``{field_name: dtype_token}`` mapping.
``_attrs``
    ``{attr_name: attr_value}`` dictionary of HDF5 dataset attributes.
    Always present when read from a file; an empty dict if the dataset has
    no attributes.
``_name``
    Optional human-readable name, read from the ``_name`` HDF5 attribute
    if it was previously stored by :class:`~jaff.drivers.HDF5`.

HDF5 groups are represented as plain nested dictionaries.  A group's HDF5
attributes (if any) are stored under its ``_attrs`` key.  File-level
attributes are stored under the top-level ``_attrs`` key.

Dtype token mapping
-------------------
:meth:`HDF5Dict._to_np` maps JAFF token strings (``"f64"``, ``"i32"``, etc.)
to NumPy dtypes.  :meth:`HDF5Dict._from_np` provides the inverse mapping.

Examples
--------
>>> import h5py
>>> import numpy as np
>>> with h5py.File("rates.hdf5", "w") as f:
...     f["temperatures"] = np.logspace(1, 4, 5)
...     f["rates/k1"] = np.ones(5)
>>> hd = HDF5Dict("rates.hdf5")
>>> hd["temperatures"]["_kind"]
'linear'
>>> hd["temperatures"]["_dtype"]
'f64'
>>> hd["rates"]["k1"]["_attrs"]
{}
"""

import fnmatch
from contextlib import nullcontext
from pathlib import Path
from typing import Any

import h5py
import numpy as np


class HDF5Dict(dict):
    """
    Dictionary representation of an HDF5 file hierarchy.

    Can be constructed from an HDF5 file path, an open :class:`h5py.File`
    handle, or a plain Python dictionary (which is stored as-is without
    any HDF5 parsing).

    The nested structure mirrors the HDF5 group hierarchy.  Leaf nodes
    (datasets) follow the schema described in the module docstring.

    Parameters
    ----------
    h5obj : h5py.File, h5py.Group, Path, str, or dict
        Source to build from:

        * ``h5py.File`` — read from an already-open file handle.
        * ``h5py.Group`` — parse the hierarchy from that group onward.
        * ``Path`` or ``str`` — open and parse the file at the given path.
          A ``"file.h5::/internal/group"`` string (``::`` delimiter) parses
          from the named internal group onward.
        * ``dict`` — use directly without any HDF5 parsing.
    include : str, list of str, or None, optional
        Keep only datasets whose bare name, or the bare name of any ancestor
        group, matches one of these :func:`fnmatch.fnmatch` patterns
        (default ``None``: keep all).  Ignored for a plain ``dict``.
    exclude : str, list of str, or None, optional
        Drop any dataset or group whose bare name, or any ancestor's, matches
        one of these patterns; wins over *include* (default ``None``: drop
        nothing).  Ignored for a plain ``dict``.

    Raises
    ------
    FileNotFoundError
        If *h5obj* is a path that does not exist on the filesystem.

    Examples
    --------
    >>> hd = HDF5Dict(
    ...     {"rates": {"_kind": "linear", "_data": arr, "_dtype": "f64", "_attrs": {}}}
    ... )
    >>> hd2 = HDF5Dict(Path("output.hdf5"))
    """

    def __init__(
        self,
        h5obj: "h5py.File | h5py.Group | Path | str | dict",
        *,
        include: "str | list[str] | None" = None,
        exclude: "str | list[str] | None" = None,
    ):
        """Build the HDF5Dict from an HDF5 file path, open file/group handle, or plain dict.

        Parameters
        ----------
        h5obj : h5py.File, h5py.Group, Path, str, or dict
            Source to build from.  An ``h5py.Group`` parses the hierarchy from
            that group onward.  A ``"file.h5::/internal/group"`` string does
            the same via the ``::`` delimiter.  Plain dicts are stored directly
            without any HDF5 parsing.
        include : str, list of str, or None, optional
            Keep only datasets whose name — or the name of any of their parent
            groups — matches at least one of these patterns.  Matching is by
            **bare name** (the last path component of each ancestor), using
            :func:`fnmatch.fnmatch` so shell-style wildcards (``*``, ``?``,
            ``[seq]``) and exact names both work.  ``None`` (default) keeps all
            datasets.  Datasets excluded by this filter are never read into
            memory.  Ignored when *h5obj* is a plain ``dict``.
        exclude : str, list of str, or None, optional
            Drop any object — dataset or group — whose name, or any ancestor
            group name, matches one of these patterns.  Excluding a group name
            drops its entire subtree.  ``exclude`` takes precedence over
            ``include``.  Same bare-name :func:`fnmatch.fnmatch` semantics.
            ``None`` (default) drops nothing.  Ignored when *h5obj* is a plain
            ``dict``.

        Raises
        ------
        FileNotFoundError
            If *h5obj* is a path that does not exist on the filesystem.
        """
        # Fast path: if passed a plain dict, wrap it directly.
        if isinstance(h5obj, dict):
            super().__init__(h5obj)
            return

        # A "file::/internal/group" string selects a sub-group: parse the
        # hierarchy from that group onward. The "::" delimiter avoids the
        # ambiguity between a filesystem path and an internal HDF5 path.
        group_path = None
        if isinstance(h5obj, str) and "::" in h5obj:
            h5obj, group_path = h5obj.split("::", 1)

        if isinstance(h5obj, (str, Path)):
            h5path = Path(h5obj)
            if not h5path.exists():
                raise FileNotFoundError(h5path)

        cm = (
            nullcontext(h5obj) if isinstance(h5obj, h5py.Group) else h5py.File(h5obj, "r")
        )
        result = {}

        with cm as f:
            # Descend to the requested sub-group, if one was given.
            target = f[group_path] if group_path is not None else f
            # Store the target's attributes (if any) at the top level.
            if target.attrs:
                result["_attrs"] = dict(target.attrs)
            # Walk every group and dataset, building the nested dict.
            target.visititems(
                lambda name, obj: self.__build_nested_dict(
                    result, name, obj, include, exclude
                )
            )

        if include is not None or exclude is not None:
            self.__prune_empty(result)

        super().__init__(result)

    # ------------------------------------------------------------------
    # Flatten / nest conversions
    # ------------------------------------------------------------------

    def flatten(self, d: dict = {}, parent_key: str = "") -> dict:
        """
        Convert a nested HDF5Dict to a flat ``{absolute_path: leaf_dict}`` mapping.

        Leaf dictionaries (those whose keys start with ``"_"``) are collected
        under their slash-separated absolute HDF5 path.  Intermediate groups
        are traversed recursively.  The root group is represented by the key
        ``"/"`` when it contains endpoint data.

        Parameters
        ----------
        d : dict, optional
            Sub-dictionary to flatten.  Defaults to ``self``.
        parent_key : str, optional
            Current path prefix used during recursion.

        Returns
        -------
        dict
            Flat ``{path: leaf_dict}`` mapping where every key is an absolute
            HDF5 path string.

        Examples
        --------
        >>> hd = HDF5Dict("rates.hdf5")  # the file from the module example
        >>> list(hd.flatten().keys())
        ['/rates/k1', '/temperatures']
        """
        if not d:
            d = self

        items = {}
        endpoint_data = {}
        has_endpoint = False

        # Collect all keys that start with "_" — these form the leaf payload.
        for k, v in d.items():
            if k.startswith("_"):
                endpoint_data[k] = v
                has_endpoint = True

        # Store the leaf payload under its absolute path (root uses "/").
        if has_endpoint:
            items[parent_key or "/"] = endpoint_data

        # Recurse into non-metadata child dictionaries.
        for k, v in d.items():
            if not k.startswith("_") and isinstance(v, dict):
                new_key = f"{parent_key}/{k}" if parent_key != "/" else f"/{k}"
                new_key = new_key.replace("//", "/")
                items.update(self.flatten(v, new_key))

        return items

    def nested(self, d: dict = {}) -> dict:
        """
        Convert a flat ``{absolute_path: leaf_dict}`` mapping back to a nested dict.

        This is the inverse of :meth:`flatten`.  The root key ``"/"`` populates
        the top level directly; all other keys are split on ``"/"`` to
        reconstruct the group hierarchy.

        Parameters
        ----------
        d : dict, optional
            Flat dictionary to convert.  Defaults to ``self``.

        Returns
        -------
        dict
            Nested dictionary mirroring the original HDF5 group structure.
        """
        if not d:
            d = self

        nested_dict = {}
        for path, data in d.items():
            if path == "/":
                # Root-level attributes / data go directly into the top dict.
                nested_dict.update(data)
                continue

            # Split the path and navigate/create intermediate group dicts.
            parts = [p for p in path.split("/") if p]
            current = nested_dict
            for part in parts[:-1]:
                current = current.setdefault(part, {})

            if parts:
                leaf = parts[-1]
                if leaf not in current:
                    current[leaf] = {}
                current[leaf].update(data)

        return nested_dict

    # ------------------------------------------------------------------
    # HDF5 parsing internals
    # ------------------------------------------------------------------

    def __build_nested_dict(
        self,
        root: dict,
        name: str,
        obj: "h5py.Group | h5py.Dataset",
        include: "str | list[str] | None" = None,
        exclude: "str | list[str] | None" = None,
    ) -> None:
        """
        Callback for :meth:`h5py.File.visititems` — populate *root* in-place.

        Navigates to the correct nesting level using the slash-separated
        *name*, then either inserts a dataset leaf dict or updates a group
        node with its HDF5 attributes.

        Parameters
        ----------
        root : dict
            Top-level result dictionary being built.
        name : str
            Slash-separated path of the HDF5 object (as provided by
            :meth:`h5py.File.visititems`).
        obj : h5py.Group or h5py.Dataset
            The HDF5 object at *name*.
        include : str, list of str, or None, optional
            Dataset-only allowlist of bare-name patterns; see :meth:`__init__`.
            Groups are never dropped here (dead branches are pruned afterwards),
            so a kept dataset's ancestor groups retain their attributes.
        exclude : str, list of str, or None, optional
            Denylist of bare-name patterns applied to both datasets and groups;
            see :meth:`__init__`.  Takes precedence over *include*.

        Returns
        -------
        None
        """
        parts = name.split("/")

        if exclude is not None and self.__name_match(parts, exclude):
            return

        # Navigate to the parent node, creating intermediate group dicts.
        current = root
        for part in parts[:-1]:
            current = current.setdefault(part, {})

        leaf_key = parts[-1]

        if isinstance(obj, h5py.Dataset):
            # Apply the dataset allowlist *before* touching the data, so
            # filtered-out datasets are never read into memory.
            if include is not None and not self.__name_match(parts, include):
                return
            # Read all dataset attributes; extract the optional _name tag.
            attrs = dict(obj.attrs)
            ds_name = attrs.pop("_name", None)
            node = {
                "_kind": "compound" if obj.dtype.fields else "linear",
                "_data": self.__read_data(obj),  # Read all data into memory
                "_attrs": attrs,
            }
            if ds_name is not None:
                node["_name"] = ds_name
            # Encode the dtype as a JAFF token or raw dtype string.
            if obj.dtype.fields:
                # Structured dtype: encode each field separately.
                node["_dtype"] = {
                    n: self.__decode_dtype(obj.dtype[n]) for n in obj.dtype.names
                }
            else:
                node["_dtype"] = self.__decode_dtype(obj.dtype)

            current[leaf_key] = node
        else:
            # Group node: ensure the dict key exists and attach attributes.
            group_node = current.setdefault(leaf_key, {})
            if obj.attrs:
                group_node["_attrs"] = dict(obj.attrs)

    @staticmethod
    def __name_match(parts: list[str], patterns: "str | list[str]") -> bool:
        """
        Return ``True`` if any path component matches any of *patterns*.

        Matching is by bare name using :func:`fnmatch.fnmatch`, so both exact
        names and shell-style wildcards (``*``, ``?``, ``[seq]``) are accepted.
        Testing every component means a match on an ancestor group name applies
        to its whole subtree.

        Parameters
        ----------
        parts : list of str
            Path components of the HDF5 object (``name.split("/")``).
        patterns : str or list of str
            A single pattern or a list of patterns to test against.

        Returns
        -------
        bool
            ``True`` if at least one component matches at least one pattern.
        """
        pats = [patterns] if isinstance(patterns, str) else patterns
        return any(fnmatch.fnmatch(part, pat) for part in parts for pat in pats)

    @classmethod
    def __prune_empty(cls, node: dict) -> None:
        """
        Recursively drop group sub-dicts that contain no surviving dataset.

        Walks *node* depth-first and removes any child group whose subtree holds
        no dataset leaf (no ``_kind`` anywhere).  Dataset leaves and metadata
        keys (those starting with ``"_"``) are left untouched.  Used after
        ``include``/``exclude`` filtering to discard groups emptied by the
        filter while keeping the ancestor groups of any retained dataset (and
        therefore their attributes).

        Parameters
        ----------
        node : dict
            Group dictionary to prune in-place.

        Returns
        -------
        None
        """
        for key in list(node):
            if key.startswith("_"):
                continue
            child = node[key]
            # Dataset leaves carry "_kind" and are always kept.
            if not isinstance(child, dict) or "_kind" in child:
                continue
            # Group node: prune its descendants first, then drop it if nothing
            # but metadata keys remain (i.e. it lost all its datasets/subgroups).
            cls.__prune_empty(child)
            if all(k.startswith("_") for k in child):
                del node[key]

    @staticmethod
    def __read_data(obj: h5py.Dataset) -> Any:
        """
        Read a dataset's values, decoding text to ``str``.

        h5py returns strings (variable- or fixed-length) as ``bytes``; a text
        dataset is read through :meth:`h5py.Dataset.asstr` instead, and the
        text fields of a compound dataset are decoded into ``object`` fields.

        Parameters
        ----------
        obj : h5py.Dataset
            Dataset to read.

        Returns
        -------
        Any
            The dataset's values (a NumPy array, or a scalar for a scalar
            dataset) with text as ``str``.
        """
        if h5py.check_string_dtype(obj.dtype) is not None:
            return obj.asstr()[()]

        data = obj[()]
        names = obj.dtype.names or ()
        text = {
            n: info.encoding
            for n in names
            if (info := h5py.check_string_dtype(obj.dtype[n])) is not None
        }
        if not text:
            return data

        decoded = np.empty(
            np.shape(data),
            dtype=[(n, object if n in text else obj.dtype[n]) for n in names],
        )
        for n in names:
            values = data[n]
            if n in text:
                values = np.vectorize(
                    lambda v, enc=text[n]: v.decode(enc) if isinstance(v, bytes) else v,
                    otypes=[object],
                )(values)
            decoded[n] = values

        return decoded

    def __decode_dtype(self, dtype) -> str:
        """
        Convert a NumPy dtype to a JAFF dtype token string.

        Uses the reverse mapping from :meth:`_from_np`.  Falls back to
        ``"s"`` for string/bytes/object dtypes, or the raw NumPy dtype string
        for anything else.

        Parameters
        ----------
        dtype : numpy.dtype
            NumPy dtype to encode.

        Returns
        -------
        str
            JAFF dtype token (e.g. ``"f64"``) or a raw dtype string.
        """
        dt = np.dtype(dtype)
        if dt in self._from_np():
            return self._from_np()[dt]
        # String, bytes, and object dtypes all map to the generic "s" token.
        return "s" if dt.kind in ("S", "U", "O") else str(dt)

    # ------------------------------------------------------------------
    # Dtype mapping tables
    # ------------------------------------------------------------------

    @staticmethod
    def _to_np() -> dict[str, Any]:
        """
        Return the mapping from JAFF dtype token strings to NumPy dtypes.

        Includes ``"f128"`` only on platforms where :attr:`numpy.float128`
        is available (typically Linux/macOS x86-64).

        Returns
        -------
        dict[str, Any]
            Keys are JAFF token strings; values are NumPy dtype objects or
            type constructors.

        Notes
        -----
        Supported tokens: ``"i8"``, ``"i16"``, ``"i32"``, ``"i64"``,
        ``"u8"``, ``"u16"``, ``"u32"``, ``"u64"``, ``"f16"``, ``"f32"``,
        ``"f64"``, ``"f128"`` (platform-dependent), ``"c64"``, ``"c128"``,
        ``"b"`` (bool), ``"s"`` (variable-length UTF-8 string).
        """
        mapping = {
            "i8": np.int8,
            "i16": np.int16,
            "i32": np.int32,
            "i64": np.int64,
            "u8": np.uint8,
            "u16": np.uint16,
            "u32": np.uint32,
            "u64": np.uint64,
            "f16": np.float16,
            "f32": np.float32,
            "f64": np.float64,
            "c64": np.complex64,
            "c128": np.complex128,
            "b": np.bool_,
            "s": h5py.string_dtype(),
        }

        # f128 is only available on some platforms (e.g. Linux x86-64).
        if hasattr(np, "float128"):
            mapping["f128"] = np.float128

        return mapping

    @staticmethod
    def _from_np() -> dict:
        """
        Return the inverse mapping from NumPy dtypes to JAFF token strings.

        Derives the mapping by inverting :meth:`_to_np`, omitting the ``"s"``
        (string) entry because :class:`h5py.string_dtype` is not a regular
        NumPy dtype and cannot be used as a dict key.

        Returns
        -------
        dict
            Keys are :class:`numpy.dtype` instances; values are JAFF token
            strings (e.g. ``numpy.dtype("float64")`` → ``"f64"``).
        """
        return {np.dtype(v): k for k, v in HDF5Dict._to_np().items() if k != "s"}
