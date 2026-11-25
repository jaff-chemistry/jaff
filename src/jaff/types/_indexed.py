# ABOUTME: IndexedValue/IndexedList containers for indexed array assignments used in
# ABOUTME: codegen, with normal/nested/flattened modes and conversions between them.
"""
Indexed expression containers for code generation.

This module provides data structures for representing indexed mathematical
expressions (array assignments) and collections thereof.  They are used
throughout the :mod:`jaff.codegen` package, particularly during common
subexpression elimination (CSE) and array-assignment code emission.

Classes
-------
IndexedValue
    Immutable two-element tuple of ``(indices, value)`` representing a
    single array assignment such as ``array[i][j] = expression``.
IndexedList
    Type-safe :class:`list` subclass that only accepts :class:`IndexedValue`
    elements.  Supports three structural modes (``normal``, ``nested``,
    ``flattened``) and provides conversion methods between them.

Structural modes
----------------
normal
    All :class:`IndexedValue` objects have single-element index lists, and no
    value contains :class:`IndexedValue` objects.  Values may be scalars or
    plain (possibly nested) lists, e.g. ``IndexedList([[1, 2], [3, 4]])``.
nested
    Each :class:`IndexedValue` value may itself be an :class:`IndexedList`
    of :class:`IndexedValue` objects, representing a 2-D structure.
flattened
    All :class:`IndexedValue` objects have multi-element index lists (e.g.
    ``[i, j]``) and scalar values.  This is the canonical form for
    generating nested array-access code such as ``k[i][j] = expr``.
"""

import logging
from functools import cached_property
from typing import Any, Iterable


class IndexedValue(tuple):
    """
    Immutable container for an indexed value with multi-dimensional indexing.

    Represents a value associated with one or more array indices, typically
    used for generating array assignment statements like:
        array[i][j] = expression

    This is an immutable tuple subclass with named properties for clarity.

    Attributes:
        indices: List of integer indices (e.g., [0], [2, 3], etc.)
        value: The value/expression associated with these indices

    Example:
        >>> iv = IndexedValue([0], "x + y")
        >>> iv.indices
        [0]
        >>> iv.value
        'x + y'
        >>> iv = IndexedValue([2, 3], "sin(theta)")  # For 2D array
    """

    def __new__(cls, indices: list[int], value: Any):
        """
        Create a new IndexedValue.

        Args:
            indices: List of integer indices for array access
            value: The value/expression to associate with these indices

        Raises:
            TypeError: If indices is not a list
        """
        if not isinstance(indices, list):
            raise TypeError(
                f"{indices}: indices must be of type list\n"
                f"Current index type detected: {type(indices)}"
            )
        return super().__new__(cls, (indices, value))

    @cached_property
    def indices(self) -> list[int]:
        """Get the list of indices."""
        return self[0]

    @cached_property
    def value(self) -> Any:
        """Get the associated value/expression."""
        return self[1]

    def __repr__(self):
        """Return detailed string representation of this IndexedValue.

        Returns
        -------
        str
            String including indices and value.
        """
        return f"IndexedValue(indices={self.indices!r}, value={self.value!r})"

    def __str__(self):
        """Return human-readable string representation showing the index mapping.

        Returns
        -------
        str
            String of form ``"[i] -> value"``.
        """
        if isinstance(self.value, list):
            value_str = "[" + ", ".join(str(v) for v in self.value) + "]"
        else:
            value_str = str(self.value)

        return f"{self.indices} -> {value_str}"


class IndexedList(list):
    """
    Type-safe list that only accepts IndexedValue elements.

    Provides a specialized list container ensuring all elements are
    IndexedValue instances, used to collect indexed expressions for
    code generation.

    Example:
        >>> items = IndexedList()
        >>> items.append(IndexedValue([0], "2*x"))
        >>> items.append(IndexedValue([1], "y**2"))
        >>> len(items)
        2
    """

    def __init__(
        self,
        items: Iterable | None = None,
        nested: bool = False,
        flatten: bool = False,
        verbose: bool = False,
    ):
        """
        Initialize an IndexedList from an iterable.

        Args:
            items: Optional iterable to convert to IndexedList. If None, creates empty list.
                Can contain IndexedValue objects or raw values to be wrapped.
            nested: If True, nested iterables will be preserved in IndexedValue.value
                as nested structures. Cannot be used with flatten=True.
            flatten: If True, nested iterables will be flattened with multi-dimensional
                indices (e.g., [0, 1] for nested element). Cannot be used with nested=True.
            verbose: If True, emit warnings for no-op conversions. Default: False.

        Raises:
            ValueError: If both nested=True and flatten=True
            TypeError: If items contain a mix of IndexedValue and non-IndexedValue objects

        Example:
            >>> print(IndexedList([1, 2, 3]))  # IndexedValue objects with single indices
            IndexedList[
            [0] -> 1
            [1] -> 2
            [2] -> 3
            ]
            >>> print(IndexedList([[1, 2], [3, 4]], flatten=True))  # 2D indices
            IndexedList[
            [0, 0] -> 1
            [0, 1] -> 2
            [1, 0] -> 3
            [1, 1] -> 4
            ]
            >>> print(IndexedList([[1, 2], [3, 4]], nested=True))  # nested structure
            IndexedList[
            [0] -> [[0] -> 1, [1] -> 2]
            [1] -> [[0] -> 3, [1] -> 4]
            ]
        """
        self.verbose = verbose
        if flatten and nested:
            raise ValueError("Cannot have both nested=True and flatten=True")

        if items is None:
            items = []

        # if not isinstance(items, list):
        items = list(items)

        if flatten:
            flat: list[IndexedValue] = []
            self.__convert_to_indexed_list(
                items,
                nested=nested,
                flatten=True,
                index_prefix=[],
                out=flat,
            )
            items = flat
        else:
            self.__convert_to_indexed_list(
                items,
                nested=nested,
                flatten=False,
                index_prefix=[],
            )

        super().__init__(items)

    @cached_property
    def _logger(self) -> logging.Logger:
        """Lazy-loaded JAFF logger for warning messages.

        Returns
        -------
        logging.Logger
        """
        from ..io._logger import JaffLogger

        return JaffLogger().get_logger()

    def __repr__(self):
        """Return detailed string representation of this IndexedList.

        Returns
        -------
        str
            Wraps the standard list repr in ``"IndexedList(...)"``
        """
        return f"IndexedList({list.__repr__(self)})"

    def __str__(self):
        """Return multi-line human-readable string of all indexed entries.

        Returns
        -------
        str
            Each entry on its own line, enclosed in ``"IndexedList[...]"``.
        """
        return "IndexedList[\n" + "\n".join(str(x) for x in self) + "\n]"

    def append(self, item: IndexedValue):
        """
        Append an IndexedValue to the list.

        Args:
            item: IndexedValue to append

        Raises:
            TypeError: If item is not an IndexedValue

        Warning:
            This method returns None (standard Python list behavior).
            Do NOT assign the return value:
                WRONG: out = out.append(item)  # out becomes None!
                RIGHT: out.append(item)         # out is modified in-place

        Example:
            >>> lst = IndexedList()
            >>> lst.append(IndexedValue([0], "value"))
            >>> len(lst)
            1
        """
        if not isinstance(item, IndexedValue):
            raise TypeError(f"{item} must be of type IndexedValue")

        super().append(item)

    def extend(self, items: Iterable[IndexedValue]):
        """
        Extend the list with multiple IndexedValue objects.

        Args:
            items: Iterable of IndexedValue objects

        Raises:
            TypeError: If any item is not an IndexedValue

        Warning:
            This method returns None (standard Python list behavior).
            Do NOT assign the return value:
                WRONG: out = out.extend(items)  # out becomes None!
                RIGHT: out.extend(items)         # out is modified in-place

        Example:
            >>> lst = IndexedList()
            >>> lst.extend([IndexedValue([0], "a"), IndexedValue([1], "b")])
            >>> len(lst)
            2
        """
        items = list(items)
        if any(not isinstance(item, IndexedValue) for item in items):
            raise TypeError(f"All items are not of type IndexedValue in: {items}")

        super().extend(items)

    @staticmethod
    def __is_iterable(obj: Any) -> bool:
        """
        Check if an object is iterable (but not a string).

        Args:
            obj: Object to check

        Returns:
            True if obj is iterable and not a string, False otherwise
        """
        return isinstance(obj, Iterable) and not isinstance(obj, (str, bytes))

    def __convert_to_indexed_list(
        self,
        items: list[Any],
        nested: bool,
        flatten: bool,
        index_prefix: list[int],
        out: list[IndexedValue] | None = None,
    ) -> None:
        """
        Recursively convert items to IndexedValue objects.

        This method handles conversion of raw values to IndexedValue objects,
        supporting both nested and flattened structures.

        Args:
            items: List of items to convert (modified in-place unless flatten=True)
            nested: If True, preserve nested iterables in IndexedValue.value
            flatten: If True, flatten nested structures with multi-dimensional indices
            index_prefix: Current index path for nested structures (used with flatten)
            out: Output list for flattened results (used only when flatten=True)
        """
        any_indexed = any(isinstance(x, IndexedValue) for x in items)
        all_indexed = all(isinstance(x, IndexedValue) for x in items)

        if any_indexed and not all_indexed:
            raise TypeError(
                "Items must be either ALL IndexedValue or NONE IndexedValue\n"
                f"{[(item, type(item)) for item in items]}"
            )

        if all_indexed:
            # Normalize IndexedValue objects to ensure nested iterables are IndexedList
            for i, item in enumerate(items):
                if (
                    isinstance(item, IndexedValue)
                    and self.__is_iterable(item.value)
                    and not isinstance(item.value, IndexedList)
                ):
                    # Materialize once: probing a generator would exhaust it.
                    payload = list(item.value)
                    if any(isinstance(v, IndexedValue) for v in payload):
                        items[i] = IndexedValue(item.indices, IndexedList(payload))
                    else:
                        # Retain the materialized payload, not the consumed original.
                        items[i] = IndexedValue(item.indices, payload)
            if flatten and out is not None:
                out.extend(items)
            return

        for i, item in enumerate(items):
            idx = index_prefix + [i] if flatten else [i]

            if self.__is_iterable(item):
                materialized = list(item)
                if flatten:
                    self.__convert_to_indexed_list(
                        materialized,
                        nested=nested,
                        flatten=True,
                        index_prefix=idx,
                        out=out,
                    )
                elif nested:
                    items[i] = IndexedValue(idx, IndexedList(materialized, nested=True))
                else:
                    items[i] = IndexedValue(idx, materialized)
            else:
                iv = IndexedValue(idx, item)
                if flatten and out is not None:
                    out.append(iv)
                else:
                    items[i] = iv

    def type(self) -> str:
        """
        Detect whether the IndexedList is 'normal', 'nested', or 'flattened'.

        Returns:
            str: One of 'normal', 'nested', or 'flattened'
                - 'flattened': Any IndexedValue has multi-dimensional indices (len > 1)
                - 'nested': Any IndexedValue has an iterable value containing IndexedValue objects
                - 'normal': All IndexedValues have single indices and no value contains IndexedValues

        Example:
            >>> lst = IndexedList([1, 2, 3])
            >>> lst.type()
            'normal'
            >>> lst = IndexedList([[1, 2], [3, 4]], nested=True)
            >>> lst.type()
            'nested'
            >>> lst = IndexedList([[1, 2], [3, 4]], flatten=True)
            >>> lst.type()
            'flattened'
        """
        if not self:
            return "normal"

        # Check for flattened: any multi-dimensional indices
        has_multidim_indices = any(len(item.indices) > 1 for item in self)
        if has_multidim_indices:
            return "flattened"

        # Check for nested: any iterable values containing IndexedValue objects
        for item in self:
            if self.__is_iterable(item.value):
                # Check if the iterable contains IndexedValue objects
                if any(isinstance(v, IndexedValue) for v in item.value):
                    return "nested"

        return "normal"

    def __normal_to_nested(self) -> "IndexedList":
        """Convert normal list to nested format."""
        # Normal with iterable values: convert to nested format
        # by wrapping list items in IndexedValue objects
        result = []
        for item in self:
            if self.__is_iterable(item.value):
                # Convert to list if needed
                nested_val = (
                    list(item.value) if not isinstance(item.value, list) else item.value
                )
                # Wrap each item in the list with IndexedValue
                nested_indexed_vals = []
                for i, val in enumerate(nested_val):
                    nested_indexed_vals.append(IndexedValue([i], val))
                # Wrap the list of IndexedValues in an IndexedList
                result.append(
                    IndexedValue(item.indices, IndexedList(nested_indexed_vals))
                )
            else:
                result.append(item)
        return IndexedList(result)

    def __to_flattened(self) -> "IndexedList":
        """Convert a normal or nested list to flattened format.

        Preserves every explicit :class:`IndexedValue` child index and appends
        every raw nesting dimension, so the flattened indices are the full
        multi-dimensional coordinate of each leaf value.
        """
        result: list[IndexedValue] = []

        def flatten_recursive(item: IndexedValue, prefix_indices: list[int]) -> None:
            """Recursively flatten an IndexedValue into ``result``."""
            if not self.__is_iterable(item.value):
                # Leaf value: emit at the accumulated coordinate.
                result.append(IndexedValue(prefix_indices, item.value))
                return

            for i, child in enumerate(item.value):
                if isinstance(child, IndexedValue):
                    # Explicit child index: keep it, do not renumber.
                    flatten_recursive(child, prefix_indices + child.indices)
                else:
                    # Raw child: its position is the next dimension.
                    flatten_recursive(IndexedValue([i], child), prefix_indices + [i])

        for item in self:
            flatten_recursive(item, item.indices)

        return IndexedList(result)

    def __nested_to_normal(self) -> "IndexedList":
        """Convert nested list to normal format by reconstructing list values from IndexedValues."""

        def rebuild(value: Any) -> Any:
            """Collapse nested IndexedValues back into plain (possibly nested) lists."""
            if self.__is_iterable(value) and any(
                isinstance(v, IndexedValue) for v in value
            ):
                return [
                    rebuild(v.value) if isinstance(v, IndexedValue) else rebuild(v)
                    for v in value
                ]
            return value

        # Preserve each item's own outer index; only reshape the value.
        return IndexedList(
            [IndexedValue(item.indices, rebuild(item.value)) for item in self]
        )

    @staticmethod
    def __reconstruct_dense(entries: list[tuple[list[int], Any]]) -> Any:
        """Rebuild a dense nested list from ``(remaining_indices, value)`` entries.

        Args:
            entries: Pairs of trailing index coordinates and their leaf values,
                all sharing the same coordinate length.

        Returns:
            The scalar value when coordinates are exhausted, otherwise a nested
            list indexed contiguously from 0.

        Raises:
            ValueError: If coordinates are ragged, collide, or are sparse
                (non-contiguous), i.e. cannot be represented as a dense list.
        """
        from collections import defaultdict

        if all(len(remaining) == 0 for remaining, _ in entries):
            if len(entries) != 1:
                raise ValueError(
                    f"Colliding flattened indices map to the same coordinate: {entries}"
                )
            return entries[0][1]

        groups: dict[int, list[tuple[list[int], Any]]] = defaultdict(list)
        for remaining, value in entries:
            if not remaining:
                raise ValueError(
                    "Ragged flattened indices: mixed coordinate depths under the "
                    f"same parent: {entries}"
                )
            groups[remaining[0]].append((remaining[1:], value))

        keys = sorted(groups.keys())
        if keys != list(range(len(keys))):
            raise ValueError(
                f"Sparse flattened indices cannot be reconstructed to a dense "
                f"nested list; got dimension keys {keys}"
            )
        return [IndexedList.__reconstruct_dense(groups[k]) for k in keys]

    def __flattened_to_normal(self) -> "IndexedList":
        """Convert flattened list to normal format by reconstructing nested list values.

        The leading index becomes the outer (normal) index and may be sparse;
        the trailing dimensions are reconstructed into dense nested lists.
        """
        from collections import defaultdict

        # Group by leading (outer) index; trailing dims rebuild the value.
        grouped: dict[int, list[tuple[list[int], Any]]] = defaultdict(list)
        for item in self:
            first_idx = item.indices[0]
            remaining_indices = item.indices[1:]
            grouped[first_idx].append((remaining_indices, item.value))

        result = []
        for idx in sorted(grouped.keys()):
            value = self.__reconstruct_dense(grouped[idx])
            result.append(IndexedValue([idx], value))

        return IndexedList(result)

    def __flattened_to_nested(self) -> "IndexedList":
        """Convert flattened list to nested format."""
        from collections import defaultdict

        grouped = defaultdict(list)

        for item in self:
            first_idx = item.indices[0]
            remaining_indices = item.indices[1:]
            grouped[first_idx].append((remaining_indices, item.value))

        result = []
        for idx in sorted(grouped.keys()):
            items_at_idx = grouped[idx]

            # If only one item with no remaining indices, keep as scalar value
            if len(items_at_idx) == 1 and not items_at_idx[0][0]:
                result.append(IndexedValue([idx], items_at_idx[0][1]))
            else:
                # Create nested structure with IndexedValues
                nested_values = []
                for remaining_idx, value in items_at_idx:
                    if remaining_idx:  # Still has indices, create IndexedValue
                        nested_values.append(IndexedValue(remaining_idx, value))
                    else:
                        nested_values.append(value)
                # Wrap the list of IndexedValues in an IndexedList
                result.append(IndexedValue([idx], IndexedList(nested_values)))

        return IndexedList(result)

    def nested(self) -> "IndexedList":
        """
        Convert the IndexedList to nested format.

        Returns:
            IndexedList: New IndexedList in nested format

        Raises:
            ValueError: If conversion to nested is not possible

        Example:
            >>> lst = IndexedList([[1, 2], [3, 4]], flatten=True)
            >>> nested = lst.nested()
            >>> nested.type()
            'nested'
        """
        list_type = self.type()

        # If already nested, return a copy
        if list_type == "nested":
            return IndexedList(list(self))

        # Convert from normal
        if list_type == "normal":
            if not any(self.__is_iterable(item.value) for item in self):
                if self.verbose:
                    self._logger.warning(
                        "Cannot convert to nested format: no iterable values. Returning copy."
                    )
                return IndexedList(list(self))
            return self.__normal_to_nested()

        # Convert from flattened
        if list_type == "flattened":
            return self.__flattened_to_nested()

        # Should never reach here
        raise ValueError(f"Unknown list type: {list_type}")

    def flatten(self) -> "IndexedList":
        """
        Convert the IndexedList to flattened format.

        Returns:
            IndexedList: New IndexedList in flattened format

        Raises:
            ValueError: If conversion to flattened is not possible

        Example:
            >>> lst = IndexedList([[1, 2], [3, 4]], nested=True)
            >>> flat = lst.flatten()
            >>> flat.type()
            'flattened'
        """
        list_type = self.type()

        # If already flattened, return a copy
        if list_type == "flattened":
            return IndexedList(list(self))

        # Convert from normal
        if list_type == "normal":
            if not any(self.__is_iterable(item.value) for item in self):
                if self.verbose:
                    self._logger.warning(
                        "Cannot convert to flattened format: no iterable values. Returning copy."
                    )
                return IndexedList(list(self))
            return self.__to_flattened()

        # Convert from nested
        if list_type == "nested":
            return self.__to_flattened()

        # Should never reach here
        raise ValueError(f"Unknown list type: {list_type}")

    def normal(self) -> "IndexedList":
        """
        Convert the IndexedList to normal format.

        Returns:
            IndexedList: New IndexedList in normal format

        Raises:
            ValueError: If conversion to normal is not possible

        Example:
            >>> lst = IndexedList([[1], [2]], nested=True)
            >>> normal = lst.normal()
            >>> normal.type()
            'normal'
        """
        list_type = self.type()

        # If already normal, return a copy
        if list_type == "normal":
            return IndexedList(list(self))

        # Convert from flattened
        if list_type == "flattened":
            return self.__flattened_to_normal()

        # Convert from nested
        if list_type == "nested":
            return self.__nested_to_normal()

        # Should never reach here
        raise ValueError(f"Unknown list type: {list_type}")
