# ABOUTME: Catalogue: generic list + name-keyed collection base for Species, Reactions and
# ABOUTME: Elements, with O(1) lookup by index, slice, name or serialized alias.
"""
Generic indexed catalogue collection.

This module provides :class:`Catalogue`, the abstract base for all JAFF
domain-object collections — :class:`~jaff.Species`, :class:`~jaff.Reactions`,
and :class:`~jaff.Elements`.

A :class:`Catalogue` stores items in both a positional list and a name-keyed
dictionary so that lookup by index, slice, or name all run in amortised O(1)
time.  A secondary ``_by_serialized`` dict allows alternative (e.g.
serialized/canonical) name aliases to be registered post-construction.
"""

from typing import Any, Generic, Iterator, SupportsIndex, TypeVar, overload

T = TypeVar("T")


class Catalogue(Generic[T]):
    """
    Generic indexed collection supporting lookup by name, index, or slice.

    Serves as the base class for domain collections such as
    :class:`~jaff.Species`, :class:`~jaff.Reactions`, and
    :class:`~jaff.Elements`.  Items are stored in parallel in a list
    (for positional access) and a dictionary (for name-based lookup).

    A second dictionary, ``_by_serialized``, can be populated manually
    after construction to support alternative lookup keys (e.g. serialized
    or canonicalized species names) without duplicating the item.

    Parameters
    ----------
    items : list[T] or None, optional
        Ordered list of items.  Must be provided together with *items_dict*.
    items_dict : dict[Any, T | list[T]] or None, optional
        Name-keyed mapping of the same items.  Must be provided together
        with *items*.  Values may be a single item or a list of items when a
        subclass allows several items to share one key.
    check_length : bool, optional
        If ``True`` (default), raise :exc:`ValueError` when the lengths of
        *items* and *items_dict* differ.

    Attributes
    ----------
    _list : list[T]
        Positionally ordered list of all items.
    _by_prop : dict[Any, T | list[T]]
        Primary name-to-item mapping (values may be lists in subclasses that
        permit key collisions).
    _by_serialized : dict[str, T]
        Secondary alias mapping (e.g. serialized names).  Populated
        externally by subclasses or application code.
    count : int
        Total number of items in the collection.

    Raises
    ------
    ValueError
        If exactly one of *items* / *items_dict* is supplied (both must be
        given or both omitted), or if *check_length* is ``True`` and the
        two sequences have different lengths.

    Examples
    --------
    >>> cat = Catalogue(items=["H", "He"], items_dict={"H": "H", "He": "He"})
    >>> cat[0]
    'H'
    >>> cat["He"]
    'He'
    >>> cat.count
    2
    """

    def __init__(
        self,
        items: list[T] | None = None,
        items_dict: dict[Any, T | list[T]] | None = None,
        check_length: bool = True,
    ):
        """Initialise the catalogue from a parallel list and dictionary.

        Parameters
        ----------
        items : list[T] or None, optional
            Ordered list of items.  Must be provided together with *items_dict*.
        items_dict : dict[Any, T | list[T]] or None, optional
            Name-keyed mapping of the same items (values may be lists in
            subclasses that permit key collisions).
        check_length : bool, optional
            If ``True`` (default), raise :exc:`ValueError` when the lengths of
            *items* and *items_dict* differ.

        Raises
        ------
        ValueError
            If exactly one of *items*/*items_dict* is supplied, or if
            *check_length* is ``True`` and their lengths differ.
        """
        # Enforce that both or neither sequence is provided.
        if (items is None) ^ (items_dict is None):
            raise ValueError("Both an list and a dict must be supplied")

        if items is not None and items_dict is not None and check_length:
            if len(items) != len(items_dict.keys()):
                raise ValueError("Length of both list and dict must be same")

        self._by_prop: dict[Any, T] = {} if items_dict is None else items_dict
        self._list: list[T] = [] if items is None else items
        # Alternative keys (e.g. serialized names) can be added post-init.
        self._by_serialized: dict[str, T] = {}
        self.count: int = len(self._list)

    # ------------------------------------------------------------------
    # Item access
    # ------------------------------------------------------------------

    @overload
    def __getitem__(self, key: str) -> T: ...

    @overload
    def __getitem__(self, key: slice) -> list[T]: ...

    @overload
    def __getitem__(self, key: SupportsIndex) -> T: ...

    def __getitem__(self, key: str | SupportsIndex | slice) -> T | list[T]:
        """
        Retrieve an item by name, integer index, or slice.

        Lookup priority for string keys: ``_by_prop`` first, then
        ``_by_serialized``.

        Parameters
        ----------
        key : str, int-like, or slice
            * ``str`` — look up by primary or serialized name.
            * ``int``/``SupportsIndex`` — positional lookup in ``_list``.
            * ``slice`` — return a sub-list from ``_list``.

        Returns
        -------
        T or list[T]
            The matching item, or a list of items for slice access.

        Raises
        ------
        KeyError
            If *key* is a string and is not found in either name dict.
        IndexError
            If *key* is an integer index out of range.
        """
        if isinstance(key, str):
            if key not in self._by_prop and key not in self._by_serialized:
                raise KeyError(f"'{key}' not found in catalogue")

            if key in self._by_prop:
                return self._by_prop[key]

            if key in self._by_serialized:
                return self._by_serialized[key]

            raise KeyError(f"{key!r} not found in catalogue")

        elif isinstance(key, slice):
            return self._list[key]

        # Integer / SupportsIndex — delegate to the list.
        return self._list[int(key)]

    # ------------------------------------------------------------------
    # Iteration and membership
    # ------------------------------------------------------------------

    def __iter__(self) -> Iterator[T]:
        """
        Iterate over items in insertion order.

        Returns
        -------
        Iterator[T]
            Iterator over the positional list.
        """
        return iter(self._list)

    def __len__(self) -> int:
        """
        Return the number of items in the catalogue.

        Returns
        -------
        int
            Same value as :attr:`count`.
        """
        return len(self._list)

    def __contains__(self, item) -> bool:
        """
        Test membership by item value or string name.

        Parameters
        ----------
        item : T or str
            If a string, checks both ``_by_prop`` and ``_by_serialized``.
            Otherwise, checks ``_list``.

        Returns
        -------
        bool
            ``True`` if *item* is present in the catalogue.
        """
        if isinstance(item, str):
            return item in self._by_prop or item in self._by_serialized

        return item in self._list
