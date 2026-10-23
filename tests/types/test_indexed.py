# ABOUTME: Tests for IndexedList/IndexedValue structural conversions
# ABOUTME: Exact-index coverage for flatten <-> normal <-> nested round trips

import pytest

from jaff.types import IndexedList, IndexedValue


def indices_of(lst):
    """Extract the list of index-lists from an IndexedList."""
    return [iv.indices for iv in lst]


def pairs_of(lst):
    """Extract (indices, value) pairs from an IndexedList."""
    return [(iv.indices, iv.value) for iv in lst]


# --------------------------------------------------------------------------- #
# flatten: raw nested iterables                                               #
# --------------------------------------------------------------------------- #
def test_flatten_3d_raw_iterable_exact_indices():
    """Every raw nesting dimension must appear in the flattened index."""
    flat = IndexedList([[[1, 2], [3, 4]]]).flatten()
    assert pairs_of(flat) == [
        ([0, 0, 0], 1),
        ([0, 0, 1], 2),
        ([0, 1, 0], 3),
        ([0, 1, 1], 4),
    ]


def test_flatten_no_duplicate_index_tuples():
    """Flatten must never emit colliding index tuples (silent overwrite guard)."""
    flat = IndexedList([[[1, 2], [3, 4]]]).flatten()
    idxs = [tuple(i) for i in indices_of(flat)]
    assert len(idxs) == len(set(idxs))


def test_flatten_preserves_explicit_child_index():
    """Explicit IndexedValue child indices survive flatten (not renumbered)."""
    items = IndexedList(
        [IndexedValue([3], IndexedList([IndexedValue([7], 42)]))]
    )
    flat = items.flatten()
    assert pairs_of(flat) == [([3, 7], 42)]


def test_flatten_2d_raw_iterable_exact_indices():
    flat = IndexedList([[1, 2], [3, 4]], nested=True).flatten()
    assert pairs_of(flat) == [
        ([0, 0], 1),
        ([0, 1], 2),
        ([1, 0], 3),
        ([1, 1], 4),
    ]


# --------------------------------------------------------------------------- #
# nested -> normal: preserve outer index                                      #
# --------------------------------------------------------------------------- #
def test_nested_to_normal_preserves_outer_index():
    """Outer index must be kept, not overwritten with enumeration position."""
    items = IndexedList(
        [IndexedValue([3], IndexedList([IndexedValue([0], 1), IndexedValue([1], 2)]))]
    )
    normal = items.normal()
    assert pairs_of(normal) == [([3], [1, 2])]


# --------------------------------------------------------------------------- #
# flattened -> normal: reconstruct inner dimensions                           #
# --------------------------------------------------------------------------- #
def test_flattened_to_normal_reconstructs_nested_dims():
    """3-D flattened must round-trip back to nested lists, not a flat list."""
    n = IndexedList([[[1, 2], [3, 4]]])
    normal = n.flatten().normal()
    assert pairs_of(normal) == [([0], [[1, 2], [3, 4]])]


def test_flatten_normal_flatten_round_trip_idempotent():
    """normal -> flatten -> normal -> flatten reproduces the flattened form."""
    once = IndexedList([[[1, 2], [3, 4]]]).flatten()
    twice = once.normal().flatten()
    assert pairs_of(twice) == pairs_of(once)


def test_flattened_to_normal_rejects_sparse_indices():
    """Non-dense inner coordinates must raise, not silently renumber."""
    sparse = IndexedList(
        [IndexedValue([0, 0], 1), IndexedValue([0, 2], 3)]  # missing [0, 1]
    )
    with pytest.raises(ValueError):
        sparse.normal()


# --------------------------------------------------------------------------- #
# construction: iterable materialization (tuples, generators, no mutation)    #
# --------------------------------------------------------------------------- #
def test_nested_tuple_child():
    """Tuple children must not fail on item assignment."""
    lst = IndexedList([(1, 2)], nested=True)
    child = lst[0].value
    assert isinstance(child, IndexedList)
    assert pairs_of(child) == [([0], 1), ([1], 2)]


def test_nested_generator_child_not_exhausted():
    """Generator children must be materialized, not consumed to empty."""
    lst = IndexedList([(i for i in [1, 2])], nested=True)
    child = lst[0].value
    assert isinstance(child, IndexedList)
    assert pairs_of(child) == [([0], 1), ([1], 2)]


def test_nested_list_child():
    lst = IndexedList([[1, 2]], nested=True)
    assert pairs_of(lst[0].value) == [([0], 1), ([1], 2)]


def test_flatten_tuple_child():
    flat = IndexedList([(1, 2)], flatten=True)
    assert pairs_of(flat) == [([0, 0], 1), ([0, 1], 2)]


def test_flatten_generator_child_not_exhausted():
    flat = IndexedList([(i for i in [1, 2])], flatten=True)
    assert pairs_of(flat) == [([0, 0], 1), ([0, 1], 2)]


def test_nested_does_not_mutate_caller_list():
    """A caller's nested list must be left untouched by construction."""
    inner = [1, 2]
    outer = [inner]
    IndexedList(outer, nested=True)
    assert inner == [1, 2]
    assert outer == [[1, 2]]
    assert outer[0] is inner


def test_flatten_does_not_mutate_caller_list():
    inner = [1, 2]
    outer = [inner]
    IndexedList(outer, flatten=True)
    assert inner == [1, 2]
    assert outer == [[1, 2]]
