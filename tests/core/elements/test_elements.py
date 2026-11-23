# ABOUTME: Tests for the Elements flyweight and its composition matrices
# ABOUTME: Column order must follow the caller's species order, not a sorted key

import pytest

from jaff import Elements


def test_flyweight_distinguishes_species_order():
    """Different species orders must be distinct flyweight instances."""
    first = Elements(["H", "H2"])
    second = Elements(["H2", "H"])
    assert first is not second


def test_density_matrix_follows_caller_order():
    """density_matrix columns follow the caller's species order."""
    assert Elements(["H", "H2"]).density_matrix() == [[1, 2]]
    assert Elements(["H2", "H"]).density_matrix() == [[2, 1]]


def test_truth_matrix_follows_caller_order():
    """truth_matrix columns follow the caller's species order."""
    # H present in both H and H2 -> both columns 1; verify O ordering instead.
    assert Elements(["H2O", "H"]).truth_matrix() == [[1, 1], [1, 0]]
    assert Elements(["H", "H2O"]).truth_matrix() == [[1, 1], [0, 1]]


def test_same_order_is_shared_flyweight():
    """Identical species order still reuses the cached instance."""
    assert Elements(["H", "H2"]) is Elements(["H", "H2"])


def test_a_species_catalogue_is_accepted():
    from jaff.core.species import Species

    species = Species(["CO", "H2O"])

    assert Elements(species).symbols() == ["C", "H", "O"]
    assert Elements(species) is Elements(list(species))


def test_a_tuple_is_accepted():
    assert Elements(("CO", "CH4")).symbols() == ["C", "H", "O"]


def test_a_mixed_list_of_names_and_species_is_accepted():
    from jaff.core.species import Specie

    assert Elements(["H2", Specie("CO")]).symbols() == ["C", "H", "O"]


def test_a_single_species_or_name_is_accepted():
    from jaff.core.species import Specie

    assert Elements("CO").symbols() == ["C", "O"]
    assert Elements(Specie("NH3")).symbols() == ["H", "N"]


@pytest.mark.parametrize("method", ["truth_matrix", "density_matrix"])
def test_matrices_are_fresh_copies(method):
    """Editing a returned matrix must not change what the next call returns."""
    elements = Elements(["CO", "H2O", "CH4"])
    first = getattr(elements, method)()
    expected = [row[:] for row in first]
    first[0][0] = 99
    first.append([7])

    assert getattr(elements, method)() == expected
