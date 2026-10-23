# ABOUTME: Tests for the Elements flyweight and its composition matrices
# ABOUTME: Column order must follow the caller's species order, not a sorted key

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
