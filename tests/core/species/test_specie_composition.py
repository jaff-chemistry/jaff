# ABOUTME: Tests for Specie composition parsing: exploded atom list and mass from the name,
# ABOUTME: including explicit counts of one (C1) and multi-digit counts

import pytest

from jaff.core.species import Specie


@pytest.mark.parametrize(
    ("name", "atoms"),
    [
        ("C", ["C"]),
        ("C1", ["C"]),
        ("CH1", ["C", "H"]),
        ("C1H4", ["C", "H", "H", "H", "H"]),
        ("H2O", ["H", "H", "O"]),
        ("C12", ["C"] * 12),
        ("He", ["He"]),
    ],
)
def test_exploded_atoms(name, atoms):
    assert Specie(name).exploded == atoms


def test_an_explicit_count_of_one_keeps_the_mass():
    assert Specie("C1").mass == pytest.approx(Specie("C").mass, rel=1e-12, abs=0)
