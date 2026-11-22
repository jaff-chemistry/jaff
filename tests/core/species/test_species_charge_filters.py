# ABOUTME: Species.neutral / Species.charged: the species themselves, or one attribute of
# ABOUTME: each (names, masses, ...); ne=True drops the electron; bad attribute names raise

import pytest

from jaff.core.species import Specie, Species


@pytest.fixture
def species():
    return Species(["H2", "H+", "e-", "O"])


def test_without_attr_they_are_species(species):
    assert all(isinstance(s, Specie) for s in species.charged())
    assert [s.name for s in species.neutral()] == ["H2", "O"]
    assert [s.name for s in species.charged()] == ["H+", "e-"]


def test_with_attr_they_are_the_values(species):
    assert species.charged("name") == ["H+", "e-"]
    assert species.charged("charge") == [1, -1]
    assert species.neutral("name") == ["H2", "O"]
    assert all(isinstance(m, float) for m in species.neutral("mass"))


def test_ne_drops_the_electron(species):
    assert species.charged("name", ne=True) == ["H+"]
    assert [s.name for s in species.charged(ne=True)] == ["H+"]


@pytest.mark.parametrize("method", ["neutral", "charged"])
def test_unknown_attr_raises(species, method):
    with pytest.raises(ValueError, match="Invalid attribute"):
        getattr(species, method)("nope")
