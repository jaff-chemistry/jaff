# ABOUTME: Species, Elements and Reactions catalogues repr as their own class name
# ABOUTME: (not the shared base name Catalogue) followed by their items

import pytest

_LINES = ["H + H -> H2 [10,1000] 1e-10", "H2 + O -> H2O [10,1000] 1e-10"]


@pytest.fixture
def net(make_network):
    return make_network(_LINES)


@pytest.mark.parametrize(
    ("attr", "expected"),
    [
        ("species", "Species(['H', 'H2', 'O', 'H2O'])"),
        ("elements", "Elements(['H', 'O'])"),
        ("reactions", "Reactions(['H + H -> H2', 'H2 + O -> H2O'])"),
    ],
)
def test_repr_names_the_catalogue(net, attr, expected):
    assert repr(getattr(net, attr)) == expected


def test_a_species_elements_repr(net):
    assert repr(net.species["H2O"].elements) == "Elements(['H', 'O'])"
