# ABOUTME: jaff.chemistry.periodic_table: the bundled element table read once, as Element
# ABOUTME: objects; species and elements use it with or without a network

import subprocess
import sys

import pytest

from jaff.chemistry import periodic_table
from jaff.core import Element, Elements, Specie, Species


class TestPeriodicTable:
    def test_whole_table_holds_elements(self):
        table = periodic_table()

        assert isinstance(table["H"], Element)
        assert table["C"].protons == 6

    def test_is_read_only(self):
        with pytest.raises(TypeError):
            periodic_table()["Xx"] = periodic_table("H")  # type: ignore[index]

    def test_is_read_once(self):
        assert periodic_table() is periodic_table()

    def test_one_element_by_symbol(self):
        hydrogen = periodic_table("H")

        assert hydrogen is periodic_table()["H"]
        assert hydrogen.name == "Hydrogen"
        assert hydrogen.mass == pytest.approx(1.673773e-24)

    def test_unknown_symbol_is_named(self):
        with pytest.raises(KeyError, match="Xx"):
            periodic_table("Xx")

    @pytest.mark.parametrize(
        "symbol, mass",
        [
            ("H", 1.673773e-24),
            ("He", 6.646473e-24),
            ("C", 1.994473e-23),
            ("O", 2.656763e-23),
        ],
    )
    def test_masses(self, symbol, mass):
        assert periodic_table(symbol).mass == pytest.approx(mass, rel=1e-6, abs=0)


class TestElement:
    def test_constructor_is_a_table_lookup(self):
        assert Element("H") is periodic_table("H")

    def test_unknown_symbol_is_named(self):
        with pytest.raises(KeyError, match="Xx"):
            Element("Xx")

    def test_compares_by_symbol(self):
        assert sorted([Element("O"), Element("C")]) == [Element("C"), Element("O")]
        assert str(Element("He")) == "He"


class TestHome:
    def test_element_lives_in_core(self):
        import jaff
        import jaff.chemistry
        import jaff.core.elements

        assert jaff.Element is jaff.core.elements.Element is Element
        assert not hasattr(jaff.chemistry, "Element")

    def test_old_names_are_gone(self):
        import jaff.core

        for name in ("ElementProps", "chemistry"):
            assert not hasattr(jaff.core, name)

    def test_network_does_not_carry_the_table(self, make_network):
        net = make_network("H + H -> H2  []  1e-17\n", name="c.jet")

        assert not hasattr(net, "mass_dict")

    @pytest.mark.parametrize("cls", [Specie, Species, Elements])
    def test_no_configure_chain(self, cls):
        assert not hasattr(cls, "configure")
        assert not hasattr(cls, "_mass_dict")


def test_species_and_elements_work_without_a_network():
    code = (
        "from jaff.core import Element, Elements, Specie\n"
        "co = Specie('CO')\n"
        "print(co.mass, Element('C').name, Elements(['CO']).count)\n"
    )
    out = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    ).stdout.split()

    assert float(out[0]) == pytest.approx(4.651236e-23, rel=1e-6)
    assert out[1:] == ["Carbon", "2"]
