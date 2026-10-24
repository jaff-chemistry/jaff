# ABOUTME: jaff.core.chemistry: chemistry reference data as module functions; the
# ABOUTME: element table is read once and used by species/elements with or without a network

import subprocess
import sys

import pytest

from jaff.core import Element, Elements, Specie, Species, chemistry


class TestElementTable:
    def test_holds_the_element_data(self):
        table = chemistry.element_table()

        assert table["H"]["mass"] == pytest.approx(1.673773e-24)
        assert table["C"]["protons"] == 6

    def test_is_read_once(self):
        assert chemistry.element_table() is chemistry.element_table()

    def test_is_the_network_table(self, make_network):
        net = make_network("H + H -> H2  []  1e-17\n", name="c.jet")

        assert net.mass_dict is chemistry.element_table()

    def test_old_loader_is_gone(self):
        import jaff.common

        assert not hasattr(jaff.common, "load_mass_dict")


class TestNoClassLevelTable:
    @pytest.mark.parametrize("cls", [Specie, Species, Element, Elements])
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
