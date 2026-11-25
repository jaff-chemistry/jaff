# ABOUTME: Tests for rate tables written by Network.to_hdf5 / to_txt / write_table: the
# ABOUTME: network label is recorded in both formats; the quokka reaction_coeff layout stays

import h5py
import pytest

from jaff.drivers import HDF5

_LINES = ["H + H -> H2 [10,1000] 1e-10*(tgas/300)**0.5", "H2 -> H + H [10,1000] 1"]


@pytest.fixture
def net(make_network):
    return make_network(_LINES, name="tiny.dat")


def _root_attrs(path):
    with h5py.File(path) as f:
        return dict(f.attrs)


class TestHdf5Label:
    def test_explicit_label_is_a_root_attribute(self, net, tmp_path):
        out = tmp_path / "rates.hdf5"
        net.to_hdf5(out, label="my network", err_tol=None)

        attrs = _root_attrs(out)
        assert attrs["network_name"] == "my network"
        assert attrs["description"] == "JAFF auto-generated rate coefficient table"

    def test_label_defaults_to_the_network_label(self, net, tmp_path):
        out = tmp_path / "rates.hdf5"
        net.to_hdf5(out, err_tol=None)

        assert _root_attrs(out)["network_name"] == net.label

    def test_write_table_records_the_label(self, net, tmp_path):
        out = tmp_path / "rates.hdf5"
        net.write_table(out, label="via write_table", err_tol=None)

        assert _root_attrs(out)["network_name"] == "via write_table"

    def test_reaction_coeff_layout_is_unchanged(self, net, tmp_path):
        out = tmp_path / "rates.hdf5"
        net.to_hdf5(out, err_tol=None)

        group = HDF5().to_dict(out)["reaction_coeff"]
        assert set(group["_attrs"]) == {
            "input_names",
            "input_units",
            "xlo",
            "xhi",
            "spacing",
        }
        assert {k for k in group if k != "_attrs"} == {
            "output_names",
            "output_units",
            "data",
        }
        assert group["output_names"]["_data"].tolist() == [
            "2_body rate coefficient: {'H': 2} --> {'H2': 1}"
        ]


def test_text_header_names_the_network(net, tmp_path):
    out = tmp_path / "rates.txt"
    net.to_txt(out, label="my network", err_tol=None)

    assert "# Network name: my network\n" in out.read_text()
