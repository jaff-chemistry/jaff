# ABOUTME: Tests for ConfigTable HDF5 conversions, incl. composite datasets
# ABOUTME: Covers single-move passthrough, compound (named/stem) and ndarray composites

from pathlib import Path

import numpy as np
import pytest

from jaff.cli.jaffgen._config_table import ConfigTable
from jaff.drivers import HDF5


def _leaf(data: np.ndarray, dtype: str = "f64") -> dict:
    return {"_kind": "linear", "_data": data, "_dtype": dtype, "_attrs": {}}


@pytest.fixture
def source_h5(tmp_path) -> Path:
    """A small source HDF5 with three data columns and two axis datasets."""
    src = tmp_path / "source.hdf5"
    tree = {
        "co": {
            "A": _leaf(np.arange(6, dtype=np.float64)),
            "B": _leaf(np.arange(6, dtype=np.float64) * 2.0),
            "C": _leaf(np.arange(6, dtype=np.float64) * 3.0),
            "x0": _leaf(np.array([1.0, 2.0, 3.0])),  # len 3
            "x1": _leaf(np.array([10.0, 20.0])),  # len 2  -> 3*2 = 6
        }
    }
    HDF5().from_dict(src, tree, mode="w")
    return src


def _make(table_dict, source_h5, tmp_path):
    # network_file only resolves the "default" alias / stem; unused here.
    return ConfigTable(table_dict, tmp_path / "cfg.toml", tmp_path / "net.jet")


def test_single_h5path_move_unchanged(source_h5, tmp_path):
    """A bare-string h5path keeps the original single-move behaviour."""
    ct = _make(
        {
            "source": {"path": str(source_h5)},
            "target": {"path": "out.hdf5", "/temperature": {"h5path": "/co/x0"}},
        },
        source_h5,
        tmp_path,
    )
    out = ct.parse()
    assert "/temperature" in out.flatten() or "temperature" in out
    flat = out.flatten()
    np.testing.assert_array_equal(flat["/temperature"]["_data"], [1.0, 2.0, 3.0])
    # Source path was moved away.
    assert "/co/x0" not in flat


def test_compound_stem_names(source_h5, tmp_path):
    """List h5path with default type -> compound; columns named by path stem."""
    ct = _make(
        {
            "source": {"path": str(source_h5)},
            "target": {
                "path": "out.hdf5",
                "/c0/data": {"h5path": ["/co/A", "/co/B", "/co/C"]},
            },
        },
        source_h5,
        tmp_path,
    )
    flat = ct.parse().flatten()
    leaf = flat["/c0/data"]
    assert leaf["_kind"] == "compound"
    assert set(leaf["_data"].dtype.names) == {"A", "B", "C"}
    np.testing.assert_array_equal(leaf["_data"]["B"], np.arange(6) * 2.0)
    assert leaf["_dtype"] == {"A": "f64", "B": "f64", "C": "f64"}
    # Consumed sources are dropped.
    for p in ("/co/A", "/co/B", "/co/C"):
        assert p not in flat


def test_compound_explicit_names(source_h5, tmp_path):
    """`names` overrides the stem-derived column names."""
    ct = _make(
        {
            "source": {"path": str(source_h5)},
            "target": {
                "path": "out.hdf5",
                "/c0/data": {
                    "h5path": ["/co/A", "/co/B"],
                    "type": "compound",
                    "names": ["alpha", "beta"],
                },
            },
        },
        source_h5,
        tmp_path,
    )
    leaf = ct.parse().flatten()["/c0/data"]
    assert list(leaf["_data"].dtype.names) == ["alpha", "beta"]


def test_ndarray_shape_and_stack(source_h5, tmp_path):
    """ndarray composite -> (Ncols, *xlens) reshaped and stacked."""
    ct = _make(
        {
            "source": {"path": str(source_h5)},
            "target": {
                "path": "out.hdf5",
                "/c0/data": {
                    "h5path": ["/co/A", "/co/B"],
                    "type": "ndarray",
                    "regrid": {"x": ["/co/x0", "/co/x1"]},
                },
            },
        },
        source_h5,
        tmp_path,
    )
    leaf = ct.parse().flatten()["/c0/data"]
    assert leaf["_kind"] == "linear"
    assert leaf["_data"].shape == (2, 3, 2)  # (Ncols, len x0, len x1)
    np.testing.assert_array_equal(
        leaf["_data"][0], np.arange(6, dtype=np.float64).reshape(3, 2)
    )
    np.testing.assert_array_equal(
        leaf["_data"][1], (np.arange(6, dtype=np.float64) * 2.0).reshape(3, 2)
    )


def test_ndarray_single_axis_string(source_h5, tmp_path):
    """A scalar regrid.x string is accepted as a single axis."""
    ct = _make(
        {
            "source": {"path": str(source_h5)},
            "target": {
                "path": "out.hdf5",
                # 3-element column over a single length-3 axis.
                "/c0/data": {
                    "h5path": ["/co/x0"],
                    "type": "ndarray",
                    "regrid": {"x": "/co/x0"},
                },
            },
        },
        source_h5,
        tmp_path,
    )
    leaf = ct.parse().flatten()["/c0/data"]
    assert leaf["_data"].shape == (1, 3)


def test_missing_leading_slash_normalised(source_h5, tmp_path):
    """Source paths without a leading slash still resolve."""
    ct = _make(
        {
            "source": {"path": str(source_h5)},
            "target": {
                "path": "out.hdf5",
                "/c0/data": {"h5path": ["co/A", "co/B"]},
            },
        },
        source_h5,
        tmp_path,
    )
    leaf = ct.parse().flatten()["/c0/data"]
    assert set(leaf["_data"].dtype.names) == {"A", "B"}


def test_unequal_length_raises(tmp_path):
    src = tmp_path / "source.hdf5"
    HDF5().from_dict(
        src,
        {
            "co": {
                "A": _leaf(np.arange(4, dtype=np.float64)),
                "B": _leaf(np.arange(3, dtype=np.float64)),
            }
        },
        mode="w",
    )
    ct = _make(
        {
            "source": {"path": str(src)},
            "target": {
                "path": "out.hdf5",
                "/c0/data": {"h5path": ["/co/A", "/co/B"]},
            },
        },
        src,
        tmp_path,
    )
    with pytest.raises(ValueError, match="share a length"):
        ct.parse()


def test_names_length_mismatch_raises(source_h5, tmp_path):
    ct = _make(
        {
            "source": {"path": str(source_h5)},
            "target": {
                "path": "out.hdf5",
                "/c0/data": {
                    "h5path": ["/co/A", "/co/B"],
                    "names": ["only_one"],
                },
            },
        },
        source_h5,
        tmp_path,
    )
    with pytest.raises(ValueError, match="names length"):
        ct.parse()


def test_ndarray_missing_regrid_raises(source_h5, tmp_path):
    ct = _make(
        {
            "source": {"path": str(source_h5)},
            "target": {
                "path": "out.hdf5",
                "/c0/data": {"h5path": ["/co/A"], "type": "ndarray"},
            },
        },
        source_h5,
        tmp_path,
    )
    with pytest.raises(ValueError, match="regrid.x"):
        ct.parse()


def test_ndarray_bad_reshape_raises(source_h5, tmp_path):
    ct = _make(
        {
            "source": {"path": str(source_h5)},
            "target": {
                "path": "out.hdf5",
                # A (len-6) over axes 3*3=9 -> mismatch.
                "/c0/data": {
                    "h5path": ["/co/A"],
                    "type": "ndarray",
                    "regrid": {"x": ["/co/x0", "/co/x0"]},
                },
            },
        },
        source_h5,
        tmp_path,
    )
    with pytest.raises(ValueError, match="axis grid"):
        ct.parse()


def test_attrs_scalar_ref_and_literal(source_h5, tmp_path):
    """A single computed ref resolves; scalar literals pass through."""
    ct = _make(
        {
            "source": {"path": str(source_h5)},
            "target": {
                "path": "out.hdf5",
                "/temperature": {
                    "h5path": "/co/x0",
                    "attrs": {"tmax": "/temperature.max", "units": "K", "ndim": 1},
                },
            },
        },
        source_h5,
        tmp_path,
    )
    attrs = ct.parse().flatten()["/temperature"]["_attrs"]
    assert attrs["tmax"] == 3.0
    assert attrs["units"] == "K"
    assert attrs["ndim"] == 1


def test_attrs_list_mixed(source_h5, tmp_path):
    """List attrs resolve element-wise: refs computed, literals verbatim."""
    ct = _make(
        {
            "source": {"path": str(source_h5)},
            "target": {
                "path": "out.hdf5",
                # /co is a group heading holding only attrs; axes stay in-tree.
                "/co/data": {"h5path": ["/co/A", "/co/B"], "type": "ndarray",
                             "regrid": {"x": ["/co/x0", "/co/x1"]}},
                "/co": {
                    "attrs": {
                        "Ndim": 2,
                        "Nx": ["/co/x0.length", "/co/x1.length"],
                        "spacing": ["log", "log"],
                        "xlo": ["/co/x0.min", "/co/x1.min"],
                        "xhigh": ["/co/x0.max", "/co/x1.max"],
                    }
                },
            },
        },
        source_h5,
        tmp_path,
    )
    attrs = ct.parse().nested()["co"]["_attrs"]
    assert attrs["Ndim"] == 2
    assert list(attrs["Nx"]) == [3, 2]
    assert list(attrs["spacing"]) == ["log", "log"]
    assert list(attrs["xlo"]) == [1.0, 10.0]
    assert list(attrs["xhigh"]) == [3.0, 20.0]


def test_attrs_bad_ref_raises(source_h5, tmp_path):
    ct = _make(
        {
            "source": {"path": str(source_h5)},
            "target": {
                "path": "out.hdf5",
                "/temperature": {
                    "h5path": "/co/x0",
                    "attrs": {"bad": "/nope/missing.max"},
                },
            },
        },
        source_h5,
        tmp_path,
    )
    with pytest.raises(ValueError, match="not found in target tree"):
        ct.parse()


def test_unknown_type_raises(source_h5, tmp_path):
    ct = _make(
        {
            "source": {"path": str(source_h5)},
            "target": {
                "path": "out.hdf5",
                "/c0/data": {"h5path": ["/co/A"], "type": "bogus"},
            },
        },
        source_h5,
        tmp_path,
    )
    with pytest.raises(ValueError, match="Unknown composite type"):
        ct.parse()


def test_relative_source_path_resolves_against_config_dir(
    source_h5, tmp_path, monkeypatch
):
    """A relative source.path is config-relative, not CWD-relative."""
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    # Same-named decoy in the CWD with different data must not be picked up.
    decoy = {"co": {"x0": _leaf(np.array([-1.0, -2.0, -3.0]))}}
    HDF5().from_dict(elsewhere / source_h5.name, decoy, mode="w")
    monkeypatch.chdir(elsewhere)

    ct = _make(
        {
            "source": {"path": source_h5.name},
            "target": {"path": "out.hdf5", "/temperature": {"h5path": "/co/x0"}},
        },
        source_h5,
        tmp_path,
    )
    assert ct.source_props["path"] == source_h5
    flat = ct.parse().flatten()
    np.testing.assert_array_equal(flat["/temperature"]["_data"], [1.0, 2.0, 3.0])


def test_default_source_with_dotted_network_name(source_h5, tmp_path):
    """"default" keeps every dot of the network basename: net.v2.jet -> net.v2.hdf5."""
    net_dir = tmp_path / "nets"
    net_dir.mkdir()
    default_h5 = net_dir / "net.v2.hdf5"
    default_h5.write_bytes(source_h5.read_bytes())

    ct = ConfigTable(
        {
            "source": {"path": "default"},
            "target": {"path": "out.hdf5", "/temperature": {"h5path": "/co/x0"}},
        },
        tmp_path / "cfg.toml",
        net_dir / "net.v2.jet",
    )
    assert ct.source_props["path"] == default_h5
    flat = ct.parse().flatten()
    np.testing.assert_array_equal(flat["/temperature"]["_data"], [1.0, 2.0, 3.0])
