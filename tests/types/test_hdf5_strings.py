# ABOUTME: HDF5Dict reads text datasets and compound text fields back as str
# ABOUTME: (variable- and fixed-length), so from_dict -> to_dict round-trips

import h5py
import numpy as np
import pytest

from jaff.drivers import HDF5
from jaff.types import HDF5Dict


@pytest.mark.parametrize(
    "dtype", [h5py.string_dtype(), np.dtype("S4")], ids=["vlen", "fixed"]
)
def test_linear_text_reads_as_str(tmp_path, dtype):
    path = tmp_path / "s.h5"
    with h5py.File(path, "w") as f:
        f.create_dataset("names", data=np.array(["H2", "CO"]).astype(dtype))

    node = HDF5Dict(path)["names"]

    assert node["_data"].tolist() == ["H2", "CO"]
    assert node["_dtype"] == "s"


@pytest.mark.parametrize(
    "dtype", [h5py.string_dtype(), np.dtype("S4")], ids=["vlen", "fixed"]
)
def test_compound_text_field_reads_as_str(tmp_path, dtype):
    path = tmp_path / "s.h5"
    table = np.array([("H2", 1.0), ("CO", 2.0)], dtype=[("sp", dtype), ("k", "f8")])
    with h5py.File(path, "w") as f:
        f.create_dataset("tab", data=table)

    node = HDF5Dict(path)["tab"]

    assert node["_data"]["sp"].tolist() == ["H2", "CO"]
    assert node["_data"]["k"].tolist() == [1.0, 2.0]
    assert node["_dtype"] == {"sp": "s", "k": "f64"}


def test_text_round_trips(tmp_path):
    path = tmp_path / "s.h5"
    HDF5().from_dict(
        path,
        {"names": {"_kind": "linear", "_data": ["H2", "CO"], "_dtype": "s"}},
        mode="w",
    )

    assert HDF5().to_dict(path)["names"]["_data"].tolist() == ["H2", "CO"]
