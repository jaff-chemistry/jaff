# ABOUTME: Tests for drivers.csv: csv_to_df reads a CSV given as a str or Path; csv_to_hdf5
# ABOUTME: checks flatten_columns and existing datasets before writing anything

import inspect

import h5py
import pytest

from jaff.drivers.csv import csv_to_df, csv_to_hdf5


@pytest.fixture
def grid(tmp_path):
    path = tmp_path / "grid.csv"
    path.write_text("T,n,k\n10,1,0.1\n20,2,0.2\n", encoding="utf-8")

    return path


def test_a_str_path_reads_like_a_path(grid):
    assert csv_to_df(str(grid)).equals(csv_to_df(grid))
    assert csv_to_df(str(grid))["k"].tolist() == pytest.approx([0.1, 0.2])


def test_a_missing_str_path_is_file_not_found(tmp_path):
    with pytest.raises(FileNotFoundError):
        csv_to_df(str(tmp_path / "absent.csv"))


def _contents(path):
    with h5py.File(path) as f:
        names = []
        f.visit(names.append)

        return {n: f[n][()].tolist() for n in names if isinstance(f[n], h5py.Dataset)}


class TestCsvToHdf5:
    def test_unknown_flatten_column_writes_no_file(self, grid, tmp_path):
        out = tmp_path / "out.h5"

        with pytest.raises(ValueError, match="flatten_columns"):
            csv_to_hdf5(grid, out, flatten_columns=["T", "x"])

        assert not out.exists()

    def test_unknown_flatten_column_leaves_an_existing_file_alone(self, grid, tmp_path):
        out = tmp_path / "out.h5"
        with h5py.File(out, "w") as f:
            f["keep"] = [1, 2]

        with pytest.raises(ValueError):
            csv_to_hdf5(grid, out, flatten_columns=["x"])

        with h5py.File(out) as f:
            assert list(f) == ["keep"]

    def test_unknown_flatten_column_message_names_only_the_missing(self, grid, tmp_path):
        with pytest.raises(ValueError) as info:
            csv_to_hdf5(grid, tmp_path / "out.h5", flatten_columns=["T", "x"])

        assert str(info.value) == "flatten_columns not in the CSV: x (columns: T, n, k)"

    def test_converting_twice_raises_and_keeps_the_first_write(self, grid, tmp_path):
        out = tmp_path / "out.h5"
        csv_to_hdf5(grid, out)
        first = _contents(out)

        with pytest.raises(ValueError, match="T, n, k"):
            csv_to_hdf5(grid, out)

        assert _contents(out) == first

    def test_a_partial_clash_writes_none_of_the_columns(self, grid, tmp_path):
        out = tmp_path / "out.h5"
        with h5py.File(out, "w") as f:
            f["grid/n"] = [9, 9]

        with pytest.raises(ValueError, match="already holds: n"):
            csv_to_hdf5(grid, out)

        assert _contents(out) == {"grid/n": [9, 9]}

    def test_flatten_columns_has_no_mutable_default(self):
        default = inspect.signature(csv_to_hdf5).parameters["flatten_columns"].default

        assert default is None

    def test_there_is_no_pytables_mode(self):
        # as_table needed PyTables (not a dependency) and nothing in jaff reads it.
        assert "as_table" not in inspect.signature(csv_to_hdf5).parameters
