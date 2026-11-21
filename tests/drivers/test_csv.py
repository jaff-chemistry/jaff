# ABOUTME: Tests for drivers.csv: csv_to_df reads a CSV given as a str or Path, and a
# ABOUTME: missing file is a FileNotFoundError

import pytest

from jaff.drivers.csv import csv_to_df


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
