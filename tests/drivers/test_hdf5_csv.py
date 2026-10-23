# ABOUTME: Tests for HDF5.to_csv output-path injectivity
# ABOUTME: Same-named sibling descendants and mixed linear/compound must not collide

from pathlib import Path

import h5py
import numpy as np
import pytest

from jaff.drivers import HDF5


def _read(path: Path) -> str:
    return path.read_text()


def test_same_named_sibling_groups_do_not_overwrite(tmp_path):
    """a/rates/k and b/rates/k must both survive, not overwrite each other."""
    src = tmp_path / "in.h5"
    with h5py.File(src, "w") as f:
        f["a/rates/k"] = [1.0, 2.0]
        f["b/rates/k"] = [3.0, 4.0]

    out = tmp_path / "csv"
    HDF5().to_csv(src, out)

    files = sorted(p.relative_to(out).as_posix() for p in out.rglob("*.csv"))
    assert len(files) == 2, files

    a_text = "".join(_read(p) for p in out.rglob("*.csv") if "a/" in p.as_posix())
    b_text = "".join(_read(p) for p in out.rglob("*.csv") if "b/" in p.as_posix())
    assert "1.0" in a_text and "2.0" in a_text
    assert "3.0" in b_text and "4.0" in b_text


def test_mixed_linear_and_compound_in_same_group(tmp_path):
    """A group with both a linear and a compound dataset yields two files."""
    src = tmp_path / "in.h5"
    comp = np.array([(1, 2.0), (3, 4.0)], dtype=[("i", "i4"), ("x", "f8")])
    with h5py.File(src, "w") as f:
        f["g/y"] = [10.0, 20.0]
        f["g/tab"] = comp

    out = tmp_path / "csv"
    HDF5().to_csv(src, out)

    csvs = list(out.rglob("*.csv"))
    assert len(csvs) == 2, [p.as_posix() for p in csvs]
    # Compound survives with its own file.
    tab_text = "".join(_read(p) for p in csvs if "tab" in p.name)
    assert "2.0" in tab_text and "4.0" in tab_text
    # Linear group table survives.
    lin_text = "".join(_read(p) for p in csvs if "tab" not in p.name)
    assert "10.0" in lin_text and "20.0" in lin_text


def test_compound_siblings_across_groups_do_not_overwrite(tmp_path):
    """Same-named compound datasets in different groups must both survive."""
    src = tmp_path / "in.h5"
    comp1 = np.array([(1,)], dtype=[("v", "i4")])
    comp2 = np.array([(2,)], dtype=[("v", "i4")])
    with h5py.File(src, "w") as f:
        f["a/tab"] = comp1
        f["b/tab"] = comp2

    out = tmp_path / "csv"
    HDF5().to_csv(src, out)

    csvs = list(out.rglob("*.csv"))
    assert len(csvs) == 2, [p.as_posix() for p in csvs]


def test_compound_keyed_as_reserved_stem_collides(tmp_path):
    """A compound dataset keyed like the reserved linear table must be caught."""
    stem = HDF5._LINEAR_TABLE_STEM
    src = tmp_path / "in.h5"
    comp = np.array([(1,)], dtype=[("v", "i4")])
    with h5py.File(src, "w") as f:
        f["g/y"] = [1.0]           # linear -> g/_group.csv
        f[f"g/{stem}"] = comp      # compound -> g/_group.csv  (collision)

    out = tmp_path / "csv"
    with pytest.raises(ValueError):
        HDF5().to_csv(src, out)
