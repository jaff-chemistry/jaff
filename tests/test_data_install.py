# ABOUTME: Tests for jaff.drivers._install: uncompressed HDF5 copies and the
# ABOUTME: ownership rules that decide whether an installed data file is replaced.
import logging
import os
from pathlib import Path
from unittest.mock import patch

import h5py
import numpy as np
import pytest

from jaff.drivers._install import SOURCE_ATTR, decompress_hdf5, file_sha256, install_file

STR = h5py.string_dtype(encoding="utf-8")


def _make_compressed(path: Path) -> None:
    """HDF5 file exercising gzip/chunked datasets and awkward attribute types."""
    with h5py.File(path, "w") as f:
        f.attrs["database"] = "leiden"
        f.attrs.create("tags", np.array(["a", "bb"], dtype=object), dtype=STR)
        g = f.create_group("CO._PHOTON__C.O")
        g.attrs.create("reactants", np.array(["CO", "_PHOTON"], dtype=object), dtype=STR)
        g.attrs["n"] = np.int64(3)
        e = g.create_dataset(
            "photon_energy",
            data=np.linspace(1, 10, 1000),
            compression="gzip",
            compression_opts=4,
            chunks=True,
        )
        e.attrs["unit"] = "eV"
        g.create_dataset("scalar", data=np.float64(2.5))
        g.create_dataset("empty", data=np.zeros(0))
        g.create_dataset("names", data=np.array(["x", "yy"], dtype=object), dtype=STR)
        cdt = np.dtype([("a", "f8"), ("b", "i4")])
        g.create_dataset(
            "table",
            data=np.array([(1.0, 2), (3.0, 4)], dtype=cdt),
            compression="gzip",
            chunks=True,
        )
        g.create_dataset(
            "grow", data=np.arange(5.0), maxshape=(None,), chunks=True, compression="gzip"
        )


def _walk(f: h5py.File) -> dict:
    out = {}

    def visit(name, obj):
        out[name] = obj

    f.visititems(visit)
    return out


class TestDecompress:
    def test_copy_is_identical_and_uncompressed(self, tmp_path):
        src, dst = tmp_path / "src.hdf5", tmp_path / "dst.hdf5"
        _make_compressed(src)
        decompress_hdf5(src, dst, "abc123")
        with h5py.File(src) as s, h5py.File(dst) as d:
            assert d.attrs[SOURCE_ATTR] == "abc123"
            for k in s.attrs:
                np.testing.assert_array_equal(s.attrs[k], d.attrs[k])
            so, do = _walk(s), _walk(d)
            assert set(so) == set(do)
            for name, sobj in so.items():
                dobj = do[name]
                for k in sobj.attrs:
                    np.testing.assert_array_equal(sobj.attrs[k], dobj.attrs[k])
                    assert sobj.attrs.get_id(k).dtype == dobj.attrs.get_id(k).dtype
                if isinstance(sobj, h5py.Dataset):
                    assert dobj.dtype == sobj.dtype and dobj.shape == sobj.shape
                    np.testing.assert_array_equal(dobj[()], sobj[()])
                    assert dobj.compression is None
                    if name.endswith("grow"):
                        assert dobj.maxshape == (None,)
                    else:
                        assert dobj.chunks is None  # contiguous

    def test_soft_link_rejected(self, tmp_path):
        src = tmp_path / "src.hdf5"
        with h5py.File(src, "w") as f:
            f["a"] = np.arange(3)
            f["b"] = h5py.SoftLink("/a")
        with pytest.raises(ValueError, match="link"):
            decompress_hdf5(src, tmp_path / "dst.hdf5", "x")


def _logger():
    lg = logging.getLogger("install-test")
    return lg


class TestInstallHdf5:
    def test_missing_target_installed(self, tmp_path):
        src, dst = tmp_path / "dl" / "f.hdf5", tmp_path / "out" / "f.hdf5"
        src.parent.mkdir()
        _make_compressed(src)
        assert install_file(src, dst, "s1", _logger()) is True
        with h5py.File(dst) as d:
            assert d.attrs[SOURCE_ATTR] == "s1"

    def test_up_to_date_not_rewritten(self, tmp_path):
        src, dst = tmp_path / "f.hdf5", tmp_path / "out.hdf5"
        _make_compressed(src)
        install_file(src, dst, "s1", _logger())
        mtime = dst.stat().st_mtime_ns
        assert install_file(src, dst, "s1", _logger()) is False
        assert dst.stat().st_mtime_ns == mtime

    def test_stale_reinstalled(self, tmp_path):
        src, dst = tmp_path / "f.hdf5", tmp_path / "out.hdf5"
        _make_compressed(src)
        install_file(src, dst, "s1", _logger())
        assert install_file(src, dst, "s2", _logger()) is True
        with h5py.File(dst) as d:
            assert d.attrs[SOURCE_ATTR] == "s2"

    def test_user_supplied_kept_with_warning(self, tmp_path):
        src, dst = tmp_path / "f.hdf5", tmp_path / "out.hdf5"
        _make_compressed(src)
        with h5py.File(dst, "w") as f:
            f["mine"] = np.arange(2)
        lg = _logger()
        with patch.object(lg, "warning") as warn:
            assert install_file(src, dst, "s1", lg) is False
        warn.assert_called_once()
        assert "out.hdf5" in warn.call_args.args[0]
        with h5py.File(dst) as d:
            assert list(d) == ["mine"]

    def test_failed_install_keeps_old_file_and_no_temp(self, tmp_path):
        src, dst = tmp_path / "f.hdf5", tmp_path / "out.hdf5"
        _make_compressed(src)
        install_file(src, dst, "s1", _logger())
        before = dst.read_bytes()
        with patch("jaff.drivers._install._copy_group", side_effect=RuntimeError("boom")):
            with pytest.raises(RuntimeError):
                install_file(src, dst, "s2", _logger())
        assert dst.read_bytes() == before
        assert sorted(p.name for p in tmp_path.iterdir()) == ["f.hdf5", "out.hdf5"]


class TestInstallOther:
    def test_csv_copied_and_sidecar_written(self, tmp_path):
        src, dst = tmp_path / "dl.csv", tmp_path / "out" / "v.csv"
        src.write_text("a b\n1 2\n")
        assert install_file(src, dst, "s1", _logger()) is True
        assert dst.read_text() == "a b\n1 2\n"
        sidecar = src.with_name(src.name + ".installed-sha256")
        assert sidecar.read_text().split() == ["s1", file_sha256(dst)]

    def test_csv_up_to_date_and_stale(self, tmp_path):
        src, dst = tmp_path / "dl.csv", tmp_path / "v.csv"
        src.write_text("1\n")
        install_file(src, dst, "s1", _logger())
        assert install_file(src, dst, "s1", _logger()) is False
        src.write_text("2\n")
        assert install_file(src, dst, "s2", _logger()) is True
        assert dst.read_text() == "2\n"

    def test_csv_edited_by_user_kept(self, tmp_path):
        src, dst = tmp_path / "dl.csv", tmp_path / "v.csv"
        src.write_text("1\n")
        install_file(src, dst, "s1", _logger())
        dst.write_text("my edit\n")
        lg = _logger()
        with patch.object(lg, "warning") as warn:
            assert install_file(src, dst, "s2", lg) is False
        warn.assert_called_once()
        assert dst.read_text() == "my edit\n"

    def test_csv_without_sidecar_kept(self, tmp_path):
        src, dst = tmp_path / "dl.csv", tmp_path / "v.csv"
        src.write_text("1\n")
        dst.write_text("pre-existing\n")
        assert install_file(src, dst, "s1", _logger()) is False
        assert dst.read_text() == "pre-existing\n"


def test_file_sha256_matches_hashlib(tmp_path):
    import hashlib

    p = tmp_path / "x.bin"
    p.write_bytes(os.urandom(3_000_000))
    assert file_sha256(p) == hashlib.sha256(p.read_bytes()).hexdigest()
