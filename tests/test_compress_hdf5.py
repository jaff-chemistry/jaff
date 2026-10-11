# ABOUTME: Tests for jaff._utils.compress_hdf5: gzip copies of local data files for
# ABOUTME: upload to the mirror, sha256 output lines and registry.txt updates.
import hashlib
from pathlib import Path

import h5py
import numpy as np
import pytest

from jaff._utils.compress_hdf5 import compress_file, main, update_registry, verify_hdf5
from jaff.drivers._install import SOURCE_ATTR, copy_hdf5, decompress_hdf5

STR = h5py.string_dtype(encoding="utf-8")
CDT = np.dtype([("a", "f8"), ("b", "i4")])
MIN = 4096


def _make_uncompressed(path: Path) -> None:
    """Uncompressed HDF5 file exercising every dataset/attribute kind we care about."""
    with h5py.File(path, "w") as f:
        f.attrs["database"] = "leiden"
        f.attrs[SOURCE_ATTR] = "deadbeef"
        f.attrs.create("tags", np.array(["a", "bb"], dtype=object), dtype=STR)
        f.attrs.create("rec", np.array([(1.0, 2)], dtype=CDT), dtype=CDT)
        g = f.create_group("CO._PHOTON__C.O")
        g.attrs.create("reactants", np.array(["CO", "_PHOTON"], dtype=object), dtype=STR)
        g.attrs["n"] = np.int64(3)
        big = g.create_dataset("big", data=np.linspace(1, 10, 4 * MIN))
        big.attrs["unit"] = "eV"
        big.attrs.create("rec", np.array([(3.0, 4)], dtype=CDT), dtype=CDT)
        g.create_dataset("small", data=np.arange(4.0))
        g.create_dataset("scalar", data=np.float64(2.5))
        g.create_dataset("empty", data=np.zeros(0))
        g.create_dataset("null", data=h5py.Empty("f8"))
        g.create_dataset(
            "names", data=np.array(["x" * 50] * 1000, dtype=object), dtype=STR
        )
        g.create_dataset(
            "table", data=np.array([(float(i), i) for i in range(MIN)], dtype=CDT)
        )
        g.create_dataset("grow", data=np.arange(5.0), maxshape=(None,), chunks=True)


def _walk(f: h5py.File) -> dict:
    out = {}
    f.visititems(lambda name, obj: out.__setitem__(name, obj))
    return out


def _assert_attrs_equal(a, b, skip=()) -> None:
    assert set(a.attrs) - set(skip) == set(b.attrs) - set(skip)
    for k in set(b.attrs) - set(skip):
        np.testing.assert_array_equal(a.attrs[k], b.attrs[k])
        assert a.attrs.get_id(k).dtype == b.attrs.get_id(k).dtype


def _assert_same_content(src: Path, dst: Path, skip_root=()) -> None:
    with h5py.File(src) as s, h5py.File(dst) as d:
        _assert_attrs_equal(s, d, skip_root)
        so, do = _walk(s), _walk(d)
        assert set(so) == set(do)
        for name, sobj in so.items():
            dobj = do[name]
            _assert_attrs_equal(sobj, dobj)
            if isinstance(sobj, h5py.Dataset):
                assert dobj.dtype == sobj.dtype and dobj.shape == sobj.shape
                if sobj.shape is not None:
                    np.testing.assert_array_equal(dobj[()], sobj[()])


class TestCopyHdf5:
    def test_compresses_only_large_fixed_size_datasets(self, tmp_path):
        src, dst = tmp_path / "src.hdf5", tmp_path / "dst.hdf5"
        _make_uncompressed(src)
        copy_hdf5(
            src, dst, compress_min_bytes=MIN, level=4, drop_root_attrs=(SOURCE_ATTR,)
        )
        _assert_same_content(src, dst, skip_root=(SOURCE_ATTR,))
        with h5py.File(dst) as d:
            assert SOURCE_ATTR not in d.attrs
            g = d["CO._PHOTON__C.O"]
            for name in ("big", "table"):
                assert g[name].compression == "gzip"
                assert g[name].compression_opts == 4
                assert g[name].shuffle
            for name in ("small", "scalar", "empty", "null", "names"):
                assert g[name].compression is None
                assert g[name].chunks is None
            assert g["grow"].compression is None
            assert g["grow"].maxshape == (None,)

    def test_root_attrs_written(self, tmp_path):
        src, dst = tmp_path / "src.hdf5", tmp_path / "dst.hdf5"
        _make_uncompressed(src)
        copy_hdf5(src, dst, root_attrs={"extra": "yes"})
        with h5py.File(dst) as d:
            assert d.attrs["extra"] == "yes"
            assert d.attrs[SOURCE_ATTR] == "deadbeef"

    def test_full_cycle_back_to_uncompressed(self, tmp_path):
        src, gz, back = (tmp_path / n for n in ("src.hdf5", "gz.hdf5", "back.hdf5"))
        _make_uncompressed(src)
        copy_hdf5(src, gz, compress_min_bytes=MIN, drop_root_attrs=(SOURCE_ATTR,))
        decompress_hdf5(gz, back, "cafe")
        _assert_same_content(src, back, skip_root=(SOURCE_ATTR,))
        with h5py.File(back) as b:
            assert b.attrs[SOURCE_ATTR] == "cafe"
            for name, obj in _walk(b).items():
                if isinstance(obj, h5py.Dataset):
                    assert obj.compression is None
                    if not name.endswith("grow"):
                        assert obj.chunks is None


class TestCompressFile:
    def test_hdf5_returns_output_sha(self, tmp_path):
        src, dst = tmp_path / "src.hdf5", tmp_path / "out" / "sub" / "src.hdf5"
        _make_uncompressed(src)
        sha = compress_file(src, dst, min_size=MIN, level=4, verify=True)
        assert sha == hashlib.sha256(dst.read_bytes()).hexdigest()
        assert [p.name for p in dst.parent.iterdir()] == ["src.hdf5"]  # no temp left

    def test_failure_leaves_no_temp(self, tmp_path):
        src, dst = tmp_path / "src.hdf5", tmp_path / "out" / "src.hdf5"
        with h5py.File(src, "w") as f:
            f["a"] = np.arange(3)
            f["b"] = h5py.SoftLink("/a")
        with pytest.raises(ValueError, match="link"):
            compress_file(src, dst, min_size=MIN, level=4, verify=True)
        assert list(dst.parent.iterdir()) == []

    def test_verify_detects_mismatch(self, tmp_path):
        a, b = tmp_path / "a.hdf5", tmp_path / "b.hdf5"
        _make_uncompressed(a)
        _make_uncompressed(b)
        with h5py.File(b, "a") as f:
            f["CO._PHOTON__C.O/small"][0] = -1.0
        with pytest.raises(ValueError, match="small"):
            verify_hdf5(a, b)


class TestMain:
    def test_csv_passthrough_prints_sha(self, tmp_path, capsys):
        root, out = tmp_path / "data", tmp_path / "upload"
        csv = root / "xsecs" / "verner.csv"
        csv.parent.mkdir(parents=True)
        csv.write_text("a,b\n1,2\n")
        main([str(csv), "--outdir", str(out), "--root", str(root)])
        dst = out / "xsecs" / "verner.csv"
        assert dst.read_bytes() == csv.read_bytes()
        sha = hashlib.sha256(dst.read_bytes()).hexdigest()
        assert capsys.readouterr().out.strip() == f"xsecs/verner.csv {sha}"

    def test_hdf5_and_registry(self, tmp_path, capsys):
        root, out = tmp_path / "data", tmp_path / "upload"
        h5 = root / "xsecs" / "leiden.hdf5"
        h5.parent.mkdir(parents=True)
        _make_uncompressed(h5)
        reg = out / "registry.txt"
        main(
            [str(h5), "--outdir", str(out), "--root", str(root)]
            + ["--min-size", str(MIN), "--registry", str(reg)]
        )
        sha = hashlib.sha256((out / "xsecs" / "leiden.hdf5").read_bytes()).hexdigest()
        assert capsys.readouterr().out.strip() == f"xsecs/leiden.hdf5 {sha}"
        assert reg.read_text() == f"xsecs/leiden.hdf5 {sha}\n"

    def test_file_outside_root_errors(self, tmp_path):
        root = tmp_path / "data"
        root.mkdir()
        other = tmp_path / "x.csv"
        other.write_text("x")
        with pytest.raises(SystemExit):
            main([str(other), "--outdir", str(tmp_path / "o"), "--root", str(root)])

    @pytest.mark.parametrize("level", ["10", "-1", "x"])
    def test_invalid_level_errors(self, tmp_path, level, capsys):
        f = tmp_path / "x.csv"
        f.write_text("x")
        with pytest.raises(SystemExit) as exc:
            main(
                [str(f), "--outdir", str(tmp_path / "o"), "--root", str(tmp_path)]
                + ["--level", level]
            )
        assert exc.value.code == 2
        assert "level" in capsys.readouterr().err


class TestUpdateRegistry:
    def test_replace_in_place_and_append(self, tmp_path):
        reg = tmp_path / "registry.txt"
        reg.write_text("a.csv 111\nxsecs/leiden.hdf5 222\nb.csv 333\n")
        update_registry(reg, {"xsecs/leiden.hdf5": "999", "dust/dust.hdf5": "444"})
        assert reg.read_text() == (
            "a.csv 111\nxsecs/leiden.hdf5 999\nb.csv 333\ndust/dust.hdf5 444\n"
        )

    def test_creates_missing(self, tmp_path):
        reg = tmp_path / "new" / "registry.txt"
        update_registry(reg, {"a.csv": "1"})
        assert reg.read_text() == "a.csv 1\n"
