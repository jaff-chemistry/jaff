# ABOUTME: End-to-end tests for jaff.drivers.pooch.Pooch against a local HTTP server:
# ABOUTME: compressed downloads cached in .downloads/, uncompressed installs, migration.
import http.server
import shutil
import threading
from pathlib import Path

import h5py
import numpy as np
import pytest

from jaff.drivers._install import SOURCE_ATTR, file_sha256
from jaff.drivers.pooch import Pooch


@pytest.fixture
def server(tmp_path):
    """Serve tmp_path/'remote' over HTTP; record requested paths."""
    remote = tmp_path / "remote"
    (remote / "xsecs").mkdir(parents=True)
    requests: list[str] = []

    class Handler(http.server.SimpleHTTPRequestHandler):
        def __init__(self, *a, **k):
            super().__init__(*a, directory=str(remote), **k)

        def log_message(self, *a):
            pass

        def do_GET(self):
            requests.append(self.path)
            super().do_GET()

    httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    t = threading.Thread(target=httpd.serve_forever, daemon=True)
    t.start()
    url = f"http://127.0.0.1:{httpd.server_address[1]}"
    yield remote, url, requests
    httpd.shutdown()


def _publish(remote: Path, values) -> str:
    """Write a gzip HDF5 + CSV to the remote and a registry; return the h5 sha."""
    h5 = remote / "xsecs" / "leiden.hdf5"
    with h5py.File(h5, "w") as f:
        f.create_dataset(
            "g/photon_energy",
            data=np.asarray(values, float),
            compression="gzip",
            chunks=True,
        )
    (remote / "xsecs" / "v.csv").write_text("a b\n1 2\n")
    (remote / "registry.txt").write_text(
        f"xsecs/leiden.hdf5 {file_sha256(h5)}\n"
        f"xsecs/v.csv {file_sha256(remote / 'xsecs' / 'v.csv')}\n"
    )
    return file_sha256(h5)


@pytest.fixture
def online(monkeypatch):
    monkeypatch.delenv("JAFF_OFFLINE", raising=False)


def test_first_fetch_installs_uncompressed(server, tmp_path, online):
    remote, url, _ = server
    sha = _publish(remote, [1, 2, 3])
    root = tmp_path / "data"
    p = Pooch(url, root)
    p.fetch_file("xsecs/leiden.hdf5")
    p.fetch_file("xsecs/v.csv")
    assert (root / ".downloads" / "xsecs" / "leiden.hdf5").exists()
    assert (root / ".downloads" / "registry.txt").exists()
    with h5py.File(root / "xsecs" / "leiden.hdf5") as f:
        assert f.attrs[SOURCE_ATTR] == sha
        assert f["g/photon_energy"].compression is None
        np.testing.assert_array_equal(f["g/photon_energy"][()], [1, 2, 3])
    assert (root / "xsecs" / "v.csv").read_text() == "a b\n1 2\n"


def test_second_fetch_no_rewrite_no_download(server, tmp_path, online):
    remote, url, requests = server
    _publish(remote, [1, 2, 3])
    root = tmp_path / "data"
    p = Pooch(url, root)
    p.fetch_file("xsecs/leiden.hdf5")
    mtime = (root / "xsecs" / "leiden.hdf5").stat().st_mtime_ns
    n = len(requests)
    p.fetch_file("xsecs/leiden.hdf5")
    assert (root / "xsecs" / "leiden.hdf5").stat().st_mtime_ns == mtime
    assert len(requests) == n


def test_registry_update_reinstalls(server, tmp_path, online):
    remote, url, _ = server
    _publish(remote, [1, 2, 3])
    root = tmp_path / "data"
    Pooch(url, root).fetch_file("xsecs/leiden.hdf5")
    sha2 = _publish(remote, [4, 5])
    Pooch._registry.clear()  # force a fresh instance that re-reads the registry
    Pooch(url, root).fetch_file("xsecs/leiden.hdf5")
    with h5py.File(root / "xsecs" / "leiden.hdf5") as f:
        assert f.attrs[SOURCE_ATTR] == sha2
        np.testing.assert_array_equal(f["g/photon_energy"][()], [4, 5])


def test_migration_moves_existing_verified_file(server, tmp_path, online):
    remote, url, requests = server
    sha = _publish(remote, [1, 2, 3])
    root = tmp_path / "data"
    (root / "xsecs").mkdir(parents=True)
    shutil.copy(remote / "xsecs" / "leiden.hdf5", root / "xsecs" / "leiden.hdf5")
    (root / "registry.txt").write_text((remote / "registry.txt").read_text())
    Pooch(url, root).fetch_file("xsecs/leiden.hdf5")
    assert "/xsecs/leiden.hdf5" not in requests  # moved, not downloaded
    assert file_sha256(root / ".downloads" / "xsecs" / "leiden.hdf5") == sha
    with h5py.File(root / "xsecs" / "leiden.hdf5") as f:
        assert f.attrs[SOURCE_ATTR] == sha
    assert not (root / "registry.txt").exists()


def test_user_supplied_file_kept(server, tmp_path, online):
    remote, url, _ = server
    _publish(remote, [1, 2, 3])
    root = tmp_path / "data"
    (root / "xsecs").mkdir(parents=True)
    with h5py.File(root / "xsecs" / "leiden.hdf5", "w") as f:
        f["mine"] = np.arange(4)
    Pooch(url, root).fetch_file("xsecs/leiden.hdf5")
    with h5py.File(root / "xsecs" / "leiden.hdf5") as f:
        assert list(f) == ["mine"]
    assert (root / ".downloads" / "xsecs" / "leiden.hdf5").exists()


def test_offline_does_nothing(server, tmp_path, monkeypatch):
    remote, url, requests = server
    _publish(remote, [1])
    monkeypatch.setenv("JAFF_OFFLINE", "1")
    root = tmp_path / "data"
    Pooch(url, root).fetch_file("xsecs/leiden.hdf5")
    assert requests == []
    assert not (root / "xsecs").exists()
