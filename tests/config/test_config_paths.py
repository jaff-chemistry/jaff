# ABOUTME: Tests for resolve_network_path: network files, predefined network names
# ABOUTME: and the errors for a missing file or an ambiguous predefined directory

from pathlib import Path

import pytest

import jaff.config._config as paths
from jaff.config import NETWORKS_DIR, resolve_network_path
from jaff.errors import ParserError


def test_file_resolves_to_absolute_path(tmp_path, monkeypatch) -> None:
    (tmp_path / "net.dat").write_text("")
    monkeypatch.chdir(tmp_path)

    assert resolve_network_path("net.dat") == tmp_path / "net.dat"


def test_predefined_name_resolves_to_its_jet_file() -> None:
    assert resolve_network_path("GOW") == (NETWORKS_DIR / "GOW" / "GOW.jet").resolve()


def test_path_ending_in_predefined_name_is_the_predefined_network() -> None:
    expected = (NETWORKS_DIR / "GOW" / "GOW.jet").resolve()

    assert resolve_network_path(Path("/some/dir/GOW")) == expected


def test_missing_file(tmp_path) -> None:
    with pytest.raises(FileNotFoundError, match="nope.dat"):
        resolve_network_path(tmp_path / "nope.dat")


@pytest.mark.parametrize(("jets", "error"), [(0, FileNotFoundError), (2, ParserError)])
def test_predefined_dir_needs_exactly_one_jet(tmp_path, monkeypatch, jets, error) -> None:
    (tmp_path / "mynet").mkdir()
    for i in range(jets):
        (tmp_path / "mynet" / f"n{i}.jet").write_text("")
    monkeypatch.setattr(paths, "NETWORKS_DIR", tmp_path)

    with pytest.raises(error, match="mynet"):
        resolve_network_path("mynet")
