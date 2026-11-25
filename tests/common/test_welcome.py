# ABOUTME: The startup banner must print on any stdout encoding (e.g. Windows cp1252
# ABOUTME: when output is piped), degrading block-art characters instead of crashing

import builtins
import io

import pytest

from jaff.common import print_motd

# conftest's autouse _quiet fixture replaces print with a no-op for every test;
# grab the real one at import (collection) time, before that fixture runs.
REAL_PRINT = builtins.print


def _cp1252_stdout(raw: io.BytesIO) -> io.TextIOWrapper:
    """A stdout like Windows gives a piped or redirected process."""
    return io.TextIOWrapper(raw, encoding="cp1252")


@pytest.fixture(autouse=True)
def _banner_not_shown(monkeypatch):
    """Each test starts as a fresh process: no banner printed yet."""
    monkeypatch.setattr("jaff.common._welcome._SHOWN", False)


@pytest.mark.parametrize("mode", ["default", "jaffgen", "jaffx"])
def test_banner_prints_on_cp1252_stdout(mode, monkeypatch):
    raw = io.BytesIO()
    stream = _cp1252_stdout(raw)
    monkeypatch.setattr("builtins.print", REAL_PRINT)
    monkeypatch.setattr("sys.stdout", stream)

    print_motd(mode)

    stream.flush()
    assert raw.getvalue()


def test_network_loads_on_cp1252_stdout(make_network, monkeypatch):
    monkeypatch.setattr("builtins.print", REAL_PRINT)
    monkeypatch.setattr("sys.stdout", _cp1252_stdout(io.BytesIO()))

    net = make_network(["H + H -> H2 [10,1000] 1e-10"])

    assert net.reactions.count == 1


def test_banner_prints_once_per_process(monkeypatch):
    raw = io.BytesIO()
    stream = _cp1252_stdout(raw)
    monkeypatch.setattr("builtins.print", REAL_PRINT)
    monkeypatch.setattr("sys.stdout", stream)

    print_motd("jaffgen")
    stream.flush()
    first = raw.getvalue()
    print_motd()
    stream.flush()

    assert first
    assert raw.getvalue() == first


def test_network_has_no_cli_flag():
    import inspect

    from jaff import Network

    assert "_from_cli" not in inspect.signature(Network.__init__).parameters
