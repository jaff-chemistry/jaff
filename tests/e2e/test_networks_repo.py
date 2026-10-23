# ABOUTME: Opt-in tests for round-tripping repo-provided networks via JSON
# ABOUTME: Runs only when JAFF_TEST_REPO_NETWORKS=1 is set (can be slow)

import io
import os
import tempfile
from contextlib import redirect_stderr, redirect_stdout

import pytest
import sympy

import jaff.core.network as jn
from jaff import Network
from tests._paths import REPO_ROOT


@pytest.mark.skipif(
    os.environ.get("JAFF_TEST_REPO_NETWORKS") != "1",
    reason="Set JAFF_TEST_REPO_NETWORKS=1 to run (slow).",
)
def test_repo_networks_roundtrip_json():
    # Keep test output clean / fast.
    jn.tqdm = lambda x: x

    repo_root = str(REPO_ROOT)
    networks_dir = os.path.join(repo_root, "networks")

    network_files = []
    for name in sorted(os.listdir(networks_dir)):
        if name.endswith("_functions"):
            continue
        path = os.path.join(networks_dir, name)
        if os.path.isfile(path):
            network_files.append(path)

    loaded = []
    unparsable = []
    unserializable = []

    for path in network_files:
        try:
            with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                net = Network(path)
            loaded.append((path, net))
        except Exception as e:
            unparsable.append((path, type(e).__name__, str(e)))

    # By default, only validate the networks the current parser can load.
    # If you want this test to be strict about *all* files in `networks/`,
    # set JAFF_TEST_REPO_NETWORKS_STRICT=1.
    if unparsable and os.environ.get("JAFF_TEST_REPO_NETWORKS_STRICT") == "1":
        msg = ["Some repo network files could not be parsed:"]
        for p, t, s in unparsable:
            msg.append(f"- {p}: {t}: {s.splitlines()[0]}")
        pytest.fail("\n".join(msg))

    if not loaded:
        pytest.skip("No repo network files could be parsed; nothing to test.")

    for path, net in loaded:
        # Treat networks containing undefined SymPy functions as unserializable.
        has_undef = False
        for r in net.reactions:
            if isinstance(r.rate_coefficient, sympy.Basic):
                if any(
                    type(f.func) is sympy.core.function.UndefinedFunction
                    for f in r.rate_coefficient.atoms(sympy.Function)
                ):
                    has_undef = True
                    break
            if isinstance(r.dE, sympy.Basic):
                if any(
                    type(f.func) is sympy.core.function.UndefinedFunction
                    for f in r.dE.atoms(sympy.Function)
                ):
                    has_undef = True
                    break
        if has_undef:
            unserializable.append(path)
            continue

        fd, json_path = tempfile.mkstemp(suffix=".jaff")
        os.close(fd)
        try:
            net.to_jaff_file(json_path)
            net2 = Network(json_path)
        finally:
            os.unlink(json_path)

        assert [sp.name for sp in net2.species] == [sp.name for sp in net.species]
        assert len(net2.reactions) == len(net.reactions)

        for r1, r2 in zip(net.reactions, net2.reactions):
            assert r2.get_verbatim() == r1.get_verbatim()
            assert r2.tmin == r1.tmin
            assert r2.tmax == r1.tmax

            if isinstance(r1.rate_coefficient, sympy.Basic):
                assert isinstance(r2.rate_coefficient, sympy.Basic)
                assert r2.rate_coefficient == r1.rate_coefficient
            else:
                assert r2.rate_coefficient == r1.rate_coefficient

            assert r2.dE == r1.dE

    if not unserializable:
        return
    # If you want this test to be strict about undefined functions, set
    # JAFF_TEST_REPO_NETWORKS_STRICT=1.
    if os.environ.get("JAFF_TEST_REPO_NETWORKS_STRICT") == "1":
        msg = [
            "Some repo networks are unserializable due to undefined SymPy function calls:"
        ]
        for p in unserializable:
            msg.append(f"- {p}")
        pytest.fail("\n".join(msg))
