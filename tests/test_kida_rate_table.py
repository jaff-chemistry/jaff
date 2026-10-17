# ABOUTME: Tabulating the full KIDA network must converge to err_tol within 5 minutes
# ABOUTME: Currently fails: 1 K gaps between rate ranges need ~64k grid points (#217)

import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
KIDA = REPO / "networks" / "kida_uva_2024" / "gas_reactions_kida.uva.2024.jet"
TIMEOUT_SECONDS = 300

# Settings of the original export_table example. nT_max=None removes the grid cap,
# so finishing in time means the table actually reached err_tol (with the cap,
# refinement would stop early with only a warning).
SCRIPT = textwrap.dedent(
    """
    import sys
    from jaff import Network

    net = Network(sys.argv[1])
    net.write_table(
        sys.argv[2],
        T_min=5.0,
        T_max=5000.0,
        nT=128,
        err_tol=0.1,
        rate_min=1e-40,
        rate_max=1e-5,
        fast_log=True,
        format="hdf5",
        include_all=False,
        nT_max=None,
    )
    """
)


def test_kida_rate_table_converges_within_timeout(tmp_path):
    cmd = [sys.executable, "-c", SCRIPT, str(KIDA), str(tmp_path / "kida_rates.hdf5")]

    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=TIMEOUT_SECONDS
        )
    except subprocess.TimeoutExpired:
        pytest.fail(f"KIDA rate table did not converge within {TIMEOUT_SECONDS} s")

    assert proc.returncode == 0, proc.stdout + proc.stderr
