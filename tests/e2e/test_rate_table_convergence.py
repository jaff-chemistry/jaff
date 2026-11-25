# ABOUTME: Tabulating a full gas-phase network must converge to err_tol within 5 minutes
# ABOUTME: Uses UCLCHEM small gas network; KIDA's 1 K range gaps make refinement unbounded (#217)

import subprocess
import sys
import textwrap

import pytest
from tests._paths import REPO_ROOT

REPO = REPO_ROOT
NETWORK = REPO / "networks" / "uclchem_small_gas" / "uclchem_small_gas_network.jet"
TIMEOUT_SECONDS = 300

# Settings of the original export_table example. Refinement has no grid cap, so
# finishing in time means the table actually reached err_tol.
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
    )
    """
)


def test_rate_table_converges_within_timeout(tmp_path):
    cmd = [sys.executable, "-c", SCRIPT, str(NETWORK), str(tmp_path / "rates.hdf5")]

    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=TIMEOUT_SECONDS
        )
    except subprocess.TimeoutExpired:
        pytest.fail(f"Rate table did not converge within {TIMEOUT_SECONDS} s")

    assert proc.returncode == 0, proc.stdout + proc.stderr
