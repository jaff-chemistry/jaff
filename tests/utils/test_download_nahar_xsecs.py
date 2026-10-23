# ABOUTME: Tests for the NORAD photoionisation downloader/parser (download_nahar_xsecs).
# ABOUTME: Guards the .dat writer against precision loss that merges close energy points.
import logging
from pathlib import Path

import numpy as np

from jaff._utils.download_nahar_xsecs import RY_TO_EV, combine, parse_local

# Minimal FS-format Fe V ground-state block.  The two middle points are 7e-6 Ry
# (~1e-4 eV) apart, as inside the narrow resonances of the real fe5 file.
FE5_RAW = """\
header text ignored by the parser
-----------------------------------------------------------------------
   26   21    2
  0.000000E+00  0.164747E+00
    5    2    0    1
 -0.551002E+01    4
    0.010000
  5.261083E+00 9.236E+00
  8.594656E+00 5.895E+01
  8.594663E+00 4.718E+01
  9.000000E+00 1.000E+00
    0    0    0    0
"""

RAW_RY = np.array([5.261083, 8.594656, 8.594663, 9.0])


def _read_dat(path: Path) -> np.ndarray:
    return np.loadtxt(path, comments="#")


def test_parse_local_keeps_close_energy_points_distinct(tmp_path: Path) -> None:
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    (raw_dir / "fe5.px.gd.fs.txt").write_text(FE5_RAW)

    written = parse_local(raw_dir, tmp_path, logging.getLogger("test"))

    assert written == 1
    data = _read_dat(next(tmp_path.glob("*.dat")))
    assert np.all(np.diff(data[:, 0]) > 0)
    np.testing.assert_allclose(data[:, 0], RAW_RY * RY_TO_EV, rtol=1e-9)
    np.testing.assert_allclose(data[:, 1], [9.236e-18, 5.895e-17, 4.718e-17, 1e-18])


def test_combine_keeps_first_of_duplicate_energies() -> None:
    # NORAD raw files (e.g. fe21) repeat some photon energies with different sigma.
    e = np.array([1.0, 2.0, 2.0, 3.0, 2.0])
    x = np.array([10.0, 20.0, 21.0, 30.0, 22.0])

    energy, xsec = combine([(1.0, e, x)])

    np.testing.assert_array_equal(energy, [1.0, 2.0, 3.0])
    np.testing.assert_array_equal(xsec, [10.0, 20.0, 30.0])


def test_combine_dedupes_each_component_before_averaging() -> None:
    blocks = [
        (1.0, np.array([1.0, 2.0, 2.0]), np.array([10.0, 20.0, 99.0])),
        (3.0, np.array([1.0, 2.0]), np.array([30.0, 40.0])),
    ]

    energy, xsec = combine(blocks)

    np.testing.assert_array_equal(energy, [1.0, 2.0])
    np.testing.assert_allclose(xsec, [(10 + 3 * 30) / 4, (20 + 3 * 40) / 4])
