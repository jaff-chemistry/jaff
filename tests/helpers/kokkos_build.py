# ABOUTME: Shared helpers for tests of generated kokkos_ode code: configure the CMake
# ABOUTME: project (fetching header-only-integrators) and locate the library's headers
"""Configure a generated ``kokkos_ode`` folder and find header-only-integrators.

CMake fetches the library from GitHub, or uses a local copy named by the
environment variable ``FETCHCONTENT_SOURCE_DIR_HEADER_ONLY_INTEGRATORS``.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

LOCAL_SOURCE = "FETCHCONTENT_SOURCE_DIR_HEADER_ONLY_INTEGRATORS"


def configure(folder: Path) -> subprocess.CompletedProcess:
    """Configure *folder*'s CMake project into ``build/`` (Release)."""
    command = ["cmake", "-S", ".", "-B", "build", "-DCMAKE_BUILD_TYPE=Release"]
    if os.environ.get(LOCAL_SOURCE):
        command.append(f"-D{LOCAL_SOURCE}={os.environ[LOCAL_SOURCE]}")

    return subprocess.run(command, cwd=folder, capture_output=True, text=True)


def include_dir(folder: Path) -> Path | None:
    """The header-only-integrators ``include`` folder for a configured *folder*."""
    fetched = next(folder.glob("build/_deps/header_only_integrators-src/include"), None)
    if fetched is None and os.environ.get(LOCAL_SOURCE):
        return Path(os.environ[LOCAL_SOURCE]) / "include"

    return fetched
