# ABOUTME: Directories the tests read from, independent of where a test file sits,
# ABOUTME: so tests can move between folders without rewriting their paths

from pathlib import Path

#: The ``tests/`` directory.
TESTS_DIR: Path = Path(__file__).resolve().parent

#: The repository root (``networks/``, ``examples/``, ``src/``).
REPO_ROOT: Path = TESTS_DIR.parent

#: Small input files: networks, configs, tables.
FIXTURES_DIR: Path = TESTS_DIR / "fixtures"

#: Committed snapshot outputs compared by characterization tests.
SNAPSHOTS_DIR: Path = TESTS_DIR / "snapshots"

#: Committed generated code and numeric references for the codegen tests.
GOLDEN_DIR: Path = TESTS_DIR / "golden"
