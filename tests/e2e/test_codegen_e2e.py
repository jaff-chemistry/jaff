# ABOUTME: End-to-end smoke test: every built-in jaffgen template renders the GOW network
# ABOUTME: Numeric golden comparisons live in test_codegen_golden.py

from types import SimpleNamespace

import pytest

from jaff.cli import JaffGen
from tests._paths import FIXTURES_DIR, REPO_ROOT


def _jaffgen(network_path, template, outdir):
    """Drive the ``jaffgen`` engine in-process (same path as the CLI).

    All CLI options default to ``None`` except the ones supplied here, matching
    ``jaffgen --network <net> --template <name> --outdir <dir>``.
    """
    args = SimpleNamespace(
        network=str(network_path),
        config=str(FIXTURE_CONFIG),
        label=None,
        funcfile=None,
        duplicate_policy=None,
        expand_nuclei=None,
        strict=None,
        network_config=None,
        outdir=str(outdir),
        indir=None,
        files=None,
        template=template,
    )
    JaffGen(args)


REPO = REPO_ROOT
FIXTURE_CONFIG = FIXTURES_DIR / "jaffgen.toml"
GOW = REPO / "networks" / "GOW" / "GOW.jet"
# Built-in template collections; each file names its language ($$JAFF LANG or extension).
ALL_TEMPLATES = ["python_solve_ivp", "microphysics", "fortran_dlsodes", "kokkos_ode"]


# --------------------------------------------------------------------------- #
# GOW generates for every template                                           #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("template", ALL_TEMPLATES)
def test_gow_generates_for_all_templates(template, tmp_path):
    """The GOW network builds without error and emits non-empty files."""
    _jaffgen(GOW, template=template, outdir=tmp_path)
    files = [p for p in tmp_path.iterdir() if p.is_file()]
    assert files, f"template {template} produced no files"
    assert all(p.stat().st_size > 0 for p in files), f"empty output file for {template}"
