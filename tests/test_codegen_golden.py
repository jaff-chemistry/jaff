# ABOUTME: Golden tests for code generation: render Python test templates with jaffgen and
# ABOUTME: evaluate them against committed golden renders on identical random inputs

# The templates in tests/fixtures/templates/ expose every expression generator
# (rates, fluxes, ODEs, RHS, Jacobian variants, dE/dt, auxiliary functions) as a
# Python function.  Each run renders them for every network, loads the fresh
# render and the committed golden render (tests/golden/<network>/) as modules,
# feeds both the same randomly drawn inputs and compares the outputs entry by
# entry.  Inputs are drawn from a fresh seed each session; reproduce a failure
# with JAFF_GOLDEN_SEED=<seed>.  After an intentional codegen change, refresh
# the goldens with JAFF_UPDATE_GOLDEN=1 uv run pytest tests/test_codegen_golden.py
# and review the diff before committing.

import math
import re
from pathlib import Path
from typing import Any, Callable, Dict, List, Tuple

import pytest

from jaff.core.network._spec import NetworkSpec
from tests.codegen_render import (
    GOLDEN,
    NETWORKS,
    TEMPLATES,
    UPDATE,
    draw_inputs,
    evaluate,
    free_names,
    load,
    template_files,
    template_functions,
)

REL_TOL = 1e-9
ABS_TOL = 1e-300


def _cases() -> List[Any]:
    return [
        pytest.param(name, path.name, func, id=f"{name}-{path.stem}-{func}")
        for name in NETWORKS
        for path in template_files(name)
        for func in template_functions(path)
    ]


def _format_rows(rows: List[Tuple[Any, float, float]], limit: int = 20) -> str:
    lines = [f"  {k!r}: generated={a!r} golden={e!r}" for k, a, e in rows[:limit]]
    if len(rows) > limit:
        lines.append(f"  ... {len(rows) - limit} more")
    return "\n".join(lines)


def _row_scales(*tables: dict) -> Callable[[Any], float]:
    """Largest |value| per row (first key element) across *tables*.

    Scalar keys (1-D outputs) share a single scale: the largest |value| overall.
    """
    rows: Dict[Any, float] = {}
    for table in tables:
        for key, value in table.items():
            if not math.isfinite(value):
                continue
            row = key[0] if isinstance(key, tuple) else None
            rows[row] = max(rows.get(row, 0.0), abs(value))
    return lambda key: rows.get(key[0] if isinstance(key, tuple) else None, 0.0)


@pytest.mark.parametrize("name", list(NETWORKS))
def test_golden_file_set_matches(name: str, rendered: Callable[[str], Path]) -> None:
    generated = {p.name for p in rendered(name).glob("*.py")}
    golden = {p.name for p in (GOLDEN / name).glob("*.py")}
    assert generated == golden, f"no/stale golden for {name}; run JAFF_UPDATE_GOLDEN=1"


@pytest.mark.parametrize("name", list(NETWORKS))
def test_aux_template_covers_every_aux_function(name: str) -> None:
    spec = NetworkSpec(
        NETWORKS[name],
        config=None,
        errors=False,
        label=None,
        funcfile=True,
        duplicate_policy=None,
        expand_nuclei=True,
        _from_cli=False,
        _metadata={},
    )
    template = (TEMPLATES / f"aux_{name}.py").read_text()
    covered = set(re.findall(r"GET aux_func FOR (\w+)", template))
    assert covered == set(spec.aux_funcs), (
        f"aux_{name}.py is out of sync with the network's .jfunc: "
        f"missing {sorted(set(spec.aux_funcs) - covered)}, "
        f"unknown {sorted(covered - set(spec.aux_funcs))}"
    )


@pytest.mark.parametrize("name, filename, func", _cases())
def test_generated_matches_golden(
    name: str,
    filename: str,
    func: str,
    rendered: Callable[[str], Path],
    render_seed: int,
) -> None:
    if UPDATE:
        pytest.skip(f"golden refreshed: {GOLDEN / name}")

    gen_path = rendered(name) / filename
    gold_path = GOLDEN / name / filename
    assert gold_path.exists(), f"no golden {gold_path}; run JAFF_UPDATE_GOLDEN=1"

    free = free_names(gen_path.read_text())
    gold_free = free_names(gold_path.read_text())
    assert free == gold_free, (
        f"free names differ: generated-only {sorted(free.keys() - gold_free.keys())}, "
        f"golden-only {sorted(gold_free.keys() - free.keys())}, changed "
        f"{sorted(k for k in free.keys() & gold_free.keys() if free[k] != gold_free[k])}"
    )
    inputs = draw_inputs(free, render_seed)

    tag = f"{name}_{Path(filename).stem}"
    actual = evaluate(load(gen_path, f"gen_{tag}"), func, inputs, render_seed)
    expected = evaluate(load(gold_path, f"gold_{tag}"), func, inputs, render_seed)

    where = f"{name}/{filename}::{func}() (seed={render_seed})"

    non_finite = [
        (k, actual.get(k), v) for k, v in expected.items() if not math.isfinite(v)
    ]
    assert not non_finite, (
        f"{where}: golden is non-finite at these inputs; adjust the sampling "
        f"ranges in tests/codegen_render.py\n{_format_rows(non_finite)}"
    )

    # Entries missing on one side count as 0 (codegen omits exact zeros); each
    # entry gets an absolute tolerance scaled to its row (whole vector for 1-D
    # outputs), so platform-dependent round-off terms that are negligible within
    # their row do not flip the comparison, while tiny but genuine entries in a
    # row of tiny values are still compared strictly.
    keys = sorted(actual.keys() | expected.keys())
    scale = _row_scales(expected, actual)
    mismatched = [
        (k, actual.get(k, 0.0), expected.get(k, 0.0))
        for k in keys
        if actual.get(k, 0.0)
        != pytest.approx(
            expected.get(k, 0.0), rel=REL_TOL, abs=max(ABS_TOL, REL_TOL * scale(k))
        )
    ]
    assert not mismatched, (
        f"{where}: {len(mismatched)}/{len(expected)} entries differ "
        f"(rel={REL_TOL})\n{_format_rows(mismatched)}"
    )
