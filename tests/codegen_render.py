# ABOUTME: Shared helpers for numeric codegen tests: render Python test templates with
# ABOUTME: jaffgen, load renders as modules, and draw reproducible random inputs for them

import ast
import builtins
import importlib.util
import math
import os
import re
import shutil
import zlib
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import Any, Callable, Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

from tests.jaffgen_run import run_jaffgen
from tests._paths import FIXTURES_DIR, GOLDEN_DIR, REPO_ROOT

REPO = REPO_ROOT
TEMPLATES = FIXTURES_DIR / "templates"
GOLDEN = GOLDEN_DIR
CONFIG = FIXTURES_DIR / "jaffgen.toml"
UPDATE = os.environ.get("JAFF_UPDATE_GOLDEN") == "1"

NETWORKS: Dict[str, Path] = {
    "h_photo": REPO / "networks" / "h_photoionization" / "h_photo.jet",
    "GOW": REPO / "networks" / "GOW" / "GOW.jet",
}

# Log-uniform sampling ranges.  Anything not listed uses the defaults.  GOW's
# rendered code raises at 10-12 K: its own condition `(T < 10) | (fit < 1e-100)`
# evaluates a fit that overflows there.  (Below ~90 K it raised until split_log:
# exp(-63590/T) underflowed to 0 inside a log(); above ~3500 K until CSE temporaries
# were guarded: a 10**(polynomial) fit was hoisted out of its Piecewise branch.)
SCALAR_RANGES: Dict[str, Tuple[float, float]] = {"tgas": (15.0, 1.0e5)}
DEFAULT_SCALAR_RANGE = (0.5, 2.0)
DEFAULT_VECTOR_RANGE = (1.0e-2, 1.0e4)

_PARTIAL = re.compile(r"^(\w+)_partial_(\d+)$")


def template_files(name: str) -> List[Path]:
    return [TEMPLATES / "codegen_outputs.py", TEMPLATES / f"aux_{name}.py"]


def template_functions(path: Path) -> List[str]:
    return re.findall(r"^def (\w+)\(\):", path.read_text(), flags=re.MULTILINE)


# --------------------------------------------------------------------------- #
# Rendering                                                                    #
# --------------------------------------------------------------------------- #


def render(network: Path, files: Sequence[Path], outdir: Path) -> Path:
    """Render *files* for *network* in-process; return the run folder in *outdir*."""
    gen = run_jaffgen(
        SimpleNamespace(
            network=str(network),
            config=str(CONFIG),
            label=None,
            funcfile=None,
            duplicate_policy=None,
            expand_nuclei=None,
            strict=None,
            network_config=None,
            outdir=str(outdir),
            indir=None,
            files=[str(f) for f in files],
            template=None,
        )
    )

    return gen.output_dir


def refresh_golden(name: str, rendered: Path) -> None:
    gdir = GOLDEN / name
    gdir.mkdir(parents=True, exist_ok=True)
    for old in gdir.glob("*.py"):
        old.unlink()
    for new in rendered.glob("*.py"):
        shutil.copy(new, gdir / new.name)


def load(path: Path, tag: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(f"_jaff_render_{tag}", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# --------------------------------------------------------------------------- #
# Inputs                                                                       #
# --------------------------------------------------------------------------- #


def _constant_index(name: str, node: ast.Subscript) -> Optional[int]:
    index = node.slice
    if isinstance(index, ast.Tuple):
        raise ValueError(
            f"'{name}' is indexed as a matrix ({ast.unparse(node)}); "
            "generated state vectors must be indexed 1-D"
        )
    if isinstance(index, ast.Constant) and isinstance(index.value, int):
        return index.value
    return None


def free_names(source: str) -> Dict[str, Tuple[str, int]]:
    """Classify every unbound name in *source* as scalar, vector or function.

    Returns ``{name: (kind, size)}`` where *size* is the vector length implied
    by the largest constant index (``0`` for scalars and functions).
    """
    tree = ast.parse(source)
    bound: Set[str] = set(dir(builtins))
    loads: Set[str] = set()
    calls: Set[str] = set()
    subscripts: Dict[str, List[ast.Subscript]] = {}

    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            bound.add(node.name)
        elif isinstance(node, (ast.Import, ast.ImportFrom)):
            bound.update((a.asname or a.name).split(".")[0] for a in node.names)
        elif isinstance(node, ast.arg):
            bound.add(node.arg)
        elif isinstance(node, ast.Name):
            (loads if isinstance(node.ctx, ast.Load) else bound).add(node.id)
        elif isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            calls.add(node.func.id)
        elif isinstance(node, ast.Subscript) and isinstance(node.value, ast.Name):
            subscripts.setdefault(node.value.id, []).append(node)

    free: Dict[str, Tuple[str, int]] = {}
    for name in sorted(loads - bound):
        if name in calls and name in subscripts:
            raise ValueError(f"'{name}' is used both as a function and a vector")
        if name in calls:
            free[name] = ("function", 0)
        elif name in subscripts:
            indices = [_constant_index(name, node) for node in subscripts[name]]
            if None in indices:
                raise ValueError(f"vector '{name}' is indexed by a non-constant")
            free[name] = ("vector", max(indices) + 1)  # type: ignore[type-var]
        else:
            free[name] = ("scalar", 0)
    return free


def _name_rng(seed: int, name: str) -> np.random.Generator:
    """Per-name stream: adding a name never shifts the values of the others."""
    return np.random.default_rng([seed, zlib.crc32(name.encode())])


def _log_uniform(
    rng: np.random.Generator, lo: float, hi: float, size: Optional[int] = None
) -> Any:
    return 10.0 ** rng.uniform(math.log10(lo), math.log10(hi), size=size)


def _stub(rng: np.random.Generator, partial: Optional[int]) -> Callable[..., float]:
    """Smooth positive stand-in for an interpolation function, or its exact partial.

    ``f(a) = scale * (1 + sin(S) / 2)`` with ``S = sum_k w_k * log1p(|a_k|)``, so
    ``f_partial_k`` stubs are true derivatives and finite differences agree.
    """
    scale = float(_log_uniform(rng, 0.1, 10.0))
    weights = [float(w) for w in rng.uniform(0.1, 1.0, size=16)]

    def phase(args: Tuple[float, ...]) -> float:
        return sum(w * math.log1p(abs(float(a))) for w, a in zip(weights, args))

    if partial is None:
        return lambda *args: scale * (1.0 + 0.5 * math.sin(phase(args)))

    def derivative(*args: float) -> float:
        a = float(args[partial])
        slope = weights[partial] * math.copysign(1.0, a) / (1.0 + abs(a))
        return scale * 0.5 * math.cos(phase(args)) * slope

    return derivative


def draw_inputs(free: Dict[str, Tuple[str, int]], seed: int) -> Dict[str, Any]:
    """Draw one value per free name; vectors are 1-D tuples of floats."""
    values: Dict[str, Any] = {}
    for name, (kind, size) in free.items():
        if kind == "function":
            match = _PARTIAL.match(name)
            base, partial = (match[1], int(match[2])) if match else (name, None)
            values[name] = _stub(_name_rng(seed, base), partial)
            continue

        rng = _name_rng(seed, name)
        if kind == "vector":
            vec = _log_uniform(rng, *DEFAULT_VECTOR_RANGE, size=size)
            values[name] = tuple(float(v) for v in vec)
        else:
            lo, hi = SCALAR_RANGES.get(name, DEFAULT_SCALAR_RANGE)
            values[name] = float(_log_uniform(rng, lo, hi))
    return values


def evaluate(module: ModuleType, func: str, inputs: Dict[str, Any], seed: int) -> Dict:
    """Call ``module.func()`` with *inputs* injected as module globals."""
    vars(module).update(inputs)
    try:
        return {k: float(v) for k, v in getattr(module, func)().items()}
    except Exception as exc:
        exc.add_note(f"while evaluating {module.__file__}::{func}() (seed={seed})")
        raise
