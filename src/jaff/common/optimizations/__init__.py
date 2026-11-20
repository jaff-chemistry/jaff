# ABOUTME: Optimizations: exact SymPy rewrites that make generated code cheaper (log
# ABOUTME: derivatives, shared CSE and logs, guards on branch-only CSE temps, power rules)
"""Rewrites of symbolic expressions that make the generated code cheaper.

Each one is exact: it changes how an expression is written, never its value.
"""

from .cost import cost
from .cse import prune, shared_cse
from .guard import guard
from .log_derivative import log_derivative
from .logs import shared_logs
from .optimize import optimize
from .powers import exp_power, integer_power

__all__ = [
    "cost",
    "exp_power",
    "guard",
    "integer_power",
    "log_derivative",
    "optimize",
    "prune",
    "shared_cse",
    "shared_logs",
]
