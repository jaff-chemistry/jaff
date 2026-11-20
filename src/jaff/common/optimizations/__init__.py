# ABOUTME: Optimizations: exact SymPy rewrites that make generated code cheaper (log
# ABOUTME: derivatives, shared CSE, guards on branch-only CSE temps, integer powers)
"""Rewrites of symbolic expressions that make the generated code cheaper.

Each one is exact: it changes how an expression is written, never its value.
"""

from .cost import cost
from .cse import prune, shared_cse
from .guard import guard
from .log_derivative import log_derivative
from .powers import integer_powers

__all__ = ["cost", "guard", "integer_powers", "log_derivative", "prune", "shared_cse"]
