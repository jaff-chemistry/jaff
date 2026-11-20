# ABOUTME: Optimizations: exact SymPy rewrites that make generated code cheaper (log
# ABOUTME: derivatives, CSE shared with symbol definitions, guards on branch-only CSE temps)
"""Rewrites of symbolic expressions that make the generated code cheaper.

Each one is exact: it changes how an expression is written, never its value.
"""

from .cost import cost
from .cse import prune, shared_cse
from .guard import guard
from .log_derivative import log_derivative

__all__ = ["cost", "guard", "log_derivative", "prune", "shared_cse"]
