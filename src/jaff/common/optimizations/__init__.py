# ABOUTME: Optimizations: exact SymPy rewrites that make generated code cheaper
# ABOUTME: (log derivatives reusing a computed value, CSE shared with symbol definitions)
"""Rewrites of symbolic expressions that make the generated code cheaper.

Each one is exact: it changes how an expression is written, never its value.
"""

from .cost import cost
from .cse import prune, shared_cse
from .log_derivative import log_derivative

__all__ = ["cost", "log_derivative", "prune", "shared_cse"]
