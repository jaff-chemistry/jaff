# ABOUTME: Optimizations: exact SymPy rewrites that make generated code cheaper
# ABOUTME: (log derivatives reusing a computed value; more to come, e.g. shared CSE)
"""Rewrites of symbolic expressions that make the generated code cheaper.

Each one is exact: it changes how an expression is written, never its value.
"""

from .log_derivative import log_derivative

__all__ = ["log_derivative"]
