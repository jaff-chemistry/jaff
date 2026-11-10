# ABOUTME: Optimizations: SymPy rewrites that make generated code cheaper while
# ABOUTME: staying exact (e.g. derivatives that reuse an already computed value)
"""Rewrites of symbolic expressions that make the generated code cheaper.

Each one is exact: it changes how an expression is written, never its value.
"""

from .derivative import derivative

__all__ = ["derivative"]
