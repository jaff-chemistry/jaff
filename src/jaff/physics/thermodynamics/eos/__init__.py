# ABOUTME: Equation-of-state package: re-exports EosFactory for symbolic
# ABOUTME: internal-energy construction (the EOS config lives in .config)

from .eos_factory import EosFactory

__all__ = ["EosFactory"]
