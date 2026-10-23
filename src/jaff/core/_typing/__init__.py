# ABOUTME: Typed dicts shared by the core package
# ABOUTME: Re-exports ElementProps and ReactionProps

from ._elements import ElementProps
from ._reaction import ReactionProps

__all__ = [
    "ElementProps",
    "ReactionProps",
]
