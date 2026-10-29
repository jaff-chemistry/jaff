# ABOUTME: Rate-segment types for reactions with per-temperature-range rates
# ABOUTME: Re-exports RateSegment, RateSegments and ReactionSource

from ._rate_segment import RateSegment
from ._rate_segments import RateSegments
from ._reaction_source import ReactionSource

__all__ = ["RateSegment", "RateSegments", "ReactionSource"]
