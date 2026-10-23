# ABOUTME: parsedListProps: the per-reaction dict NetworkParser hands to Network
# ABOUTME: Typed keys for reactants, products, bounds, rate, type and provenance

from typing import NotRequired, TypedDict

parsedListProps = TypedDict(
    "parsedListProps",
    {
        "r": list[str],
        "p": list[str],
        "tmin": float | None,
        "tmax": float | None,
        "rate": str,
        "type": str,
        "string": str,
        "source_index": NotRequired[int],
        "format": NotRequired[str],
        "metadata": NotRequired[dict],
    },
)
