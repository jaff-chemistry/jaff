# ABOUTME: ResolvedPath (a file to render and where it lands) and the default
# ABOUTME: jaffgen output directory <repo_root>/generated.
from dataclasses import dataclass
from pathlib import Path

from ...config import JAFF_DIR

#: Fallback output directory (``<repo_root>/generated``) when none is supplied.
DEFAULT_OUTPUT: Path = JAFF_DIR.parent.parent / "generated"


@dataclass(frozen=True)
class ResolvedPath:
    """A path paired with the relative path used to mirror it under ``outdir``.

    Attributes
    ----------
    abspath : Path
        Absolute, resolved path to the source file or directory.
    relpath : Path
        Path used when writing the file under the output directory
        (``outdir / relpath``), so the source tree structure is preserved.
    """

    abspath: Path
    relpath: Path
