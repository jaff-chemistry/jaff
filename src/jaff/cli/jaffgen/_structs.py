from dataclasses import dataclass, field
from pathlib import Path

from ...config import JAFF_DIR
from .config import JaffgenTomlConfig

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


@dataclass
class State:
    """The loaded config and the files of a single ``jaffgen`` run.

    Attributes
    ----------
    config_file : Path or None
        The ``jaffgen.toml`` the config was loaded from, if any.
    config : JaffgenTomlConfig
        The run's validated config (jaff.toml < jaffgen.toml < command line).
    output_dir : ResolvedPath
        Directory the generated files are written to (defaults to
        :data:`DEFAULT_OUTPUT`).
    output_files : list of ResolvedPath
        Every file to render, gathered from the template, the input
        directory and files, and the config file.
    """

    config_file: Path | None = None
    config: JaffgenTomlConfig = field(default_factory=JaffgenTomlConfig)
    output_dir: ResolvedPath = field(
        default_factory=lambda: ResolvedPath(DEFAULT_OUTPUT, Path())
    )
    output_files: list[ResolvedPath] = field(default_factory=list)
