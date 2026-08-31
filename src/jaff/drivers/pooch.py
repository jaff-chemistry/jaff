import logging
from pathlib import Path

import pooch
from rich.filesize import decimal
from rich.progress import TaskID

from ..config import DATA_DIR
from ..io._logger import jaff_progress

pooch.get_logger().setLevel(logging.WARNING)


class _JaffProgressBar:
    """tqdm-compatible adapter so pooch renders downloads on ``jaff_progress``.

    pooch's :class:`~pooch.HTTPDownloader` drives any progress object through a
    tqdm-like protocol: it assigns ``.total`` (bytes) before the transfer,
    calls ``.update(n)`` per chunk, then ``.reset()`` and ``.close()`` at the
    end.  This maps those calls onto a task on the shared JAFF Rich bar.

    The shared bar's ``MofNCompleteColumn`` renders task counts as plain
    integers, which would expose raw byte totals.  To keep the bar readable the
    task is driven by completion *percent* (0-100) while the human-readable
    transferred / total size is shown in the description via
    :func:`rich.filesize.decimal`.
    """

    def __init__(self, description: str = "Downloading") -> None:
        self.description = description
        self.total: int | None = None
        self._downloaded: int = 0
        self._task_id: TaskID | None = None

    def _label(self) -> str:
        if self.total:
            return f"{self.description} ({decimal(self._downloaded)} / {decimal(self.total)})"
        return f"{self.description} ({decimal(self._downloaded)})"

    def _percent(self) -> float:
        if not self.total:
            return 0.0
        return min(100.0, 100.0 * self._downloaded / self.total)

    def _ensure_task(self) -> TaskID:
        if self._task_id is None:
            self._task_id = jaff_progress.add_task(self._label(), total=100)
        return self._task_id

    def update(self, n: int) -> None:
        task_id = self._ensure_task()
        self._downloaded += n
        jaff_progress.update(
            task_id, completed=self._percent(), description=self._label()
        )
        jaff_progress.refresh()

    def reset(self) -> None:
        if self._task_id is not None:
            jaff_progress.update(self._task_id, completed=100, description=self._label())
            jaff_progress.refresh()

    def close(self) -> None:
        if self._task_id is not None:
            jaff_progress.remove_task(self._task_id)
            self._task_id = None
            self._downloaded = 0


class Pooch:
    """Cached wrapper around a :class:`pooch.Pooch` data fetcher.

    Instances are deduplicated by ``base_url`` + ``cache_path``: constructing a
    :class:`Pooch` with the same pair returns the previously created instance
    from the class-level ``_registry``, avoiding redundant fetcher objects for
    the same remote source.

    The registry of downloadable files (names, hashes, URLs) is downloaded from
    ``registry.txt`` under ``base_url`` (fetched unverified, since its own hash
    is not known ahead of time).
    """

    _registry: dict[str, "Pooch"] = {}

    def __new__(cls, base_url: str, cache_path: Path) -> "Pooch":
        """Return the cached instance for ``base_url``/``cache_path`` if any.

        Builds a registry key from the two arguments. If the key is already
        present, the stored instance is returned; otherwise a new instance is
        created, registered under the key, and returned.
        """
        key = f"_{base_url}__{cache_path}"
        if key in cls._registry:
            return cls._registry[key]

        instance = super().__new__(cls)
        cls._registry[key] = instance

        return instance

    def __init__(self, base_url: str, cache_path: Path) -> None:
        """Create the underlying pooch fetcher and load the file registry.

        Parameters
        ----------
        base_url : str
            Root URL the registered files are downloaded from.
        cache_path : Path
            Local directory the fetched files are cached in.

        Notes
        -----
        ``__new__`` returns cached instances, but Python still re-invokes
        ``__init__`` on each construction. The ``_initialized`` guard makes
        repeat calls a no-op so the fetcher and registry are built only once.
        """
        if getattr(self, "_initialized", False):
            return
        self._initialized = True

        self.pooch: pooch.Pooch = pooch.create(
            path=cache_path,
            base_url=base_url,
            registry=None,
        )
        registry_path = pooch.retrieve(
            url=f"{base_url}/registry.txt",
            known_hash=None,
            fname="registry.txt",
            path=cache_path,
        )
        self.pooch.load_registry(registry_path)

    def fetch_file(self, filename: str) -> None:
        """Download ``filename`` from the registry, rendering a progress bar.

        The file is fetched into the cache directory (skipped if already
        present and hash-valid) and progress is shown on the shared JAFF Rich
        bar via :class:`_JaffProgressBar`.
        """
        self.pooch.fetch(
            filename,
            progressbar=_JaffProgressBar(f"Downloading {filename}"),  # ty: ignore[invalid-argument-type]
        )


def download_xsecs() -> None:
    """Fetch the photochemistry cross-section data files into ``data/xsecs``.

    Downloads the Leiden, NORAD, and Verner cross-section files from the ANU
    mirror, caching them under the package ``data/xsecs`` directory. Files
    already present and hash-valid are not re-downloaded.
    """
    pooch = Pooch(
        "https://www.mso.anu.edu.au/~anishs",
        DATA_DIR,
    )
    for file in ["xsecs/leiden.hdf5", "xsecs/norad.hdf5", "xsecs/verner_1996.csv"]:
        pooch.fetch_file(file)


def download_shielding() -> None:
    """Fetch the line-shielding data files into ``data/shielding``.

    Downloads the collapsed Leiden line-shielding HDF5 file from the ANU
    mirror, caching it under the package ``data/shielding`` directory. Files
    already present and hash-valid are not re-downloaded.
    """
    pooch = Pooch(
        "https://www.mso.anu.edu.au/~anishs",
        DATA_DIR,
    )

    for file in ["shielding/leiden.hdf5"]:
        pooch.fetch_file(file)


def download_background_radiation() -> None:
    """Fetch the background-radiation data file into ``data/background_radiation``.

    Downloads the collapsed background-radiation HDF5 file from the ANU mirror,
    caching it under the package ``data/background_radiation`` directory. Files
    already present and hash-valid are not re-downloaded.
    """
    pooch = Pooch(
        "https://www.mso.anu.edu.au/~anishs",
        DATA_DIR,
    )

    for file in ["background_radiation/radiation.hdf5"]:
        pooch.fetch_file(file)
