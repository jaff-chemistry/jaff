# ABOUTME: Installs downloaded data files at the paths readers use: HDF5 files are
# ABOUTME: rewritten uncompressed/contiguous, other files copied, owned files tracked.
"""
Install verified downloads as fast-to-read local data files.

Downloads are gzip-compressed to keep transfers small; reading them back costs
decompression on every access.  :func:`install_file` writes an uncompressed,
contiguous copy of each HDF5 download (and a byte copy of other files) to the
path readers open, and records which source it came from so a later download
can replace it -- while leaving files the user supplied untouched.
"""

from __future__ import annotations

import hashlib
import logging
import os
import shutil
from pathlib import Path

import h5py

#: Root attribute marking an HDF5 file installed by :func:`install_file`.
SOURCE_ATTR: str = "jaff_source_sha256"

#: File suffixes treated as HDF5.
HDF5_SUFFIXES: tuple[str, ...] = (".hdf5", ".h5")

_CHUNK: int = 1 << 20


def file_sha256(path: Path) -> str:
    """Return the hex sha256 digest of *path*, read in 1 MiB chunks."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while block := f.read(_CHUNK):
            h.update(block)
    return h.hexdigest()


def _copy_attrs(
    src: h5py.HLObject, dst: h5py.HLObject, skip: tuple[str, ...] = ()
) -> None:
    """Copy every attribute not in *skip*.

    Each keeps its stored dtype (strings, compound, lists).
    """
    for name in src.attrs:
        if name in skip:
            continue
        dst.attrs.create(name, src.attrs[name], dtype=src.attrs.get_id(name).dtype)


def _compressible(src: h5py.Dataset, min_bytes: int | None) -> bool:
    """Whether *src* should be gzip-compressed for threshold *min_bytes*.

    Scalars, null dataspaces, empty arrays and variable-length dtypes (vlen
    strings, ragged arrays) are never compressed: gzip gains little on vlen
    heap pointers and chunking tiny/shape-less data only adds overhead.
    """
    if min_bytes is None or src.shape is None or src.shape == () or src.size == 0:
        return False
    string_info = h5py.check_string_dtype(src.dtype)
    if string_info is not None and string_info.length is None:
        return False
    if h5py.check_vlen_dtype(src.dtype) is not None:
        return False
    return src.nbytes >= min_bytes


def _copy_dataset(
    src: h5py.Dataset, parent: h5py.Group, name: str, min_bytes: int | None, level: int
) -> None:
    """Recreate *src* under *parent*, gzip-compressed if :func:`_compressible`.

    Uncompressed datasets are contiguous unless they have unlimited dimensions,
    which require chunked storage (kept chunked, without a filter).
    """
    unlimited = src.maxshape is not None and any(m is None for m in src.maxshape)
    if src.shape is None:  # null dataspace
        ds = parent.create_dataset(name, data=h5py.Empty(src.dtype))
    elif _compressible(src, min_bytes):
        ds = parent.create_dataset(
            name,
            data=src[()],
            dtype=src.dtype,
            maxshape=src.maxshape,
            chunks=True,
            compression="gzip",
            compression_opts=level,
            shuffle=True,
        )
    elif unlimited:
        ds = parent.create_dataset(
            name, data=src[()], dtype=src.dtype, maxshape=src.maxshape, chunks=True
        )
    else:
        ds = parent.create_dataset(name, data=src[()], dtype=src.dtype)
    _copy_attrs(src, ds)


def _copy_group(
    src: h5py.Group,
    dst: h5py.Group,
    min_bytes: int | None,
    level: int,
    skip_attrs: tuple[str, ...] = (),
) -> None:
    """Recursively copy groups, datasets and attributes from *src* into *dst*."""
    _copy_attrs(src, dst, skip_attrs)
    for name in src:
        link = src.get(name, getlink=True)
        if not isinstance(link, h5py.HardLink):
            raise ValueError(
                f"Unsupported HDF5 link {src.name.rstrip('/')}/{name}: "
                f"{type(link).__name__}"
            )
        obj = src[name]
        if isinstance(obj, h5py.Group):
            _copy_group(obj, dst.create_group(name), min_bytes, level)
        else:
            _copy_dataset(obj, dst, name, min_bytes, level)


def copy_hdf5(
    src: Path,
    dst: Path,
    *,
    compress_min_bytes: int | None = None,
    level: int = 4,
    drop_root_attrs: tuple[str, ...] = (),
    root_attrs: dict[str, str] | None = None,
) -> None:
    """Copy HDF5 *src* to *dst* dataset by dataset, choosing storage per dataset.

    Parameters
    ----------
    src : Path
        Source HDF5 file.
    dst : Path
        Output path; overwritten.
    compress_min_bytes : int or None, optional
        ``None`` (default) writes every dataset uncompressed and contiguous
        (chunked without a filter only if it has unlimited dimensions).
        Otherwise datasets of at least this many bytes are written
        gzip + shuffle compressed and chunked; scalar, empty, null and
        variable-length datasets stay uncompressed regardless.
    level : int, optional
        Gzip level (0-9) for compressed datasets.
    drop_root_attrs : tuple of str, optional
        Root attributes not copied.
    root_attrs : dict of str to str, optional
        Root attributes written after copying (overriding copied ones).

    Raises
    ------
    ValueError
        If *src* contains soft or external links.
    """
    with h5py.File(src, "r") as s, h5py.File(dst, "w") as d:
        _copy_group(s, d, compress_min_bytes, level, drop_root_attrs)
        for key, value in (root_attrs or {}).items():
            d.attrs[key] = value


def decompress_hdf5(src: Path, dst: Path, source_sha256: str) -> None:
    """Write an uncompressed, contiguous copy of HDF5 *src* to *dst*.

    Parameters
    ----------
    src : Path
        Source HDF5 file (typically gzip-compressed and chunked).
    dst : Path
        Output path; overwritten.
    source_sha256 : str
        Hash of *src*, stored as the root attribute :data:`SOURCE_ATTR`.

    Raises
    ------
    ValueError
        If *src* contains soft or external links.
    """
    copy_hdf5(src, dst, root_attrs={SOURCE_ATTR: source_sha256})


def _sidecar(src: Path) -> Path:
    """Marker file recording what was installed from download *src*."""
    return src.with_name(src.name + ".installed-sha256")


def _installed_source(dst: Path, src: Path, is_hdf5: bool) -> str | None:
    """Source hash *dst* was installed from, or ``None`` if not installed by us."""
    if is_hdf5:
        try:
            with h5py.File(dst, "r") as f:
                value = f.attrs.get(SOURCE_ATTR)
        except OSError:  # not a readable HDF5 file
            return None
        if value is None:
            return None
        return value.decode() if isinstance(value, bytes) else str(value)

    sidecar = _sidecar(src)
    if not sidecar.exists():
        return None
    source, output = sidecar.read_text().split()
    return source if file_sha256(dst) == output else None


def _write_atomic(src: Path, dst: Path, source_sha256: str, is_hdf5: bool) -> None:
    """Install *src* to *dst* through a temp file + ``os.replace``."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(f".{dst.name}.tmp-{os.getpid()}")
    try:
        if is_hdf5:
            decompress_hdf5(src, tmp, source_sha256)
        else:
            shutil.copyfile(src, tmp)
        os.replace(tmp, dst)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def install_file(
    src: Path, dst: Path, source_sha256: str, logger: logging.Logger
) -> bool:
    """Install verified download *src* at *dst* unless that would clobber user data.

    Rules: a missing *dst* is installed; a *dst* we installed from the same
    source is left alone; one we installed from another source is replaced; a
    *dst* we did not install (no marker) is kept and a warning is logged.
    HDF5 files carry the marker as the root attribute :data:`SOURCE_ATTR`;
    other files use a ``<src>.installed-sha256`` sidecar next to the download.

    Parameters
    ----------
    src : Path
        Verified download.
    dst : Path
        Install path read by the rest of JAFF.
    source_sha256 : str
        Registry hash of *src*.
    logger : logging.Logger
        Receives the user-supplied-file warning.

    Returns
    -------
    bool
        ``True`` if *dst* was (re)written.
    """
    is_hdf5 = dst.suffix.lower() in HDF5_SUFFIXES
    if dst.exists():
        installed = _installed_source(dst, src, is_hdf5)
        if installed is None:
            logger.warning(
                f"{dst} was not installed by jaff (user-supplied); not overwritten"
            )
            return False
        if installed == source_sha256:
            return False

    _write_atomic(src, dst, source_sha256, is_hdf5)
    if not is_hdf5:
        _sidecar(src).write_text(f"{source_sha256} {file_sha256(dst)}\n")
    return True
