# ABOUTME: Produces gzip-compressed copies of locally edited data files for upload
# ABOUTME: to the download mirror, printing/recording their registry.txt lines.
"""Compress local data files for upload to the JAFF download mirror.

JAFF installs data files uncompressed (see :mod:`jaff.drivers._install`), so a
locally regenerated file is uncompressed too.  The mirror serves compressed
files, hash-checked against its ``registry.txt``.  This script writes, for
each FILE under ROOT, a copy at ``OUTDIR/<relpath>``:

- ``.hdf5`` / ``.h5``: every dataset of at least ``--min-size`` bytes is gzip
  (+ shuffle) compressed; scalar, empty and variable-length datasets are left
  uncompressed.  The install marker attribute ``jaff_source_sha256`` is
  dropped.  Unless ``--no-verify``, the output is re-read and compared with the
  source (every group, dataset, attribute, dtype and shape).
- anything else: byte copy (sha256-checked unless ``--no-verify``).

One registry line ``<relpath> <sha256>`` per file is printed to stdout, and
with ``--registry`` merged into that file (existing lines for the same relpath
replaced in place, new ones appended).

Usage
-----
::

    python -m jaff._utils.compress_hdf5 src/jaff/data/xsecs/leiden.hdf5 \\
        src/jaff/data/dust/dust.hdf5 --outdir upload/ --registry upload/registry.txt

then upload the contents of ``upload/`` to the mirror.

Exit status
-----------
A FILE outside ROOT or an invalid option exits via :class:`SystemExit`
(argparse error, status 2).  A verification mismatch raises
:class:`ValueError`.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

import h5py
import numpy as np

from jaff.config import DATA_DIR
from jaff.drivers._install import HDF5_SUFFIXES, SOURCE_ATTR, copy_hdf5, file_sha256
from jaff.io import JaffLogger

#: Default dataset size (bytes) from which HDF5 datasets are compressed.
DEFAULT_MIN_SIZE: int = 65536
#: Default gzip level.
DEFAULT_LEVEL: int = 4


def _attrs_equal(a: h5py.HLObject, b: h5py.HLObject, skip: tuple[str, ...]) -> bool:
    """Whether *a* and *b* carry the same attributes (ignoring *skip* on *a*)."""
    if set(a.attrs) - set(skip) != set(b.attrs):
        return False
    for k in b.attrs:
        if a.attrs.get_id(k).dtype != b.attrs.get_id(k).dtype:
            return False
        if not np.array_equal(np.asarray(a.attrs[k]), np.asarray(b.attrs[k])):
            return False
    return True


def _datasets_equal(a: h5py.Dataset, b: h5py.Dataset) -> bool:
    """Whether *a* and *b* have the same dtype, shape and values."""
    if a.dtype != b.dtype or a.shape != b.shape:
        return False
    if a.shape is None:  # null dataspace: nothing to compare
        return True
    return bool(np.array_equal(np.asarray(a[()]), np.asarray(b[()])))


def verify_hdf5(src: Path, dst: Path, skip_root_attrs: tuple[str, ...] = ()) -> None:
    """Check *dst* holds exactly the groups, datasets and attributes of *src*.

    Parameters
    ----------
    src, dst : Path
        Files to compare.
    skip_root_attrs : tuple of str, optional
        Root attributes of *src* expected to be absent from *dst*.

    Raises
    ------
    ValueError
        Naming the first object that differs.
    """
    with h5py.File(src, "r") as s, h5py.File(dst, "r") as d:
        if not _attrs_equal(s, d, skip_root_attrs):
            raise ValueError(f"{dst}: root attributes differ from {src}")
        names: List[str] = []
        s.visit(names.append)
        other: List[str] = []
        d.visit(other.append)
        if set(names) != set(other):
            raise ValueError(f"{dst}: object names differ from {src}")
        for name in names:
            so, do = s[name], d[name]
            if type(so) is not type(do) or not _attrs_equal(so, do, ()):
                raise ValueError(f"{dst}: {name} differs from {src}")
            if isinstance(so, h5py.Dataset) and not _datasets_equal(so, do):
                raise ValueError(f"{dst}: dataset {name} differs from {src}")


def compress_file(
    src: Path, dst: Path, *, min_size: int, level: int, verify: bool
) -> str:
    """Write a compressed (HDF5) or byte-identical (other) copy of *src* at *dst*.

    The copy goes through a temp file next to *dst* and ``os.replace``; the
    temp file is removed on failure.

    Parameters
    ----------
    src : Path
        Input file.
    dst : Path
        Output path (parent directories created).
    min_size : int
        Bytes from which HDF5 datasets are compressed.
    level : int
        Gzip level (0-9).
    verify : bool
        Re-read the output and compare it with *src*.

    Returns
    -------
    str
        sha256 of the written file.

    Raises
    ------
    ValueError
        On verification mismatch, or if *src* has soft/external links.
    """
    is_hdf5 = src.suffix.lower() in HDF5_SUFFIXES
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(f".{dst.name}.tmp-{os.getpid()}")
    try:
        if is_hdf5:
            copy_hdf5(
                src,
                tmp,
                compress_min_bytes=min_size,
                level=level,
                drop_root_attrs=(SOURCE_ATTR,),
            )
            if verify:
                verify_hdf5(src, tmp, skip_root_attrs=(SOURCE_ATTR,))
        else:
            shutil.copyfile(src, tmp)
            if verify and file_sha256(src) != file_sha256(tmp):
                raise ValueError(f"{dst}: copy of {src} differs (sha256)")
        os.replace(tmp, dst)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    return file_sha256(dst)


def update_registry(path: Path, entries: Dict[str, str]) -> None:
    """Merge ``<relpath> <sha256>`` *entries* into registry file *path*.

    A line whose first token is a relpath in *entries* is replaced in place;
    remaining entries are appended in order.  Other lines keep their order.
    The file is created if missing and written atomically.

    Parameters
    ----------
    path : Path
        Registry file.
    entries : dict of str to str
        Relpath (posix) to sha256.
    """
    lines = path.read_text().splitlines() if path.exists() else []
    pending = dict(entries)
    out: List[str] = []
    for line in lines:
        tokens = line.split()
        if tokens and tokens[0] in pending:
            key = tokens[0]
            out.append(f"{key} {pending.pop(key)}")
        else:
            out.append(line)
    out.extend(f"{k} {v}" for k, v in pending.items())

    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        tmp.write_text("".join(f"{line}\n" for line in out))
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def _level(value: str) -> int:
    """argparse type for ``--level``: an int in 0..9."""
    try:
        n = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"invalid level {value!r} (expected 0-9)")
    if not 0 <= n <= 9:
        raise argparse.ArgumentTypeError(f"invalid level {n} (expected 0-9)")
    return n


def _parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="python -m jaff._utils.compress_hdf5",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("files", nargs="+", type=Path, metavar="FILE")
    ap.add_argument("--outdir", required=True, type=Path, help="output root")
    ap.add_argument(
        "--root", type=Path, default=DATA_DIR, help="data root (default: %(default)s)"
    )
    ap.add_argument(
        "--min-size",
        type=int,
        default=DEFAULT_MIN_SIZE,
        help="compress HDF5 datasets of at least this many bytes (default: %(default)s)",
    )
    ap.add_argument(
        "--level", type=_level, default=DEFAULT_LEVEL, help="gzip level 0-9 (default: 4)"
    )
    ap.add_argument("--registry", type=Path, help="registry.txt to create/update")
    ap.add_argument("--no-verify", action="store_true", help="skip output check")
    return ap


def main(argv: Optional[List[str]] = None) -> None:
    """Command-line entry point; see the module docstring."""
    ap = _parser()
    args = ap.parse_args(argv)
    root = args.root.resolve()

    jobs: List[tuple[Path, str]] = []
    for f in args.files:
        src = f.resolve()
        if not src.is_relative_to(root):
            ap.error(f"{f} is not under --root {root}")
        jobs.append((src, src.relative_to(root).as_posix()))

    logger = JaffLogger().get_logger()
    entries: Dict[str, str] = {}
    for src, rel in jobs:
        dst = args.outdir / rel
        t0 = time.perf_counter()
        sha = compress_file(
            src, dst, min_size=args.min_size, level=args.level, verify=not args.no_verify
        )
        entries[rel] = sha
        # sys.stdout.write, not print: tests (and callers) may patch print.
        sys.stdout.write(f"{rel} {sha}\n")
        sys.stdout.flush()
        logger.info(
            f"{rel}: {src.stat().st_size / 1e6:.1f} MB -> "
            f"{dst.stat().st_size / 1e6:.1f} MB ({time.perf_counter() - t0:.1f} s)"
        )

    if args.registry is not None:
        update_registry(args.registry, entries)
        logger.info(f"updated {args.registry}")


if __name__ == "__main__":
    main()
