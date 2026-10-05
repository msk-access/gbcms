"""Atomic output files.

Every output gbcms writes goes first to a temp file beside it, named
``.<name>.partial``, which is flushed, fsynced, closed (with its errors checked)
and then renamed over the final name. A run that fails mid-write therefore
leaves nothing at the final path and removes its temp file: no truncated output
can pass for a finished one. The rename is atomic on POSIX (the temp file is in
the same directory, so on the same filesystem).

An output path that is a symlink is written through: the temp file sits beside
the link's target, which is replaced, so the link stays a link. A replaced file
keeps its permission bits. A path that exists but is not a regular file (a
device such as ``/dev/stdout``, a FIFO) is written in place: it cannot be
renamed over.
"""

from __future__ import annotations

import contextlib
import os
import shutil
from collections.abc import Iterator
from pathlib import Path
from typing import IO, Any

__all__ = [
    "atomic_output",
    "atomic_path",
    "commit_partial",
    "discard_partial",
    "partial_path",
    "write_target",
    "writes_in_place",
]


def _target(path: Path) -> Path:
    """The file an output path names: a symlink's target, else the path."""
    return Path(os.path.realpath(path))


def writes_in_place(path: Path) -> bool:
    """Whether ``path`` is written directly: it exists and is not a regular file."""
    target = _target(path)
    return target.exists() and not target.is_file()


def partial_path(path: Path) -> Path:
    """The temp file an output is written to before its rename."""
    target = _target(path)
    return target.with_name(f".{target.name}.partial")


def write_target(path: Path) -> Path:
    """Where a writer opens ``path``: its temp file, or ``path`` itself when it is
    written in place."""
    return Path(path) if writes_in_place(path) else partial_path(path)


def commit_partial(path: Path) -> None:
    """Rename a finished temp file over the final name (fsynced first; the
    replaced file's permission bits kept). No-op when written in place."""
    if writes_in_place(path):
        return
    tmp, target = partial_path(path), _target(path)
    with open(tmp, "rb") as fh:
        os.fsync(fh.fileno())
    if target.exists():
        shutil.copymode(target, tmp)
    os.replace(tmp, target)


def discard_partial(path: Path) -> None:
    """Remove a temp file after a failure; the final path is not touched."""
    if writes_in_place(path):
        return
    with contextlib.suppress(FileNotFoundError):
        partial_path(path).unlink()


@contextlib.contextmanager
def atomic_output(path: Path, mode: str = "w", **open_kwargs: Any) -> Iterator[IO[Any]]:
    """Open ``path`` for writing through its temp file; renamed on success,
    removed on any exception (which is re-raised)."""
    path = Path(path)
    try:
        with open(write_target(path), mode, **open_kwargs) as fh:
            yield fh
        commit_partial(path)
    except BaseException:
        discard_partial(path)
        raise


@contextlib.contextmanager
def atomic_path(path: Path) -> Iterator[Path]:
    """Yield the path a writer that takes a path (the native parquet writers,
    polars) should write; renamed over ``path`` on success, removed on failure."""
    path = Path(path)
    try:
        yield write_target(path)
        commit_partial(path)
    except BaseException:
        discard_partial(path)
        raise
