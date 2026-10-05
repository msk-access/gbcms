"""Atomic output files.

Every output gbcms writes goes first to a temp file beside it, named
``.<name>.partial``, which is flushed, fsynced, closed (with its errors checked)
and then renamed over the final name. A run that fails mid-write therefore
leaves nothing at the final path and removes its temp file: no truncated output
can pass for a finished one. The rename is atomic on POSIX (the temp file is in
the same directory, so on the same filesystem).
"""

from __future__ import annotations

import contextlib
import os
from collections.abc import Iterator
from pathlib import Path
from typing import IO, Any

__all__ = ["atomic_output", "atomic_path", "partial_path", "commit_partial", "discard_partial"]


def partial_path(path: Path) -> Path:
    """The temp file an output is written to before its rename."""
    path = Path(path)
    return path.with_name(f".{path.name}.partial")


def _fsync_file(path: Path) -> None:
    with open(path, "rb") as fh:
        os.fsync(fh.fileno())


def commit_partial(path: Path) -> None:
    """Rename a finished temp file over the final name (fsynced first)."""
    tmp = partial_path(path)
    _fsync_file(tmp)
    os.replace(tmp, path)


def discard_partial(path: Path) -> None:
    """Remove a temp file after a failure; the final path is not touched."""
    with contextlib.suppress(FileNotFoundError):
        partial_path(path).unlink()


@contextlib.contextmanager
def atomic_output(path: Path, mode: str = "w", **open_kwargs: Any) -> Iterator[IO[Any]]:
    """Open ``path`` for writing through its temp file; renamed on success,
    removed on any exception (which is re-raised)."""
    path = Path(path)
    tmp = partial_path(path)
    try:
        with open(tmp, mode, **open_kwargs) as fh:
            yield fh
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
    except BaseException:
        discard_partial(path)
        raise


@contextlib.contextmanager
def atomic_path(path: Path) -> Iterator[Path]:
    """Yield the temp path for a writer that takes a path (the native parquet
    writers, polars); renamed over ``path`` on success, removed on failure."""
    path = Path(path)
    try:
        yield partial_path(path)
        commit_partial(path)
    except BaseException:
        discard_partial(path)
        raise
