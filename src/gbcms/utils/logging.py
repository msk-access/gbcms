"""
Logging utilities for gbcms.

Provides centralized logging configuration with dual output:
- Structured logging via Python logging module
- Rich console output for interactive use
"""

import logging
import time
from collections.abc import Callable
from contextlib import contextmanager
from functools import wraps
from typing import Any

from rich.console import Console
from rich.logging import RichHandler

__all__ = [
    "setup_logging",
    "get_logger",
    "timed",
    "log_call",
]

# Module-level console for rich output
_console = Console()


def setup_logging(verbose: bool = False, trace: bool = False, log_file: str | None = None) -> None:
    """
    Configure logging for gbcms.

    Args:
        verbose: If True, set log level to DEBUG. Otherwise INFO.
        trace: If True, enable Rust per-read trace logging (slow).
            Implies verbose. Sets the ``_rs`` Python logger to level 5
            (TRACE) so pyo3-log forwards Rust ``trace!()`` calls through
            the GIL. (pyo3-log names loggers after the extension module's
            lib name ``_rs``, e.g. ``_rs.counting.variant_checks`` — not
            the crate name ``gbcms_rs``.)
        log_file: Optional path to write logs to file.
    """
    if trace:
        verbose = True
    log_level = logging.DEBUG if verbose else logging.INFO

    handlers: list[logging.Handler] = [
        RichHandler(
            console=_console,
            rich_tracebacks=True,
            markup=True,
            show_path=verbose,
        )
    ]

    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(
            logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        )
        handlers.append(file_handler)

    logging.basicConfig(
        level=log_level,
        format="%(message)s",
        datefmt="[%X]",
        handlers=handlers,
        force=True,  # Override existing config
    )

    # Enable Rust trace-level logs when --trace is active.
    # pyo3-log maps Rust trace!() to Python level 5 (below DEBUG=10).
    # By default, pyo3-log caches that level 5 is disabled and
    # short-circuits all trace!() calls in Rust with near-zero GIL cost.
    # Setting the logger to level 5 lets them through (slow but
    # comprehensive per-read diagnostics). The logger is "_rs" — pyo3-log
    # derives logger names from the extension lib name (_rs.counting.…),
    # not the crate name; setting "gbcms_rs" here silently did nothing.
    if trace:
        logging.getLogger("_rs").setLevel(5)  # TRACE level
        # pyo3-log caches each Rust target's effective level at its first
        # record; if any Rust code already logged in this process (library
        # use), the level change above would be silently ignored for those
        # targets. Reset the cache so trace takes effect now.
        try:
            from gbcms._rs import reset_log_caching

            reset_log_caching()
        except ImportError:  # pragma: no cover - extension not built yet
            pass


def get_logger(name: str) -> logging.Logger:
    """Get a logger instance for the given module name."""
    return logging.getLogger(name)


@contextmanager
def timed(operation: str, logger: logging.Logger | None = None):
    """
    Context manager for timing operations.

    Args:
        operation: Description of the operation being timed.
        logger: Logger to use. If None, uses root logger.

    Example:
        with timed("Loading variants", logger):
            variants = load_variants()
    """
    log = logger or logging.getLogger(__name__)
    start = time.perf_counter()
    log.debug("Starting: %s", operation)
    try:
        yield
    finally:
        elapsed = time.perf_counter() - start
        log.debug("Completed: %s (%.3fs)", operation, elapsed)


def log_call(logger: logging.Logger | None = None) -> Callable:
    """
    Decorator to log function calls with timing.

    Args:
        logger: Logger to use. If None, uses function's module logger.

    Example:
        @log_call()
        def process_sample(sample_name: str) -> dict:
            ...
    """

    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            log = logger or logging.getLogger(func.__module__)
            log.debug("Calling %s", func.__name__)
            start = time.perf_counter()
            try:
                result = func(*args, **kwargs)
                elapsed = time.perf_counter() - start
                log.debug("%s completed (%.3fs)", func.__name__, elapsed)
                return result
            except Exception as e:
                log.error("%s failed: %s", func.__name__, e)
                raise

        return wrapper

    return decorator
