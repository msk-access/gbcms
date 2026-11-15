#![allow(unsafe_op_in_unsafe_fn)]
use pyo3::prelude::*;

mod annotation;
mod counting;
mod normalize;
mod shared;
mod types;

/// pyo3-log's reset handle, kept so Python can invalidate the per-target
/// level cache after changing logger levels (see `reset_log_caching`).
static LOG_RESET_HANDLE: std::sync::OnceLock<pyo3_log::ResetHandle> = std::sync::OnceLock::new();

/// Invalidate pyo3-log's cached per-target log levels.
///
/// pyo3-log caches each Rust target's effective Python level at its FIRST
/// record; a later `logging.getLogger("_rs").setLevel(...)` on the Python
/// side is otherwise silently ignored for targets that already logged.
/// `setup_logging` calls this after changing levels so trace can be enabled
/// mid-process (library use); the CLI path works either way because it
/// configures logging before any counting.
#[pyfunction]
fn reset_log_caching() {
    if let Some(handle) = LOG_RESET_HANDLE.get() {
        handle.reset();
        // A counting call in progress (it releases the GIL) may hold the maximum at the old
        // levels (`CountingLogCap`); open it, so the new levels decide from now on.
        log::set_max_level(log::LevelFilter::Trace);
    }
}

/// The most verbose level any Python logger asks for: the root logger and every logger
/// created so far. A Rust record goes to the logger named after its target (pyo3-log), which
/// is either one of these or a child of one, so no logger can pass a level below this one.
fn python_lowest_log_level(py: Python<'_>) -> PyResult<i64> {
    let logging = py.import("logging")?;
    let root = logging.getattr("root")?;
    let mut lowest: i64 = root.call_method0("getEffectiveLevel")?.extract()?;
    let logger_class = logging.getattr("Logger")?;
    let table = root.getattr("manager")?.getattr("loggerDict")?;
    for logger in table.call_method0("values")?.try_iter()? {
        let logger = logger?;
        if logger.is_instance(&logger_class)? {
            let level: i64 = logger.call_method0("getEffectiveLevel")?.extract()?;
            lowest = lowest.min(level);
        }
    }
    Ok(lowest)
}

/// The Rust level filter that lets through every record a Python logger at `level` would pass
/// (pyo3-log gives Rust's levels the Python levels 40, 30, 20, 10 and 5).
fn level_filter_for(level: i64) -> log::LevelFilter {
    match level {
        l if l > 40 => log::LevelFilter::Off,
        l if l > 30 => log::LevelFilter::Error,
        l if l > 20 => log::LevelFilter::Warn,
        l if l > 10 => log::LevelFilter::Info,
        l if l > 5 => log::LevelFilter::Debug,
        _ => log::LevelFilter::Trace,
    }
}

/// For the length of one counting call, the Rust `log` crate's global maximum follows the
/// Python loggers' levels, so a `trace!` or `debug!` that Python would drop is not even built
/// (its arguments cost time on every read). A record the Python side would pass is never
/// dropped: the cap is the most verbose level any logger asks for when the call starts, and a
/// failure to read the levels leaves it open. On drop the maximum goes back to `Trace`, the
/// filter installed with the module, so code outside a counting call logs as before; two
/// calls at once can only open the cap early, never close it. `reset_log_caching` opens it at
/// once, so a level lowered during a call (from another thread) is followed after the reset.
pub(crate) struct CountingLogCap {
    active: bool,
}

impl CountingLogCap {
    pub(crate) fn new(py: Python<'_>) -> Self {
        if LOG_RESET_HANDLE.get().is_none() {
            // another logger owns the process: leave its level alone
            return Self { active: false };
        }
        let filter = python_lowest_log_level(py)
            .map(level_filter_for)
            .unwrap_or(log::LevelFilter::Trace);
        log::set_max_level(filter);
        Self { active: true }
    }
}

impl Drop for CountingLogCap {
    fn drop(&mut self) {
        if self.active {
            log::set_max_level(log::LevelFilter::Trace);
        }
    }
}

/// The commit this build came from (8 hex characters), or empty when the build
/// had neither git nor `GBCMS_BUILD_COMMIT` (see build.rs). Output provenance
/// names it after the version, so two builds of one version are told apart.
#[pyfunction]
fn build_commit() -> &'static str {
    env!("GBCMS_BUILD_COMMIT")
}

#[cfg(test)]
mod log_cap_tests {
    use super::level_filter_for;
    use log::LevelFilter;

    #[test]
    fn a_python_level_maps_to_the_filter_that_passes_what_it_passes() {
        // pyo3-log gives error, warn, info, debug and trace the Python levels 40, 30, 20, 10, 5;
        // a logger at a level between two of them passes the less verbose one only
        assert_eq!(level_filter_for(50), LevelFilter::Off);
        assert_eq!(level_filter_for(41), LevelFilter::Off);
        assert_eq!(level_filter_for(40), LevelFilter::Error);
        assert_eq!(level_filter_for(30), LevelFilter::Warn);
        assert_eq!(level_filter_for(25), LevelFilter::Warn);
        assert_eq!(level_filter_for(20), LevelFilter::Info);
        assert_eq!(level_filter_for(15), LevelFilter::Info);
        assert_eq!(level_filter_for(10), LevelFilter::Debug);
        assert_eq!(level_filter_for(6), LevelFilter::Debug);
        assert_eq!(level_filter_for(5), LevelFilter::Trace);
        assert_eq!(level_filter_for(0), LevelFilter::Trace);
    }
}

/// A Python module implemented in Rust (bundled as gbcms._rs).
// gil_used: the module keeps the GIL on a free-threaded interpreter, as before
// pyo3 0.28 changed the default; gbcms is neither built nor tested free-threaded.
#[pymodule(gil_used = true)]
fn _rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Forward Rust log records to Python's `logging` module. The filter must
    // admit TRACE: pyo3_log::init() caps the global max level at DEBUG, which
    // silently discarded every per-read trace!() diagnostic behind the CLI's
    // --trace flag. With the filter open, each record still passes a cached
    // Python-side level check (Caching::LoggersAndLevels), so trace stays
    // near-free unless setup_logging(--trace) drops the `_rs` logger to
    // level 5 (pyo3-log names loggers after the extension lib: _rs.counting…).
    match pyo3_log::Logger::default()
        .filter(log::LevelFilter::Trace)
        .install()
    {
        Ok(handle) => {
            let _ = LOG_RESET_HANDLE.set(handle);
        }
        Err(_) => {
            // A logger is already installed (module re-import in the same
            // process) — keep it rather than abort the import.
            log::debug!("pyo3-log already installed; keeping existing logger");
        }
    }
    m.add_function(wrap_pyfunction!(reset_log_caching, m)?)?;
    m.add_function(wrap_pyfunction!(build_commit, m)?)?;
    m.add_function(wrap_pyfunction!(counting::count_bam_binned, m)?)?;
    m.add_function(wrap_pyfunction!(counting::count_bam_binned_observations, m)?)?;
    m.add_function(wrap_pyfunction!(normalize::prepare_variants, m)?)?;
    m.add_function(wrap_pyfunction!(counting::write_fsd_parquet, m)?)?;
    m.add_function(wrap_pyfunction!(shared::stats::fisher_exact_2x2_py, m)?)?;
    m.add_class::<types::Variant>()?;
    m.add_class::<types::BaseCounts>()?;
    m.add_class::<types::Observation>()?;
    m.add_class::<normalize::PreparedVariant>()?;
    Ok(())
}
