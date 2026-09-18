#![allow(unsafe_op_in_unsafe_fn)]
use pyo3::prelude::*;

#[allow(dead_code)] // annotation types used internally by engine.rs (not exported to Python)
mod annotation;
mod counting;
mod normalize;
mod shared;
mod types;

/// A Python module implemented in Rust (bundled as gbcms._rs).
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
    }
}

#[pymodule]
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
    #[cfg(feature = "legacy-parity")]
    m.add_function(wrap_pyfunction!(counting::count_bam, m)?)?;
    m.add_function(wrap_pyfunction!(counting::count_bam_binned, m)?)?;
    m.add_function(wrap_pyfunction!(counting::count_bam_binned_observations, m)?)?;
    m.add_function(wrap_pyfunction!(counting::build_gtf_cache, m)?)?;
    m.add_function(wrap_pyfunction!(normalize::prepare_variants, m)?)?;
    m.add_function(wrap_pyfunction!(counting::write_fsd_parquet, m)?)?;
    m.add_function(wrap_pyfunction!(shared::stats::fisher_exact_2x2_py, m)?)?;
    m.add_class::<types::Variant>()?;
    m.add_class::<types::BaseCounts>()?;
    m.add_class::<types::Observation>()?;
    m.add_class::<normalize::PreparedVariant>()?;
    Ok(())
}
