#![allow(unsafe_op_in_unsafe_fn)]
use pyo3::prelude::*;

mod counting;
mod normalize;
mod parquet_writer;
mod stats;
mod types;

/// A Python module implemented in Rust (bundled as gbcms._rs).
#[pymodule]
fn _rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    pyo3_log::init();
    m.add_function(wrap_pyfunction!(counting::count_bam, m)?)?;
    m.add_function(wrap_pyfunction!(counting::count_bam_binned, m)?)?;
    m.add_function(wrap_pyfunction!(normalize::prepare_variants, m)?)?;
    m.add_function(wrap_pyfunction!(parquet_writer::write_fsd_parquet, m)?)?;
    m.add_class::<types::Variant>()?;
    m.add_class::<types::BaseCounts>()?;
    m.add_class::<normalize::PreparedVariant>()?;
    Ok(())
}
