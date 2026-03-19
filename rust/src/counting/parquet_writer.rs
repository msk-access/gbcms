//! Parquet export for mFSD fragment size distribution data.
//!
//! Writes per-variant raw fragment size arrays (REF and ALT) to a columnar
//! Parquet file using the native Arrow/Parquet Rust crates — no pyarrow
//! dependency required.
//!
//! ## Schema
//! ```text
//! chrom     : Utf8        (chromosome name)
//! pos       : Int64       (1-based position, MAF/VCF convention)
//! ref       : Utf8        (reference allele)
//! alt       : Utf8        (alternate allele)
//! ref_sizes : List<Int32> (insert sizes of REF-classified fragments, 50–1000 bp)
//! alt_sizes : List<Int32> (insert sizes of ALT-classified fragments, 50–1000 bp)
//! ```
//!
//! ## Usage
//! Called from `pipeline.py` immediately after `count_bam()` when
//! `--mfsd-parquet` is enabled:
//! ```python
//! _get_rs().write_fsd_parquet(str(fsd_path), chroms, positions, refs, alts, counts)
//! ```

use std::fs::File;
use std::sync::Arc;

use arrow_array::{ArrayRef, Int32Array, Int64Array, ListArray, StringArray};
use arrow_buffer::OffsetBuffer;
use arrow_schema::{DataType, Field, Schema};
use parquet::arrow::ArrowWriter;
use parquet::file::properties::WriterProperties;
use pyo3::exceptions::PyIOError;
use pyo3::prelude::*;

use crate::types::BaseCounts;

/// Write per-variant mFSD fragment size arrays to a Parquet file.
///
/// Produces a file at `path` with the schema described in the module-level doc.
/// Each row corresponds to one variant; `ref_sizes` and `alt_sizes` are
/// variable-length lists of fragment insert sizes (50–1000 bp).
///
/// # Arguments
/// * `path`      - Absolute path to the output `.fsd.parquet` file.
/// * `chroms`    - Chromosome names (one per variant).
/// * `positions` - 1-based positions (one per variant).
/// * `refs`      - Reference alleles (one per variant).
/// * `alts`      - Alternate alleles (one per variant).
/// * `counts`    - `BaseCounts` list returned by `count_bam()`.
///                 `ref_sizes` and `alt_sizes` fields are read internally.
///
/// # Errors
/// Returns `PyIOError` if the file cannot be created or written.
///
/// # Notes
/// * `ref_sizes`/`alt_sizes` in `BaseCounts` are NOT exposed to Python via
///   `#[pyo3(get)]`; this function is the only consumer.
/// * Parquet compression uses Zstandard level 1 (fast, always available,
///   better ratio than Snappy for integer lists).
#[pyfunction]
pub fn write_fsd_parquet(
    path: &str,
    chroms: Vec<String>,
    positions: Vec<i64>,
    refs: Vec<String>,
    alts: Vec<String>,
    counts: Vec<BaseCounts>,
) -> PyResult<()> {
    let n = chroms.len();
    debug_assert_eq!(positions.len(), n, "positions length mismatch");
    debug_assert_eq!(refs.len(), n, "refs length mismatch");
    debug_assert_eq!(alts.len(), n, "alts length mismatch");
    debug_assert_eq!(counts.len(), n, "counts length mismatch");

    // ── Build Arrow schema ────────────────────────────────────────────────────
    let schema = Arc::new(Schema::new(vec![
        Field::new("chrom", DataType::Utf8, false),
        Field::new("pos", DataType::Int64, false),
        Field::new("ref", DataType::Utf8, false),
        Field::new("alt", DataType::Utf8, false),
        Field::new(
            "ref_sizes",
            DataType::List(Arc::new(Field::new("item", DataType::Int32, true))),
            true,
        ),
        Field::new(
            "alt_sizes",
            DataType::List(Arc::new(Field::new("item", DataType::Int32, true))),
            true,
        ),
    ]));

    // ── Build flat value + offsets arrays for ListArray ───────────────────────
    // Arrow list columns are stored as (offsets, flat values).
    // offsets[i]..offsets[i+1] is the slice of flat values for row i.

    let mut ref_offsets: Vec<i32> = Vec::with_capacity(n + 1);
    let mut ref_values: Vec<i32> = Vec::new();
    ref_offsets.push(0);
    for c in &counts {
        for &sz in &c.ref_sizes {
            ref_values.push(sz as i32);
        }
        ref_offsets.push(ref_values.len() as i32);
    }

    let mut alt_offsets: Vec<i32> = Vec::with_capacity(n + 1);
    let mut alt_values: Vec<i32> = Vec::new();
    alt_offsets.push(0);
    for c in &counts {
        for &sz in &c.alt_sizes {
            alt_values.push(sz as i32);
        }
        alt_offsets.push(alt_values.len() as i32);
    }

    // ── Construct Arrow arrays ────────────────────────────────────────────────
    let chrom_arr: ArrayRef = Arc::new(StringArray::from(chroms));
    let pos_arr: ArrayRef = Arc::new(Int64Array::from(positions));
    let ref_arr: ArrayRef = Arc::new(StringArray::from(refs));
    let alt_arr: ArrayRef = Arc::new(StringArray::from(alts));

    let ref_sizes_arr: ArrayRef = Arc::new(
        ListArray::try_new(
            Arc::new(Field::new("item", DataType::Int32, true)),
            OffsetBuffer::new(ref_offsets.into()),
            Arc::new(Int32Array::from(ref_values)),
            None, // no null bitmap — all rows are non-null lists
        )
        .map_err(|e| PyIOError::new_err(format!("Arrow list error (ref_sizes): {e}")))?,
    );

    let alt_sizes_arr: ArrayRef = Arc::new(
        ListArray::try_new(
            Arc::new(Field::new("item", DataType::Int32, true)),
            OffsetBuffer::new(alt_offsets.into()),
            Arc::new(Int32Array::from(alt_values)),
            None,
        )
        .map_err(|e| PyIOError::new_err(format!("Arrow list error (alt_sizes): {e}")))?,
    );

    // ── Build RecordBatch ─────────────────────────────────────────────────────
    let batch = arrow_array::RecordBatch::try_new(
        schema.clone(),
        vec![chrom_arr, pos_arr, ref_arr, alt_arr, ref_sizes_arr, alt_sizes_arr],
    )
    .map_err(|e| PyIOError::new_err(format!("RecordBatch construction failed: {e}")))?;

    // ── Write to Parquet ──────────────────────────────────────────────────────
    let file = File::create(path)
        .map_err(|e| PyIOError::new_err(format!("Cannot create mFSD Parquet file '{path}': {e}")))?;

    let props = WriterProperties::builder()
        .set_compression(parquet::basic::Compression::ZSTD(
            parquet::basic::ZstdLevel::try_new(1)
                .map_err(|e| PyIOError::new_err(format!("Invalid ZSTD level: {e}")))?,
        ))
        .build();

    let mut writer = ArrowWriter::try_new(file, schema, Some(props))
        .map_err(|e| PyIOError::new_err(format!("Parquet writer init failed: {e}")))?;

    writer
        .write(&batch)
        .map_err(|e| PyIOError::new_err(format!("Parquet write failed: {e}")))?;

    writer
        .close()
        .map_err(|e| PyIOError::new_err(format!("Parquet writer close failed: {e}")))?;

    log::debug!(
        "write_fsd_parquet: wrote {} variants ({} REF fragments, {} ALT fragments) → {}",
        n,
        counts.iter().map(|c| c.ref_sizes.len()).sum::<usize>(),
        counts.iter().map(|c| c.alt_sizes.len()).sum::<usize>(),
        path,
    );

    Ok(())
}
