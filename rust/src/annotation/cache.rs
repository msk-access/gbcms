//! M5a: on-disk cache for the parsed GTF intermediate.
//!
//! Parsing a full Ensembl GTF is ~8.7s (3.4M lines of text). Under a Nextflow
//! cohort that cost is paid once *per task* — N samples × 8.7s of redundant
//! parsing. This module persists the **parsed intermediate** (exon records,
//! splice sites, introns, chrom map) to a small binary file so repeated runs over
//! the same GTF skip the text parse and only rebuild the cheap COITrees.
//!
//! What is **not** cached: the [`coitrees::COITree`] structures. Their in-memory
//! layout is architecture-specific (x86-AVX vs ARM-Neon vs scalar fallback, chosen
//! at compile time), so a serialized tree would be unsafe to share across a
//! heterogeneous cluster. Rebuilding them from the cached exon records is cheap and
//! portable — see [`super::build_exon_trees`].
//!
//! The cache is strictly an optimization: every read/write error falls back to a
//! plain parse (logged, never silent, never fatal).

use std::collections::hash_map::DefaultHasher;
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};
use std::path::{Path, PathBuf};

use log::{debug, info, warn};
use serde::{Deserialize, Serialize};

use super::gtf::parse_gtf_to_bundle;
use super::{build_exon_trees, AnnotationIndex, ExonRecord, TranscriptIntrons};

/// Bumped whenever the bundle layout or any parse semantics change, so an older
/// cache file is rejected (re-parsed) instead of deserializing into stale data.
/// Also folded into the cache key, so a bump changes every file name too.
pub(crate) const CACHE_FORMAT_VERSION: u32 = 1;

/// The serializable parsed intermediate — everything an [`AnnotationIndex`] needs
/// except the (arch-specific, rebuilt-on-load) COITrees.
#[derive(Serialize, Deserialize)]
pub(crate) struct GtfIndexBundle {
    /// Format tag; validated against [`CACHE_FORMAT_VERSION`] on load.
    pub(crate) format_version: u32,
    pub(crate) exons: Vec<ExonRecord>,
    pub(crate) splice_sites: HashMap<u32, Vec<i32>>,
    pub(crate) transcript_introns: HashMap<String, TranscriptIntrons>,
    pub(crate) chrom_map: HashMap<String, u32>,
}

impl GtfIndexBundle {
    /// Rebuild the COITrees from the exon records and assemble the full index.
    pub(crate) fn into_index(self) -> AnnotationIndex {
        let exon_trees = build_exon_trees(&self.exons);
        AnnotationIndex::new(
            exon_trees,
            self.exons,
            self.splice_sites,
            self.transcript_introns,
            self.chrom_map,
        )
    }
}

/// Parse a GTF, using `cache_dir` as a disk cache for the parsed intermediate.
///
/// On a cache hit the ~8.7s text parse is skipped entirely (only the cheap tree
/// rebuild runs). On a miss the GTF is parsed and the result is written back for
/// the next process. The cache is keyed on the GTF identity (path + size + mtime)
/// and the variant-chromosome set (the parse is variant-guided), so a cohort that
/// genotypes the same variant list across N samples reuses one cache entry.
///
/// Any cache error (unreadable, corrupt, stale-version, unwritable dir) is logged
/// and falls back to a plain parse — the cache never affects correctness.
pub(crate) fn parse_gtf_cached(
    gtf_path: &str,
    variant_chroms: &HashSet<String>,
    cache_dir: &str,
) -> anyhow::Result<AnnotationIndex> {
    let cache_file: Option<PathBuf> = cache_key(gtf_path, variant_chroms)
        .map(|key| Path::new(cache_dir).join(format!("gbcms-gtf-{key}.idx")));

    if let Some(ref cf) = cache_file {
        if let Some(bundle) = try_load(cf) {
            info!(
                "GTF cache hit: reused parsed index from {} ({} exons) — skipped GTF text parse",
                cf.display(),
                bundle.exons.len(),
            );
            return Ok(bundle.into_index());
        }
    }

    // Miss (or caching disabled by an un-stattable GTF): do the full parse.
    let bundle = parse_gtf_to_bundle(gtf_path, variant_chroms)?;
    if let Some(ref cf) = cache_file {
        try_write(cf, &bundle);
    }
    Ok(bundle.into_index())
}

/// Build a stable cache key from the GTF identity + the variant-chrom set.
///
/// Returns `None` if the GTF cannot be stat-ed (caller then parses without
/// caching). Uses size + mtime rather than a content hash so we never re-read the
/// multi-GB file just to decide freshness; an in-place edit that preserves both is
/// vanishingly unlikely and would be caught by [`CACHE_FORMAT_VERSION`] only on a
/// layout change — acceptable for an optimization-only cache.
fn cache_key(gtf_path: &str, variant_chroms: &HashSet<String>) -> Option<String> {
    let meta = std::fs::metadata(gtf_path).ok()?;
    let mtime_ns = meta
        .modified()
        .ok()
        .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
        .map(|d| d.as_nanos())
        .unwrap_or(0);

    let mut hasher = DefaultHasher::new();
    CACHE_FORMAT_VERSION.hash(&mut hasher);
    // Canonicalize so equivalent paths (./gtf vs /abs/gtf) share a key when possible.
    let canon = std::fs::canonicalize(gtf_path)
        .map(|p| p.to_string_lossy().into_owned())
        .unwrap_or_else(|_| gtf_path.to_string());
    canon.hash(&mut hasher);
    meta.len().hash(&mut hasher);
    mtime_ns.hash(&mut hasher);
    // Variant chroms gate which exons are loaded → part of the identity. Sort for
    // order-independence (HashSet iteration order is non-deterministic).
    let mut chroms: Vec<&String> = variant_chroms.iter().collect();
    chroms.sort();
    for c in chroms {
        c.hash(&mut hasher);
    }
    Some(format!("{:016x}", hasher.finish()))
}

/// Try to load + validate a cached bundle. Any failure → `None` (parse instead).
fn try_load(cache_file: &Path) -> Option<GtfIndexBundle> {
    let bytes = match std::fs::read(cache_file) {
        Ok(b) => b,
        // Missing file is the normal first-run case → debug, not warn.
        Err(e) => {
            debug!("GTF cache miss ({}): {}", cache_file.display(), e);
            return None;
        }
    };
    match bincode::deserialize::<GtfIndexBundle>(&bytes) {
        Ok(bundle) if bundle.format_version == CACHE_FORMAT_VERSION => Some(bundle),
        Ok(bundle) => {
            warn!(
                "Ignoring GTF cache {}: format version {} != expected {} (will re-parse)",
                cache_file.display(),
                bundle.format_version,
                CACHE_FORMAT_VERSION,
            );
            None
        }
        Err(e) => {
            warn!(
                "Ignoring corrupt GTF cache {} ({}); will re-parse",
                cache_file.display(),
                e,
            );
            None
        }
    }
}

/// Serialize + atomically write the bundle. Best-effort: any error is warned and
/// swallowed (the run proceeds on the freshly parsed index).
fn try_write(cache_file: &Path, bundle: &GtfIndexBundle) {
    if let Some(parent) = cache_file.parent() {
        if let Err(e) = std::fs::create_dir_all(parent) {
            warn!(
                "Could not create GTF cache dir {}: {} (continuing without caching)",
                parent.display(),
                e,
            );
            return;
        }
    }
    let bytes = match bincode::serialize(bundle) {
        Ok(b) => b,
        Err(e) => {
            warn!("Could not serialize GTF cache (continuing without caching): {}", e);
            return;
        }
    };
    // Atomic publish: write a pid-unique temp file then rename into place, so
    // concurrent writers (Nextflow fan-out) never observe a partial file.
    let tmp = cache_file.with_extension(format!("idx.tmp-{}", std::process::id()));
    if let Err(e) = std::fs::write(&tmp, &bytes) {
        warn!("Could not write GTF cache temp {}: {}", tmp.display(), e);
        return;
    }
    if let Err(e) = std::fs::rename(&tmp, cache_file) {
        warn!("Could not publish GTF cache {}: {}", cache_file.display(), e);
        let _ = std::fs::remove_file(&tmp);
        return;
    }
    info!(
        "Wrote GTF index cache: {} ({} exons, {} bytes)",
        cache_file.display(),
        bundle.exons.len(),
        bytes.len(),
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tiny_bundle() -> GtfIndexBundle {
        let mut chrom_map = HashMap::new();
        chrom_map.insert("1".to_string(), 0u32);
        GtfIndexBundle {
            format_version: CACHE_FORMAT_VERSION,
            exons: vec![ExonRecord {
                transcript_id: "ENST1".to_string(),
                gene_id: "ENSG1".to_string(),
                chrom_id: 0,
                start: 100,
                end: 200,
                strand: '+',
            }],
            splice_sites: HashMap::from([(0u32, vec![100, 200])]),
            transcript_introns: HashMap::new(),
            chrom_map,
        }
    }

    fn tmp_path(tag: &str) -> PathBuf {
        std::env::temp_dir().join(format!("gbcms-m5a-{}-{}.idx", tag, std::process::id()))
    }

    #[test]
    fn write_then_load_roundtrips_and_rebuilds_trees() {
        let cf = tmp_path("roundtrip");
        try_write(&cf, &tiny_bundle());
        let loaded = try_load(&cf).expect("freshly written cache must load");
        assert_eq!(loaded.exons.len(), 1);
        assert_eq!(loaded.exons[0].start, 100);
        // into_index must rebuild the COITrees from the cached exon records.
        let idx = loaded.into_index();
        assert_eq!(idx.n_exons(), 1);
        assert_eq!(idx.n_chromosomes(), 1);
        let _ = std::fs::remove_file(&cf);
    }

    #[test]
    fn load_missing_file_is_none() {
        // Normal first-run case: no cache yet → fall back to parse.
        assert!(try_load(&std::env::temp_dir().join("gbcms-m5a-does-not-exist.idx")).is_none());
    }

    #[test]
    fn load_corrupt_file_is_none() {
        let cf = tmp_path("corrupt");
        std::fs::write(&cf, b"not a valid bincode payload").unwrap();
        assert!(try_load(&cf).is_none(), "corrupt cache must fall back, not panic");
        let _ = std::fs::remove_file(&cf);
    }

    #[test]
    fn load_wrong_format_version_is_none() {
        let cf = tmp_path("version");
        let mut b = tiny_bundle();
        b.format_version = CACHE_FORMAT_VERSION + 1; // a future, incompatible layout
        std::fs::write(&cf, bincode::serialize(&b).unwrap()).unwrap();
        assert!(try_load(&cf).is_none(), "version mismatch must invalidate the cache");
        let _ = std::fs::remove_file(&cf);
    }

    #[test]
    fn cache_key_is_chrom_order_independent() {
        // The key folds in variant chroms; HashSet order must not change it.
        // (Uses this source file as a stand-in existing path to stat.)
        let path = file!();
        let a = HashSet::from(["1".to_string(), "2".to_string(), "X".to_string()]);
        let b = HashSet::from(["X".to_string(), "1".to_string(), "2".to_string()]);
        assert_eq!(cache_key(path, &a), cache_key(path, &b));
        let c = HashSet::from(["1".to_string()]);
        assert_ne!(cache_key(path, &a), cache_key(path, &c));
    }
}
