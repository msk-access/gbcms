//! Cross-source contig (chromosome) name normalization.
//!
//! BAM, GTF, variant files and RNA-editing databases disagree on chromosome
//! naming (`chr1` vs `1`, and the mitochondrion as `M`/`chrM`/`chrMT`/`MT`).
//! Counting joins these sources by chromosome name, so every comparison must pass
//! through one canonical form. Stripping a `chr` prefix alone reconciles
//! `chr1`↔`1` but leaves the mitochondrial aliases unmatched (`chrM` strips to
//! `M`, which still never equals `MT`) — silently dropping every mitochondrial
//! variant/exon/editing join. This normalizer closes that gap.

/// Canonicalize a chromosome name for cross-source comparison.
///
/// - Strips a leading `chr` (any case): `chr1` → `1`, `CHR7` → `7`.
/// - Folds the mitochondrial aliases to `MT`: `M`, `chrM`, `chrMT` → `MT`.
///
/// The returned string is a matching key only — it is never written to output —
/// so the canonical mitochondrial form (`MT`, Ensembl-style) is purely an internal
/// convention and does not change how contigs are reported.
pub fn normalize_contig(name: &str) -> String {
    // Strip a leading "chr" prefix regardless of case (chr/Chr/CHR).
    let bare = match name.get(..3) {
        Some(prefix) if prefix.eq_ignore_ascii_case("chr") => &name[3..],
        _ => name,
    };
    // Fold mitochondrial aliases (M / MT, post-strip) to a single token.
    if bare.eq_ignore_ascii_case("m") || bare.eq_ignore_ascii_case("mt") {
        "MT".to_string()
    } else {
        bare.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strips_chr_prefix() {
        assert_eq!(normalize_contig("chr1"), "1");
        assert_eq!(normalize_contig("CHR7"), "7");
        assert_eq!(normalize_contig("chrX"), "X");
        assert_eq!(normalize_contig("1"), "1"); // already bare — unchanged
    }

    #[test]
    fn reconciles_mitochondrial_aliases() {
        // Every spelling collapses to the same key, so chrM and MT now join.
        for name in ["MT", "M", "chrM", "chrMT", "chrm", "mt"] {
            assert_eq!(normalize_contig(name), "MT", "{name} should normalize to MT");
        }
    }

    #[test]
    fn does_not_overmatch() {
        // A gene/contig that merely starts with "M" is not the mitochondrion.
        assert_eq!(normalize_contig("MT1"), "MT1");
        // "chr" inside the name (not a prefix) is untouched.
        assert_eq!(normalize_contig("GL000220"), "GL000220");
    }
}
