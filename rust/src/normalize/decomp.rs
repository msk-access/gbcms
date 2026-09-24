//! Homopolymer decomposition detection.
//!
//! Detects complex variants that may be miscollapsed in a homopolymer run and
//! produces a corrected ALT allele to dual-count against them.

use log::debug;

/// Check if a complex variant looks like a miscollapsed homopolymer event.
///
/// Returns the corrected ALT allele if the pattern is detected, `None` otherwise.
///
/// **Pattern**: The REF allele is a homopolymer run (all same base, ≥3bp),
/// the variant is a net deletion (ref_len > alt_len), and the ALT's last base
/// matches the reference base immediately AFTER the REF span — suggesting
/// D(n) + SNV were merged into a larger deletion than the reads actually support.
///
/// # Corrected allele
/// The run with its last base replaced by the base that follows it, the same
/// length as REF:
/// ```text
/// ref = CCCCCC, alt = T, next_ref_base = T
/// Homopolymer base = C, alt ends with T == next_ref_base, T ≠ C
/// → corrected alt = CCCCCT (REF[..len-1] + T)
/// ```
/// The dual-count keeps whichever of the two alleles more reads support. Reads
/// at real loci carry several forms — the called delins, this allele, a 1bp
/// deletion plus the change (C^(L-2)T), or other alleles of the run — and both
/// classifiers can claim reads of forms neither describes exactly, so the
/// arbitration is a heuristic.
pub(crate) fn check_homopolymer_decomp(
    ref_al: &str,
    alt_al: &str,
    next_ref_base: u8,
) -> Option<String> {
    let ref_bytes = ref_al.as_bytes();
    let alt_bytes = alt_al.as_bytes();

    // Must be a net deletion with a remaining allele
    if ref_bytes.len() <= alt_bytes.len() || alt_bytes.is_empty() {
        return None;
    }

    // REF must be at least 3bp (meaningful homopolymer)
    if ref_bytes.len() < 3 {
        return None;
    }

    // REF must be a homopolymer (all same base)
    let homo_base = ref_bytes[0].to_ascii_uppercase();
    if !ref_bytes
        .iter()
        .all(|&b| b.to_ascii_uppercase() == homo_base)
    {
        return None;
    }

    // ALT's last base must match the next reference base after the run
    let alt_last = alt_bytes.last()?.to_ascii_uppercase();
    if alt_last != next_ref_base.to_ascii_uppercase() {
        return None;
    }

    // ALT's last base must differ from the homopolymer base (it's the "SNV")
    if alt_last == homo_base {
        return None;
    }

    // Construct the corrected ALT allele: every base of the run but the last,
    // then the base that follows the run (REF[..len-1] + alt_last) — the same
    // length as REF.
    let mut corrected = String::with_capacity(ref_bytes.len());
    // All but the last base of the homopolymer
    for &b in &ref_bytes[..ref_bytes.len() - 1] {
        corrected.push(b as char);
    }
    // Replace the last position with the SNV base
    corrected.push(alt_last as char);

    debug!(
        "Homopolymer decomp detected: ref={} alt={} → corrected alt={}",
        ref_al, alt_al, corrected
    );

    Some(corrected)
}
