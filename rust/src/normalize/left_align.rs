//! bcftools-style left-alignment algorithm.

/// Left-align a variant using the bcftools `realign_left()` algorithm.
///
/// Three-phase algorithm:
///
/// **Phase 1 — Extend left:**
/// Prepend up to `aln_win` upstream reference bases to both alleles,
/// giving Phase 2 the full repeat context for proper right-trimming.
///
/// **Phase 2 — Right-trim:**
/// While the rightmost base is identical in ref and alt, and both alleles
/// have length > 1 (preserving VCF anchor), remove the rightmost base.
///
/// **Phase 3 — Left-trim:**
/// While the leftmost base is identical and both alleles have length > 1,
/// remove the leftmost base and advance `pos`.
///
/// # Arguments
/// * `pos` — 0-based genomic position of the variant
/// * `ref_allele` — Reference allele bytes (e.g., b"ACTCT")
/// * `alt_allele` — Alternate allele bytes (e.g., b"A")
/// * `reference` — Reference sequence covering the region
/// * `ref_offset` — 0-based genomic start position of `reference`
/// * `aln_win` — Maximum number of bases to extend left (bcftools default: 100)
///
/// # Returns
/// `(new_pos, new_ref, new_alt, was_modified)` — normalized coordinates
pub fn left_align_variant(
    pos: i64,
    ref_allele: &[u8],
    alt_allele: &[u8],
    reference: &[u8],
    ref_offset: i64,
    aln_win: usize,
) -> (i64, Vec<u8>, Vec<u8>, bool) {
    let original_pos = pos;
    let original_ref = ref_allele.to_vec();
    let original_alt = alt_allele.to_vec();

    // SNPs and same-length substitutions: cannot be left-aligned
    if ref_allele.len() == 1 && alt_allele.len() == 1 {
        return (pos, original_ref, original_alt, false);
    }

    // Phase 1: Extend both alleles to the left by up to aln_win bases.
    // This gives Phase 2 the full repeat context to right-trim through.
    let pos_in_ref = (pos - ref_offset) as usize;
    let actual_extend = aln_win.min(pos_in_ref);

    let mut r = Vec::with_capacity(actual_extend + ref_allele.len());
    let mut a = Vec::with_capacity(actual_extend + alt_allele.len());

    // Prepend upstream reference bases to both alleles
    let prefix_start = pos_in_ref - actual_extend;
    for &base in reference.iter().take(pos_in_ref).skip(prefix_start) {
        r.push(base);
        a.push(base);
    }
    let mut cur_pos = pos - actual_extend as i64;

    // Append original alleles
    r.extend_from_slice(ref_allele);
    a.extend_from_slice(alt_allele);

    // Phase 2: Right-trim — remove matching trailing bases
    // Keep at least 1 base in each allele (VCF anchor requirement)
    while r.len() > 1 && a.len() > 1
        && r.last().unwrap().eq_ignore_ascii_case(a.last().unwrap())
    {
        r.pop();
        a.pop();
    }

    // Phase 3: Left-trim — remove matching leading bases
    // Keep at least 1 base in each allele
    while r.len() > 1 && a.len() > 1
        && r[0].eq_ignore_ascii_case(&a[0])
    {
        r.remove(0);
        a.remove(0);
        cur_pos += 1;
    }

    // Case-insensitive comparison: the trim phases use eq_ignore_ascii_case,
    // so the modified check must also be case-insensitive to avoid false positives
    // when FASTA has soft-masked (lowercase) bases near the variant.
    let was_modified = cur_pos != original_pos
        || !r.eq_ignore_ascii_case(&original_ref)
        || !a.eq_ignore_ascii_case(&original_alt);
    (cur_pos, r, a, was_modified)
}
