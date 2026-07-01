/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FILTER_MAF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Pre-filter a multi-sample MAF to extract variants for a specific sample.

    Matching modes:
    - Exact: meta.id must equal Tumor_Sample_Barcode (default)
    - Regex: meta.tsb pattern(s) searched against Tumor_Sample_Barcode
    - Multi-select: comma-separated patterns in meta.tsb, any match accepted

    Deduplication via composite key prevents duplicate rows when overlapping
    patterns match the same variant.

    Outputs:
    - *.filtered.maf  — filtered variants (header + matching data rows)
    - *.filter_stats.tsv — per-sample filtering statistics for pipeline summary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FILTER_MAF {
    tag "${meta.id}${meta.suffix ?: ''}"
    label 'process_low'

    input:
    tuple val(meta), path(maf)

    output:
    tuple val(meta), path("*.filtered.maf")    , emit: maf
    tuple val(meta), path("*.filter_stats.tsv"), emit: stats
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Use explicit TSB mapping if provided, otherwise fall back to meta.id
    def suffix = meta.suffix ?: ''
    def sample_id = meta.id
    def output_id = "${meta.id}${suffix}"
    def tsb_patterns = meta.tsb ?: ''
    """
    python3 << 'PYEOF'
import re
import sys

sample_id = '${sample_id}'
tsb_raw = '${tsb_patterns}'

# ── Build matching strategy ──────────────────────────────────────────
if tsb_raw:
    # Regex mode: comma-separated patterns, match ANY
    patterns = [re.compile(p.strip()) for p in tsb_raw.split(",") if p.strip()]

    def match(tsb):
        return any(p.search(tsb) for p in patterns)

    mode = f"regex ({len(patterns)} pattern(s): {tsb_raw})"
else:
    # Exact mode: sample column must equal Tumor_Sample_Barcode
    def match(tsb):
        return tsb == sample_id

    mode = f"exact ('{sample_id}')"

# ── Counters ─────────────────────────────────────────────────────────
matched = 0
skipped = 0
duplicates = 0
seen_tsbs = set()
seen_rows = set()  # Dedup key: (Chrom, Start, Ref, Alt, TSB)

# ── Filter ───────────────────────────────────────────────────────────
with open("${maf}") as infile, open(f"${output_id}.filtered.maf", "w") as outfile:
    header_written = False
    tsb_col = None
    chrom_col = start_col = ref_col = alt_col = None

    for line in infile:
        # Preserve comment lines
        if line.startswith("#"):
            outfile.write(line)
            continue

        fields = line.rstrip("\\n").split("\\t")

        # Parse header row
        if not header_written:
            outfile.write(line)
            header_written = True
            col_map = {col: i for i, col in enumerate(fields)}
            tsb_col = col_map.get("Tumor_Sample_Barcode")
            chrom_col = col_map.get("Chromosome")
            start_col = col_map.get("Start_Position")
            ref_col = col_map.get("Reference_Allele")
            alt_col = col_map.get("Tumor_Seq_Allele2")
            if tsb_col is None:
                print(
                    "ERROR: Tumor_Sample_Barcode column not found in MAF header",
                    file=sys.stderr,
                )
                sys.exit(1)
            continue

        # Data rows: check match and dedup
        tsb_value = fields[tsb_col] if len(fields) > tsb_col else ""
        seen_tsbs.add(tsb_value)

        if match(tsb_value):
            # Build composite key for deduplication
            dedup_key = tuple(
                fields[c] if c is not None and c < len(fields) else ""
                for c in (chrom_col, start_col, ref_col, alt_col, tsb_col)
            )
            if dedup_key in seen_rows:
                duplicates += 1
                continue
            seen_rows.add(dedup_key)
            outfile.write(line)
            matched += 1
        else:
            skipped += 1

# ── Summary ──────────────────────────────────────────────────────────
total = matched + skipped + duplicates
status = "processed" if matched > 0 else "skipped"
matched_tsb_list = sorted(t for t in seen_tsbs if match(t)) if matched > 0 else []

# Write stats TSV for pipeline summary collection
with open(f"${output_id}.filter_stats.tsv", "w") as sf:
    sf.write(
        "sample\\tstatus\\tvariants_matched\\tvariants_total"
        "\\tduplicates_removed\\tmatched_tsbs\\tmode\\n"
    )
    sf.write(
        f"{sample_id}\\t{status}\\t{matched}\\t{total}"
        f"\\t{duplicates}\\t{';'.join(matched_tsb_list)}\\t{mode}\\n"
    )

# Log to stderr for Nextflow task log
print(
    f"FILTER_MAF [{sample_id}]: {matched}/{total} variants matched ({mode})",
    file=sys.stderr,
)
if duplicates:
    print(f"  Removed {duplicates} duplicate row(s)", file=sys.stderr)
if matched == 0:
    print(
        f"WARNING: 0 matches. TSBs in MAF: {sorted(seen_tsbs)}", file=sys.stderr
    )
elif matched < total:
    print(f"  Matched TSBs: {matched_tsb_list}", file=sys.stderr)
PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def output_id = "${meta.id}${meta.suffix ?: ''}"
    """
    echo -e "Hugo_Symbol\\tChromosome\\tStart_Position\\tTumor_Sample_Barcode" > ${output_id}.filtered.maf
    echo -e "sample\\tstatus\\tvariants_matched\\tvariants_total\\tduplicates_removed\\tmatched_tsbs\\tmode" > ${output_id}.filter_stats.tsv
    echo -e "${output_id}\\tskipped\\t0\\t0\\t0\\t\\texact" >> ${output_id}.filter_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
