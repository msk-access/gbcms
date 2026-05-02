#!/usr/bin/env python3
"""
Compare py-gbcms output against DMP ground truth (C++ GBCMS counts).

Each output MAF contains ALL input variants genotyped against ONE BAM.
We use --bam-map to select only results from the correct BAM per sample.

Usage:
    python compare_results.py \
        --gbcms-dir /path/to/gbcms_output \
        --truth /path/to/sample_mutations.tsv \
        --bam-map /path/to/cluster_bam_list.tsv \
        --output /path/to/concordance_report.tsv
"""

import argparse
import csv
import sys
from pathlib import Path

csv.field_size_limit(sys.maxsize)


def safe_int(v):
    try:
        return int(v) if v and str(v).strip() else 0
    except ValueError:
        return 0


def safe_float(v):
    try:
        return float(v) if v and str(v).strip() else 0.0
    except ValueError:
        return 0.0


def load_bam_map(bam_map_path: Path) -> dict[str, str]:
    """Load sample -> BAM stem mapping.

    Returns dict mapping sample ID to the BAM filename stem
    (e.g., 'P-0042513-T01-IM6' -> 'WI498341-T').
    """
    mapping = {}
    with open(bam_map_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                sample = parts[0]
                bam_path = parts[1]
                # Extract stem: /data1/.../WI498341-T.bam -> WI498341-T
                stem = Path(bam_path).stem
                mapping[sample] = stem
    return mapping


def load_gbcms_outputs(gbcms_dir: Path, sample_to_bam_stem: dict[str, str]) -> dict:
    """Load gbcms output MAFs, only reading each sample's results from its own BAM.

    Returns dict indexed by (sample, chrom, start, ref, alt).
    """
    # Invert: bam_stem -> set of samples that belong to it
    stem_to_samples: dict[str, set[str]] = {}
    for sample, stem in sample_to_bam_stem.items():
        stem_to_samples.setdefault(stem, set()).add(sample)

    results = {}
    loaded_files = 0
    for maf_file in sorted(gbcms_dir.glob("*.maf")):
        file_stem = maf_file.stem  # e.g., "WI498341-T"
        valid_samples = stem_to_samples.get(file_stem, set())
        if not valid_samples:
            continue  # No samples map to this BAM

        loaded_files += 1
        with open(maf_file) as f:
            # Skip comment lines
            while True:
                pos = f.tell()
                line = f.readline()
                if not line.startswith("#"):
                    f.seek(pos)
                    break
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                sample = row.get("Tumor_Sample_Barcode", "")
                if sample not in valid_samples:
                    continue  # Skip variants from other samples
                chrom = row.get("Chromosome", "")
                start = row.get("Start_Position", "")
                ref = row.get("Reference_Allele", "")
                alt = row.get("Tumor_Seq_Allele2", "")
                key = (sample, chrom, start, ref, alt)
                results[key] = row

    print(f"  Loaded {loaded_files} MAF files, {len(results)} sample-matched records")
    return results


def main():
    parser = argparse.ArgumentParser(description="Compare py-gbcms vs DMP ground truth")
    parser.add_argument("--gbcms-dir", required=True, help="Directory with py-gbcms output MAFs")
    parser.add_argument("--truth", required=True, help="DMP ground truth TSV (sample_mutations.tsv)")
    parser.add_argument("--bam-map", required=True,
                        help="Sample-to-BAM mapping TSV (cluster_bam_list.tsv)")
    parser.add_argument("--output", required=True, help="Output concordance TSV")
    args = parser.parse_args()

    gbcms_dir = Path(args.gbcms_dir)
    truth_path = Path(args.truth)
    bam_map_path = Path(args.bam_map)
    output_path = Path(args.output)

    # Load BAM mapping
    print(f"Loading BAM map from {bam_map_path}...")
    sample_to_bam = load_bam_map(bam_map_path)
    print(f"  {len(sample_to_bam)} sample -> BAM mappings")

    # Load gbcms outputs (only correct BAM per sample)
    print(f"Loading py-gbcms outputs from {gbcms_dir}...")
    gbcms = load_gbcms_outputs(gbcms_dir, sample_to_bam)

    # Load DMP truth
    print(f"Loading DMP ground truth from {truth_path}...")
    truth_rows = []
    with open(truth_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            truth_rows.append(row)
    print(f"  Loaded {len(truth_rows)} truth records")

    # Compare
    with open(output_path, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow([
            "sample", "gene", "chrom", "start", "ref", "alt", "variant_type",
            "dmp_ref", "dmp_alt", "dmp_depth",
            "gbcms_ref", "gbcms_alt", "gbcms_depth", "gbcms_vaf",
            "gbcms_rdf", "gbcms_adf", "gbcms_dpf",
            "alt_delta", "alt_pct_diff",
            "frag_ref_ok", "frag_alt_ok", "frag_sum_ok",
            "matched",
        ])

        stats = {"total": 0, "matched": 0, "exact_alt": 0, "within_5": 0, "within_10": 0,
                 "frag_ok": 0, "by_type": {}}

        for row in truth_rows:
            sample = row["Tumor_Sample_Barcode"]
            chrom = row["Chromosome"]
            start = row["Start_Position"]
            ref = row["Reference_Allele"]
            alt = row["Tumor_Seq_Allele2"]
            vtype = row.get("Variant_Type", "UNK")
            gene = row.get("Hugo_Symbol", "")
            dmp_ref = safe_int(row.get("t_ref_count"))
            dmp_alt = safe_int(row.get("t_alt_count"))
            dmp_depth = safe_int(row.get("t_depth"))

            key = (sample, chrom, start, ref, alt)
            gbcms_row = gbcms.get(key)

            stats["total"] += 1
            stats["by_type"].setdefault(vtype, {"total": 0, "matched": 0, "exact": 0, "w5": 0, "w10": 0})
            stats["by_type"][vtype]["total"] += 1

            if gbcms_row is None:
                writer.writerow([
                    sample, gene, chrom, start, ref, alt, vtype,
                    dmp_ref, dmp_alt, dmp_depth,
                    "", "", "", "", "", "", "",
                    "", "", "", "", "", "False",
                ])
                continue

            stats["matched"] += 1
            stats["by_type"][vtype]["matched"] += 1

            # Extract gbcms counts (column names depend on --column-prefix)
            g_ref = safe_int(gbcms_row.get("ref_count", gbcms_row.get("t_ref_count", 0)))
            g_alt = safe_int(gbcms_row.get("alt_count", gbcms_row.get("t_alt_count", 0)))
            g_dp = safe_int(gbcms_row.get("total_count", gbcms_row.get("t_depth", 0)))
            g_vaf = safe_float(gbcms_row.get("vaf", gbcms_row.get("t_alt_freq", 0)))
            g_rdf = safe_int(gbcms_row.get("ref_count_fragment", gbcms_row.get("t_ref_count_fragment", 0)))
            g_adf = safe_int(gbcms_row.get("alt_count_fragment", gbcms_row.get("t_alt_count_fragment", 0)))
            g_dpf = safe_int(gbcms_row.get("total_count_fragment", gbcms_row.get("t_depth_fragment", 0)))

            alt_delta = g_alt - dmp_alt
            alt_pct = (alt_delta / dmp_alt * 100) if dmp_alt > 0 else 0

            frag_ref_ok = g_rdf <= g_ref
            frag_alt_ok = g_adf <= g_alt
            frag_sum_ok = (g_rdf + g_adf) <= g_dpf

            if alt_delta == 0:
                stats["exact_alt"] += 1
                stats["by_type"][vtype]["exact"] += 1
            if abs(alt_delta) <= max(1, dmp_alt * 0.05):
                stats["within_5"] += 1
                stats["by_type"][vtype]["w5"] += 1
            if abs(alt_delta) <= max(1, dmp_alt * 0.10):
                stats["within_10"] += 1
                stats["by_type"][vtype]["w10"] += 1
            if all([frag_ref_ok, frag_alt_ok, frag_sum_ok]):
                stats["frag_ok"] += 1

            writer.writerow([
                sample, gene, chrom, start, ref, alt, vtype,
                dmp_ref, dmp_alt, dmp_depth,
                g_ref, g_alt, g_dp, f"{g_vaf:.4f}",
                g_rdf, g_adf, g_dpf,
                alt_delta, f"{alt_pct:.1f}",
                frag_ref_ok, frag_alt_ok, frag_sum_ok, "True",
            ])

    # Print summary
    print()
    print("=" * 80)
    print("CONCORDANCE SUMMARY: py-gbcms (MNP fix) vs DMP C++ GBCMS")
    print("=" * 80)
    m = stats["matched"]
    print(f"Total truth variants:  {stats['total']}")
    print(f"Matched in output:     {m}")
    print(f"Fragment invariants:   {stats['frag_ok']}/{m} ({100*stats['frag_ok']/m:.0f}%)" if m else "")
    print(f"Exact alt match:       {stats['exact_alt']}/{m} ({100*stats['exact_alt']/m:.0f}%)" if m else "")
    print(f"Within 5%:             {stats['within_5']}/{m} ({100*stats['within_5']/m:.0f}%)" if m else "")
    print(f"Within 10%:            {stats['within_10']}/{m} ({100*stats['within_10']/m:.0f}%)" if m else "")
    print()
    print(f"{'Type':>5s}  {'Total':>5s}  {'Match':>5s}  {'Exact':>5s}  {'≤5%':>5s}  {'≤10%':>5s}")
    print("-" * 40)
    for vtype in ["SNP", "DNP", "TNP", "ONP", "INS", "DEL"]:
        t = stats["by_type"].get(vtype, {})
        if t.get("total", 0) > 0:
            print(f"{vtype:>5s}  {t['total']:>5d}  {t['matched']:>5d}  "
                  f"{t['exact']:>5d}  {t['w5']:>5d}  {t['w10']:>5d}")
    print()
    print(f"Results written to: {output_path}")


if __name__ == "__main__":
    main()
