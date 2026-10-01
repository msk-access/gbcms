"""
Integration test for gbcms v2 pipeline.

Tests the full flow: read variants → prepare (validate+normalize) → count → write output.
Counting is checked for binning invariance (count_bam_checked compares every field
with one variant per bin).
"""

from pathlib import Path

from helpers import count_bam_checked

from gbcms import _rs as gbcms_rs
from gbcms.io.input import VcfReader
from gbcms.io.output import MafWriter


def test_pipeline_v2(tmp_path):
    base_dir = Path(__file__).parent / "testdata"
    bam_path = str(base_dir / "sample1_integration_test.bam")
    vcf_path = str(base_dir / "integration_test_variants.vcf")
    fasta_path = str(base_dir / "integration_test_reference.fa")
    output_path = tmp_path / "output.maf"

    # 1. Read Variants
    reader = VcfReader(Path(vcf_path))
    variants = list(reader)
    reader.close()

    assert len(variants) > 0

    # 2. Prepare variants (validate + normalize + ref_context)
    #    Note: The test reference FASTA is only 20kb, but the VCF variants are at
    #    chr1:11M+, so validation returns FETCH_FAILED. This is expected — this test
    #    verifies counting correctness, not validation. We use all variants.
    rs_input = [
        gbcms_rs.Variant(v.chrom, v.pos, v.ref, v.alt, v.variant_type.value) for v in variants
    ]
    prepared = gbcms_rs.prepare_variants(
        rs_input, fasta_path, context_padding=5, is_maf=False, threads=1, adaptive_context=False
    )
    rs_variants = [p.variant for p in prepared]

    assert len(rs_variants) > 0, "No variants after preparation"

    # 3. Count (production bins, checked against one variant per bin)
    results = count_bam_checked(
        bam_path,
        rs_variants,
        [None] * len(rs_variants),
        min_mapq=20,
        min_baseq=10,
        filter_duplicates=True,
        filter_secondary=False,
        filter_supplementary=False,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
        threads=1,
    )

    assert len(results) == len(rs_variants)

    # 4. Write Output
    writer = MafWriter(output_path)
    for pv, counts in zip(prepared, results, strict=True):
        # Use the original variant for output coords
        v = next(v for v in variants if v.chrom == pv.variant.chrom and v.pos == pv.original_pos)
        writer.write(v, counts, gbcms_status=pv.gbcms_status)
    writer.close()

    # 5. Verify Output
    assert output_path.exists()
    with open(output_path) as f:
        lines = [line for line in f.readlines() if not line.startswith("#")]
        assert len(lines) > 1  # Header + Data
        header = lines[0].strip().split("\t")
        # Default prefix is '' so columns are unprefixed
        assert "ref_count" in header
        assert "vaf_fragment" in header
        assert "ref_count_forward" in header
        assert "alt_count_reverse" in header
        assert "gbcms_status" in header
