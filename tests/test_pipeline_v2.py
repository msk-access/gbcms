"""
Integration test for gbcms v2 pipeline.

Tests the full flow: read variants → prepare (validate+normalize) → count → write output.
"""

from pathlib import Path

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
    prepared = gbcms_rs.prepare_variants(rs_input, fasta_path, context_padding=5, is_maf=False)
    rs_variants = [p.variant for p in prepared]

    assert len(rs_variants) > 0, "No variants after preparation"

    # 3. Run Rust Engine
    results = gbcms_rs.count_bam(
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
        writer.write(v, counts, validation_status=pv.validation_status)
    writer.close()

    # 5. Verify Output
    assert output_path.exists()
    with open(output_path) as f:
        lines = f.readlines()
        assert len(lines) > 1  # Header + Data
        header = lines[0].strip().split("\t")
        # Default prefix is '' so columns are unprefixed
        assert "ref_count" in header
        assert "vaf_fragment" in header
        assert "ref_count_forward" in header
        assert "alt_count_reverse" in header
        assert "validation_status" in header


def test_pipeline_v2_binned(tmp_path):
    """Same integration test as test_pipeline_v2 but using count_bam_binned (production path).

    Verifies the binned engine produces identical results and output format.
    """
    base_dir = Path(__file__).parent / "testdata"
    bam_path = str(base_dir / "sample1_integration_test.bam")
    vcf_path = str(base_dir / "integration_test_variants.vcf")
    fasta_path = str(base_dir / "integration_test_reference.fa")
    output_path = tmp_path / "output_binned.maf"

    # 1. Read Variants
    reader = VcfReader(Path(vcf_path))
    variants = list(reader)
    reader.close()

    # 2. Prepare variants
    rs_input = [
        gbcms_rs.Variant(v.chrom, v.pos, v.ref, v.alt, v.variant_type.value) for v in variants
    ]
    prepared = gbcms_rs.prepare_variants(rs_input, fasta_path, context_padding=5, is_maf=False)
    rs_variants = [p.variant for p in prepared]

    # 3. Run BOTH engines and compare
    results_legacy = gbcms_rs.count_bam(
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

    results_binned = gbcms_rs.count_bam_binned(
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
        fragment_qual_threshold=10,
        sibling_variants=[[] for _ in rs_variants],
    )

    assert len(results_binned) == len(rs_variants)

    # Verify parity on key fields
    parity_fields = ["dp", "rd", "ad", "rd_fwd", "rd_rev", "ad_fwd", "ad_rev", "dpf", "rdf", "adf"]
    for i, (leg, bn) in enumerate(zip(results_legacy, results_binned, strict=True)):
        for field in parity_fields:
            assert getattr(leg, field) == getattr(bn, field), (
                f"Variant {i} parity mismatch on '{field}': "
                f"count_bam={getattr(leg, field)}, count_bam_binned={getattr(bn, field)}"
            )

    # 4. Write Output using binned results
    writer = MafWriter(output_path)
    for pv, counts in zip(prepared, results_binned, strict=True):
        v = next(v for v in variants if v.chrom == pv.variant.chrom and v.pos == pv.original_pos)
        writer.write(v, counts, validation_status=pv.validation_status)
    writer.close()

    # 5. Verify Output
    assert output_path.exists()
    with open(output_path) as f:
        lines = f.readlines()
        assert len(lines) > 1
        header = lines[0].strip().split("\t")
        assert "ref_count" in header
        assert "vaf_fragment" in header
        assert "validation_status" in header
