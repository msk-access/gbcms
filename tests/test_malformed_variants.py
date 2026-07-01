"""CR-3: malformed variants with an empty REF/ALT must not panic.

VCF/MAF indels are left-anchored — REF and ALT share a leading anchor base, so
both are non-empty for a well-formed record. A record with an empty ALT
(insertion) or empty REF (deletion) would underflow ``len() - 1`` and panic the
``[1..]``/``[0]`` slice inside ``check_insertion``/``check_deletion``, surfacing
as an opaque ``PyErr`` that kills the whole sample. The entry guards classify
such records as neither instead.

``count_one_both`` exercises BOTH the legacy and binned APIs (and asserts their
parity), so each test covers both code paths.
"""

from helpers import build_bam, count_one_both, make_read

from gbcms._rs import Variant


def _bam_with_one_read(tmp_path):
    """A single 40M read at chr1:100-140 so the variant locus is covered."""
    reads = [make_read("r1", "A" * 40, 100, ((0, 40),))]
    return build_bam(tmp_path, reads)


def test_empty_alt_insertion_does_not_panic(tmp_path):
    bam = _bam_with_one_read(tmp_path)
    # Empty ALT with INSERTION type — malformed (no inserted bases after anchor).
    variant = Variant("chr1", 110, "A", "", "INSERTION")
    counts = count_one_both(bam, variant)  # both APIs; asserts parity; must not raise
    assert counts.ad == 0, "malformed insertion must contribute no ALT"
    assert counts.adf == 0


def test_empty_ref_deletion_does_not_panic(tmp_path):
    bam = _bam_with_one_read(tmp_path)
    # Empty REF with DELETION type — malformed (no deleted bases after anchor).
    variant = Variant("chr1", 110, "", "A", "DELETION")
    counts = count_one_both(bam, variant)
    assert counts.ad == 0, "malformed deletion must contribute no ALT"
    assert counts.adf == 0
