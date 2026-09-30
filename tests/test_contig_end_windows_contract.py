"""Reference windows are clamped to the contig (#142).

The FASTA reader errors on a window that passes the contig end, and prep's
windows are padded on both sides but clamped only at the contig start. So an
indel within a few hundred bases of a contig end lost every window:

- the left-align window (about 100bp): the indel was not left-aligned (a WARN);
- `ref_context`: no haplotypes, so reads needing Phase 3 went to the SW fallback;
- the shift region (256bp) and the event reference (60bp): a read aligning the
  indel elsewhere in its repeat could not be matched to it, and a complex
  variant fell back to the tolerant classifier the exact-carrier rule replaced,
  so the same reads were judged differently near a contig end than elsewhere.

Each window now keeps the part of the contig it covers. Exact fetches (REF
validation, a MAF anchor) stay exact.

Committed red (xfail-strict) before the implementation; flipped green with it.
"""

import random

import pysam
import pytest
from helpers import count_both, make_read

from gbcms import _rs as gbcms_rs

L = 700  # contig length
READ = 100
RUN = "AAAAA"


def _contig(anchor):
    """A random contig with G + AAAAA + T starting at `anchor` (0-based)."""
    rng = random.Random(142)
    seq = [rng.choice("CGT") for _ in range(L)]
    seq[anchor : anchor + 7] = list("G" + RUN + "T")
    return "".join(seq)


# (label, 0-based anchor of the left-aligned deletion G+A > G)
SITES = {
    "10bp from the contig end": L - 14,
    "50bp from the contig end": L - 54,
    "40bp from the contig start (guard)": 40,
    "mid-contig (guard)": 340,
}


def _files(tmp_path, ref, reads):
    fa = tmp_path / "ref.fa"
    fa.write_text(f">1\n{ref}\n")
    pysam.faidx(str(fa))
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "1", "LN": len(ref)}]}
    raw = tmp_path / "raw.bam"
    with pysam.AlignmentFile(str(raw), "wb", header=hdr) as fh:
        for r in reads:
            fh.write(r)
    bam = tmp_path / "s.bam"
    pysam.sort("-o", str(bam), str(raw))
    pysam.index(str(bam))
    return str(fa), str(bam)


def _right_shifted(ref, anchor):
    """The deletion of one A as a caller may report it: anchored on the run's
    second-to-last A, deleting its last."""
    p = anchor + len(RUN) - 1
    return gbcms_rs.Variant("1", p, ref[p : p + 2], ref[p], "DELETION")


def _reads(ref, anchor):
    """Five REF reads and five carriers, the deletion aligned left (after the
    G), all inside the contig and spanning the run with flank on both sides, so
    every read is informative (one ending inside the run is neither)."""
    last = min(L - READ - 1, anchor - 10)
    starts = range(last - 4, last + 1)
    hap = ref[: anchor + 1] + ref[anchor + 2 :]
    reads = [make_read(f"r{i}", ref[s : s + READ], s, ((0, READ),)) for i, s in enumerate(starts)]
    for i, s in enumerate(starts):
        left = anchor + 1 - s
        reads.append(
            make_read(f"a{i}", hap[s : s + READ], s, ((0, left), (2, 1), (0, READ - left)))
        )
    return reads


@pytest.mark.parametrize("site", list(SITES))
def test_prep_left_aligns_and_fetches_every_window(tmp_path, site):
    anchor = SITES[site]
    ref = _contig(anchor)
    fa, _ = _files(tmp_path, ref, [])
    (pv,) = gbcms_rs.prepare_variants([_right_shifted(ref, anchor)], fa, 5, False, 1, True)
    assert pv.gbcms_status == "PASS", pv.gbcms_status_reason
    v = pv.variant
    assert (v.pos, v.ref_allele, v.alt_allele) == (anchor, "GA", "G"), "left-aligned to the G"
    assert pv.was_left_aligned
    assert v.ref_context is not None and v.ref_context_start <= anchor
    assert v.shift_region == (anchor + 1, anchor + 1 + len(RUN)), "the whole A run"
    assert v.event_ref is not None


@pytest.mark.parametrize("site", list(SITES))
def test_carriers_are_counted_by_the_haplotype_matrix(tmp_path, site):
    """Guard: the carriers align the deletion left of where the input puts it;
    they count as ALT, the REF reads as REF, and no read goes to the SW
    fallback (informative reads were counted right before too)."""
    anchor = SITES[site]
    ref = _contig(anchor)
    fa, bam = _files(tmp_path, ref, _reads(ref, anchor))
    (pv,) = gbcms_rs.prepare_variants([_right_shifted(ref, anchor)], fa, 5, False, 1, True)
    (c,) = count_both(bam, [pv.variant])
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev
    assert (c.rd, c.ad) == (5, 5)
    assert c.sw_fallback_reads == 0


def test_ref_validation_stays_exact_at_the_contig_end(tmp_path):
    """Guard: a REF that runs past the contig end fails validation; it is not
    matched against the clamped part (which could 'correct' it to a shorter
    REF)."""
    ref = _contig(340)
    fa, _ = _files(tmp_path, ref, [])
    tail = ref[L - 20 :] + "ACGT"  # 20 real bases, then 4 past the end
    (pv,) = gbcms_rs.prepare_variants(
        [gbcms_rs.Variant("1", L - 20, tail, tail[0], "DELETION")], fa, 5, False, 1, True
    )
    assert pv.gbcms_status == "FAIL", (pv.gbcms_status, pv.gbcms_status_reason)


def test_a_complex_variant_is_judged_near_the_contig_end_as_mid_contig(tmp_path):
    """The same 40 bases (a delins GG>TCC in a GGTT repeat) mid-contig and at the
    contig end, with the same reads relative to it: five REF reads, five exact
    carriers and three reads carrying TCG (another allele). The reads end inside
    the event grown through the repeat, so the exact-carrier rule finds none
    informative. At the contig end the rule had no reference and the tolerant
    classifier ran instead: REF 5, ALT 7 (two TCG reads among them)."""
    base = _contig(L - 32)
    tail = base[L - 40 :]  # the delins sits at offset 28
    mid = 300
    ref = base[:mid] + tail + base[mid + 40 :]
    counts = {}
    for label, p in (("mid-contig", mid + 28), ("contig end", L - 12)):
        reads = []
        for i, k in enumerate(range(90, 95)):  # each read ends 5-9 bases past p
            s = p - k
            reads.append(make_read(f"r{i}", ref[s : s + READ], s, ((0, READ),)))
            for tag, alt in (("a", "TCC"), ("n", "TCG"))[: 2 if i < 3 else 1]:
                hap = ref[:p] + alt + ref[p + 2 :]
                cig = ((0, k), (1, 1), (0, READ - k - 1))
                reads.append(make_read(f"{tag}{i}", hap[s : s + READ], s, cig))
        d = tmp_path / label.replace(" ", "_")
        d.mkdir()
        fa, bam = _files(d, ref, reads)
        (pv,) = gbcms_rs.prepare_variants(
            [gbcms_rs.Variant("1", p, ref[p : p + 2], "TCC", "COMPLEX")], fa, 5, False, 1, True
        )
        (c,) = count_both(bam, [pv.variant])
        counts[label] = (c.rd, c.ad, c.partial_alt, c.dp)
    assert counts["contig end"] == counts["mid-contig"], counts
