"""Reads of repeats where the guard's masks take the base telling the alleles apart, and
reads the aligner placed a unit or a copy off.

In a repeat whose units differ at single bases, the base past REF's tract is often the
only one telling a one-unit insertion from REF over a short read; a germline change there
makes the other haplotype read as the insertion (HG002's 28-bp VNTR row). The guard masks
a change the REF molecules carry, and then:

- the windows read on to the next base telling the alleles apart, which prep holds past
  every base the guard can mask; a read short of it shows neither allele;
- a change giving the reads the ALT's own base, shifted onto them where the aligner wrote no
  length change, is a germline change only where reads placed apart from the event carry
  it; otherwise it is the ALT's own base, never masked, and spells no other allele;
- an error both readings share decides nothing, however the read is written;
- a read the aligner placed a unit or a copy off is read whole over all its bases (the
  whole-read reach prep holds): the allele its bases fit decides, and one fitting the
  longer allele only from inside its event is depth only.

Every shape but the last is a reference and two haplotypes with paired reads simulated from
each and placed by BWA-MEM 0.7.17 as production runs it (tests/testdata/README.md); each
read's call is checked against the haplotype it came from (read from the per-read trace).
"""

import collections
import json
import logging
import random
import re
import zlib
from pathlib import Path

import pysam
import pytest
from helpers import count_checked, make_read, write_contig

from gbcms import _rs

FIXTURE = json.loads((Path(__file__).parent / "testdata" / "masked_repeat_reads.json").read_text())
ENGINE_LOGGER = "_rs.counting.engine"
CALL = re.compile(
    r"read call \S+ \S+ read=(h[12]|ref)_(\S+) mate=(\d) ref=(\w+) alt=(\w+) .*?partial=(\d+).*?uninformative=(\w+)"
)
# The #290 limit: a read lying inside an imperfect repeat, or read from a flank the repeat
# runs on past, is placed by the aligner, not by its bases, until the windows anchor on the
# repeat's true ends.
IMPERFECT_REPEAT = "#290: windows anchor on the exact repeat's edge, inside an imperfect one"


def _counted(
    tmp_path,
    caplog,
    name,
    thin=None,
    seed=0,
    error=None,
    change=None,
    one_orientation=False,
    lowq=None,
):
    """The shape's counts, and each haplotype's read calls (REF, ALT, partial, depth only,
    other), from one call. `thin` (haplotype, share kept, whether the rest read as the
    reference) keeps that haplotype on a share of its fragments (drawn with `seed`), the
    others read as the reference (a third haplotype of its length, tallied as `ref`) or
    dropped. `error` (haplotype, genome base, a later base, the base read, reads[, F1R2
    only]) gives that many of the haplotype's reads holding both bases with no length
    change between them (of one pair orientation, read 1 forward or read 2 reverse, if so
    asked) a confident error at the first. `lowq` (haplotype, genome base, reads) has that
    many of the haplotype's reads that read another base than the reference's there read
    it at quality 2. `change` (haplotype, base, the base it reads) gives that
    haplotype another change, its reads placed where they were. Half the fragments take
    the other pair orientation, as a library's do, unless `one_orientation` (an amplicon
    library's)."""
    f = dict(FIXTURE[name])
    if change:
        hap, at, base = change
        f[hap] = f[hap][:at] + base + f[hap][at + 1 :]
    swapped, as_ref = set(), False
    if thin:
        hap, keep, as_ref = thin
        rng = random.Random(seed)
        swapped = {
            q
            for q in dict.fromkeys(r[0] for r in f["reads"] if r[5] == hap)
            if rng.random() >= keep
        }
    fa = tmp_path / f"{name}.fa"
    fa.write_text(f">1\n{f['ref']}\n")
    pysam.faidx(str(fa))
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "1", "LN": len(f["ref"])}]}
    raw, bam = tmp_path / f"{name}.raw.bam", tmp_path / f"{name}.bam"
    erred = lowered = 0
    with pysam.AlignmentFile(str(raw), "wb", header=header) as out:
        for qname, flag, pos, cigar, mapq, hap, off, length, mate_pos, tlen in f["reads"]:
            if qname in swapped and not as_ref:
                continue
            if qname in swapped:
                qname, hap = "ref" + qname[2:], "ref"
            if not one_orientation and zlib.crc32(qname[qname.index("_") :].encode()) & 1:
                flag ^= 0xC0  # half the fragments the other pair orientation, as a library holds
            a = pysam.AlignedSegment()
            a.query_name, a.flag, a.reference_id, a.reference_start = qname, flag, 0, pos
            a.cigarstring, a.mapping_quality = cigar, mapq
            a.query_sequence = seq = f[hap][off : off + length]
            f1r2 = bool(a.flag & 0x10) == bool(a.flag & 0x80)
            if error and hap == error[0] and erred < error[4] and (len(error) < 6 or f1r2):
                at, tell, base = error[1], error[2], error[3]
                q = {r: q for q, r in a.get_aligned_pairs(matches_only=True)}
                if at in q and tell in q and q[tell] - q[at] == tell - at and seq[q[at]] != base:
                    a.query_sequence = seq[: q[at]] + base + seq[q[at] + 1 :]
                    erred += 1
            quals = "I" * length
            if lowq and hap == lowq[0] and lowered < lowq[2]:
                q = {r: q for q, r in a.get_aligned_pairs(matches_only=True)}
                at = lowq[1]
                if at in q and a.query_sequence[q[at]] != f["ref"][at]:
                    quals = quals[: q[at]] + "#" + quals[q[at] + 1 :]
                    lowered += 1
            a.query_qualities = pysam.qualitystring_to_array(quals)
            a.next_reference_id, a.next_reference_start, a.template_length = 0, mate_pos, tlen
            out.write(a)
    pysam.sort("-o", str(bam), str(raw))
    pysam.index(str(bam))
    p = f["pos"]  # the base after the anchor (0-based): the insertion's junction
    anchor = f["ref"][p - 1]
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", p - 1, anchor, anchor + f["alt"], "INS")], str(fa), 5, False, 1, True
    )
    caplog.clear()
    with caplog.at_level(5, logger=ENGINE_LOGGER):
        _rs.reset_log_caching()
        (counts,) = count_checked(str(bam), [pv.variant])
    calls = {}
    for record in caplog.records:
        if m := CALL.search(record.getMessage()):
            hap, read, mate, ref, alt, partial, depth = m.groups()
            shown = {
                "REF": ref == "true",
                "ALT": alt == "true",
                "partial": partial != "0",
                "depth": depth == "true",
            }
            calls[(hap, read, mate)] = next((c for c, on in shown.items() if on), "other")
    tally = collections.Counter((hap, call) for (hap, _, _), call in calls.items())
    return counts, tally


# ── The VNTR of HG002's 28-bp row: C and T units differing at one base ──────────────────
def test_reads_of_two_haplotypes_reading_alike_are_neither(tmp_path, caplog):
    """The insertion of a C unit on one haplotype, and on the other a germline SNP turning
    the first T unit into a C unit: over 151 bp both read four C units and a T unit, first
    telling apart 160 bases on. The guard masks the SNP (reads placed apart carry it), the
    windows read on to that base, and no read reaches it: every read is depth only. Neither
    the SNP's reads as ALT (un-applying the mask did), nor the insertion's as REF (the CIGAR
    checks did, where the rule gave up)."""
    counts, calls = _counted(tmp_path, caplog, "vntr28b_ins_and_snp")
    assert (counts.rd, counts.ad) == (0, 0)
    assert not any(c in ("REF", "ALT") for (_, c) in calls.elements()), dict(calls)


def test_the_insertion_alone_counts_where_the_aligner_wrote_no_gap(tmp_path, caplog):
    """The insertion alone, its reads written without the gap: the change it shifts onto
    them past REF's tract is no germline change (reads placed apart from the event never
    carry it), so it is not masked and the insertion's reads count ALT by their bases."""
    counts, calls = _counted(tmp_path, caplog, "vntr28b_ins_alone")
    assert (counts.rd, counts.ad) == (27, 34)
    assert calls[("h1", "REF")] == 0 and calls[("h2", "ALT")] == 0, dict(calls)
    assert counts.possible_minority == "", "no read placed apart carries it"


def test_a_germline_snp_reading_as_the_insertion_never_counts_alt(tmp_path, caplog):
    """The SNP alone, no insertion: its haplotype reads as the insertion over a short read.
    Reads placed apart carry the change, so it is masked; the SNP's reads show neither
    allele and never count ALT (as un-applying the mask made them, at VAF 0.5)."""
    counts, calls = _counted(tmp_path, caplog, "vntr28b_snp_alone")
    assert (counts.rd, counts.ad) == (34, 0)
    assert calls[("h1", "ALT")] == 0, dict(calls)
    assert counts.possible_minority == "", "masked: a haplotype's"


# A haplotype at a share of the molecules (a tumour's gain or loss, a subclone): reads placed
# apart from the event carry its change beyond an error, however small its share of them.
@pytest.mark.parametrize("seed", range(3))
def test_a_minority_snp_haplotype_in_trans_never_counts_alt(tmp_path, caplog, seed):
    """The insertion on half the molecules, the SNP on a quarter of the others, the
    reference on the rest: reads placed apart carry the SNP on under a fifth of them (the
    insertion's own read the reference's base there), yet a haplotype carries it. Masked,
    the SNP's reads never count ALT."""
    _, calls = _counted(tmp_path, caplog, "vntr28b_ins_and_snp", ("h2", 0.25, True), seed)
    assert calls[("h2", "ALT")] == 0, dict(calls)


@pytest.mark.parametrize("seed", range(3))
def test_the_snp_haplotype_beside_a_gained_insertion_never_counts_alt(tmp_path, caplog, seed):
    """No reference haplotype: the insertion's gained, the SNP's on a quarter as many
    molecules."""
    _, calls = _counted(tmp_path, caplog, "vntr28b_ins_and_snp", ("h2", 0.25, False), seed)
    assert calls[("h2", "ALT")] == 0, dict(calls)


@pytest.mark.parametrize("seed", range(3))
def test_a_minority_germline_snp_reading_as_the_insertion_never_counts_alt(tmp_path, caplog, seed):
    """The SNP alone, on 15% of the molecules: too few of the REF fragments carry it for
    background, but reads placed apart carry it, so it is masked and its reads, reading
    as the insertion over a short read, never count ALT (develop counted them)."""
    _, calls = _counted(tmp_path, caplog, "vntr28b_snp_alone", ("h1", 0.15, True), seed)
    assert calls[("h1", "ALT")] == 0, dict(calls)


def test_a_germline_haplotype_reading_as_the_alt_at_the_next_telling_base_never_counts_alt(
    tmp_path, caplog
):
    """The SNP alone, its haplotype also carrying the shifted ALT's base at the next base
    telling them apart (G900T): its reads never show the reference there, but they do at
    the bases past it, so they are placed apart and carry the change; masked, the
    haplotype's reads never count ALT (with no insertion in the sample)."""
    counts, calls = _counted(tmp_path, caplog, "vntr28b_snp_alone", change=("h1", 900, "T"))
    assert counts.ad == 0 and calls[("h1", "ALT")] == 0, dict(calls)


@pytest.mark.parametrize("errors", [1, 2])
def test_an_error_at_the_shifted_base_leaves_the_insertion_counted(tmp_path, caplog, errors):
    """The insertion alone, its reads error-free but for one or two REF reads placed apart
    from the event with a confident error giving the ALT's shifted base (base 829, the C
    an insertion read written without its gap shows past REF's tract; 900 the next base
    telling them apart). Errors explain so few: still the ALT's own base, not masked, and
    the insertion's reads count ALT as without the error (masked, every one was depth
    only)."""
    counts, calls = _counted(
        tmp_path, caplog, "vntr28b_ins_alone", error=("h2", 829, 900, "C", errors)
    )
    assert (counts.rd, counts.ad) == (27, 34)
    # unmasked though carried: named from two carriers (POSSIBLE_MINORITY_HAPLOTYPE; one
    # is an error's), the counts unchanged
    named = counts.possible_minority.startswith(f"830:T>C:{errors}/")
    assert named == (errors >= 2), counts.possible_minority


def test_a_germline_snp_in_a_one_orientation_library_never_counts_alt(tmp_path, caplog):
    """The SNP alone, every fragment in one pair orientation (an amplicon library): its
    carriers are all one-sided, as every read is, so orientation tells nothing and the
    SNP is still a haplotype's change, masked; its reads never count ALT."""
    counts, calls = _counted(tmp_path, caplog, "vntr28b_snp_alone", one_orientation=True)
    assert counts.ad == 0 and calls[("h1", "ALT")] == 0, dict(calls)


@pytest.mark.parametrize("artifacts", [12, 20])
def test_an_artifact_in_one_pair_orientation_leaves_the_insertion_counted(
    tmp_path, caplog, artifacts
):
    """The insertion alone, with 12 or 20 REF reads placed apart giving the shifted ALT
    base, all of one pair orientation (an oxoG or FFPE artifact sits in one): over 3% of
    the reads placed apart, but a haplotype's molecules are read in both orientations, and
    so many carriers falling one way is too unlikely for a haplotype's (under 0.1% against
    the library's balance): no haplotype's change, it stays the ALT's own base, and the
    insertion's reads count ALT as without it. Eight carriers all one way are a
    haplotype's chance (0.4%), so they are masked."""
    error = ("h2", 829, 900, "C", artifacts, "F1R2")
    counts, _ = _counted(tmp_path, caplog, "vntr28b_ins_alone", error=error)
    assert (counts.rd, counts.ad) == (27, 34)
    assert counts.possible_minority.startswith(f"830:T>C:{artifacts}/"), counts.possible_minority


@pytest.mark.parametrize(
    "name",
    [
        pytest.param(n, marks=pytest.mark.xfail(strict=True, reason=IMPERFECT_REPEAT))
        for n in ("vntr28_ins_and_snp", "vntr28_ins_and_snp_20x")
    ]
    + ["vntr28_ins_and_snp_250", "vntr28_ins_alone", "vntr28_ins_alone_20x"],
)
def test_no_insertion_read_counts_ref_in_an_imperfect_repeat(tmp_path, caplog, name):
    """Insertion reads lying inside the repeat, read from a flank the imperfect repeat runs
    on past and placed there by the aligner a unit off: their bases fit REF as written.
    Decided past the windows, the ALT placed a unit over fits them as well: depth only.
    Those decided within the windows still count REF (#290)."""
    _, calls = _counted(tmp_path, caplog, name)
    assert calls[("h1", "REF")] == 0, dict(calls)


@pytest.mark.parametrize("name", ["vntr28_snp_alone", "vntr28_snp_alone_20x"])
@pytest.mark.xfail(strict=True, reason=IMPERFECT_REPEAT)
def test_no_snp_read_counts_alt_in_an_imperfect_repeat(tmp_path, caplog, name):
    """SNP reads lying wholly inside the repeat, the aligner's tie placing them a unit off:
    read from that flank they look like one more unit."""
    _, calls = _counted(tmp_path, caplog, name)
    assert calls[("h1", "ALT")] == 0, dict(calls)


def test_the_repeat_without_either_change_counts_no_alt(tmp_path, caplog):
    counts, calls = _counted(tmp_path, caplog, "vntr28_ref_only")
    assert counts.ad == 0 and calls[("h1", "ALT")] + calls[("h2", "ALT")] == 0


# ── Reads the aligner placed a unit or a copy off ──────────────────────────────────────
def test_a_read_placed_a_unit_off_is_judged_whole(tmp_path, caplog):
    """A VNTR of 40-bp units differing at single bases, a compound het: the given unit on
    one haplotype, another unit upstream on the other (the same length). BWA writes the
    insertion's reads ungapped a unit off, where they fit REF over the windows: read whole
    over all their bases, they are the ALT."""
    counts, calls = _counted(tmp_path, caplog, "vntr40_compound")
    assert (counts.rd, counts.ad) == (66, 33)
    assert calls[("h1", "ALT")] == 33 and calls[("h2", "ALT")] == 0, dict(calls)


@pytest.mark.xfail(strict=True, reason=IMPERFECT_REPEAT)
def test_no_insertion_read_placed_a_unit_off_counts_ref(tmp_path, caplog):
    """Those whose bases, read whole, reach no unique sequence on either side."""
    _, calls = _counted(tmp_path, caplog, "vntr40_compound")
    assert calls[("h1", "REF")] == 0, dict(calls)


@pytest.mark.parametrize(
    ("name", "rd", "ad"), [("grch38_dup100", 18, 17), ("grch38_dup100_20x", 7, 5)]
)
def test_a_tandem_duplication_written_without_its_gap_counts_alt(tmp_path, caplog, name, rd, ad):
    """A 100-bp tandem duplication (a GRCh38 segment): the copy's bases the aligner shifts
    onto reads written without the gap are the ALT's own, never masked."""
    counts, calls = _counted(tmp_path, caplog, name)
    assert (counts.rd, counts.ad) == (rd, ad)
    assert calls[("h1", "REF")] == 0 and calls[("h2", "ALT")] == 0, dict(calls)


@pytest.mark.parametrize(
    ("name", "ref_reads"), [("grch38_dup89_cis", 32), ("grch38_dup89_cis_20x", 16)]
)
def test_no_insertion_read_counts_ref_where_it_fits_the_alt_a_copy_over(
    tmp_path, caplog, name, ref_reads
):
    """HG002's chr21 89-bp duplication (a GRCh38 segment): the inserted copy differs from
    the next 89 bases only at its last, where a germline SNP in cis gives that copy the
    same base. The guard masks the SNP (reads placed apart carry it) and the windows read
    on past it where the aligner wrote each read. The insertion's reads written through
    the copy then miss REF only at the anchor, the ALT's own last base shifted onto it;
    read there at low quality, they fit REF where they were written and the ALT exactly
    a copy over. A call decided past the windows stands only where the other allele,
    the whole read placed anywhere on it, misses more: those reads are depth only (HG002
    at this row: 2 counted REF). The REF haplotype's reads keep their calls."""
    f = FIXTURE[name]
    _, calls = _counted(tmp_path, caplog, name, lowq=("h1", f["pos"] - 1, len(f["reads"])))
    assert calls[("h1", "REF")] == 0 and calls[("h2", "ALT")] == 0, dict(calls)
    assert calls[("h2", "REF")] == ref_reads, dict(calls)


def test_shifted_alt_bases_spell_no_other_allele_of_refs_length(tmp_path, caplog):
    """A 23-bp tandem duplication in a CT-rich repeat (a GRCh38 segment): the ALT's bases
    shifted onto reads written without the gap are not masked, and they spell no other
    allele of REF's length (reads placed through the event carrying them are the ALT's):
    the REF haplotype's reads stay REF."""
    counts, calls = _counted(tmp_path, caplog, "grch38_dup23_20x")
    assert (counts.rd, counts.ad) == (25, 29)
    assert calls[("h2", "REF")] == 24 and calls[("h2", "other")] == 0, dict(calls)


def test_a_plain_insertion_is_unchanged(tmp_path, caplog):
    counts, calls = _counted(tmp_path, caplog, "plain_insertion")
    assert (counts.rd, counts.ad) == (83, 70)
    assert calls[("h1", "REF")] + calls[("h2", "ALT")] == 0


# ── An error both readings share ───────────────────────────────────────────────────────
def test_an_error_both_readings_share_decides_nothing(tmp_path):
    """A REF read ending on the base telling a 1-bp insertion from REF, that base unreadable
    (Q5) and a confident error on the anchor: both readings miss the anchor alone. Written
    without an indel, REF used to tolerate the error and count it; its bases tell nothing,
    so it is neither allele."""
    rng = __import__("random").Random(5)
    contig = list("".join(rng.choice("ACGT") for _ in range(400)))
    a = 200  # the insertion's junction; the anchor at a - 1
    contig[a - 1], contig[a] = "G", "C"
    contig = "".join(contig)
    reads = []
    for i in range(5):
        s = a - 90 + i
        seq, quals = list(contig[s : a + 1]), [30] * (a + 1 - s)
        if i == 0:
            seq[a - 1 - s] = "T"
            quals[a - s] = 5
        reads.append(make_read(f"r{i}", "".join(seq), s, ((0, len(seq)),), quals=quals))
    fa, bam = write_contig(tmp_path, contig, reads, "c")
    (pv,) = _rs.prepare_variants([_rs.Variant("1", a - 1, "G", "GT", "INS")], fa, 5, False, 1, True)
    (counts,) = count_checked(bam, [pv.variant])
    assert (counts.rd, counts.ad) == (4, 0)


@pytest.mark.parametrize("second", [5, 20, 40])
def test_an_error_both_readings_share_spends_no_tolerance_past_the_windows(tmp_path, second):
    """A REF read of a 1-bp deletion in an A10 run, its last A misread as T (both readings
    miss it alone) and one more error `second` bases past the run: read on past the
    windows with the shared error unread, its bases there fit REF but for that one error,
    so it counts REF as the twenty clean REF reads do. At 5 and 20 bases the error gives
    the ALT's base there (the ALT shifted by its length): past the windows one error the
    other allele's reading does not share is forgiven where that reading misses three or
    more (2026-10-21)."""
    rng = __import__("random").Random(3)
    left = "".join(rng.choice("ACGT") for _ in range(300))[:-1] + "G"
    right = "C" + "".join(rng.choice("ACGT") for _ in range(299))
    contig = left + "A" * 10 + right
    p = len(left)  # the run's first A
    alt = contig[:p] + contig[p + 1 :]
    reads = []
    for i in range(20):
        s = p - 100 + i
        reads.append(make_read(f"ref{i}", contig[s : s + 151], s, ((0, 151),)))
        k = p - s
        reads.append(make_read(f"alt{i}", alt[s : s + 150], s, ((0, k), (2, 1), (0, 150 - k))))
    s = p - 60
    seq = list(contig[s : s + 151])
    seq[p + 9 - s] = "T"
    q = p + 10 + second - s
    seq[q] = "A" if seq[q] != "A" else "C"
    reads.append(make_read("probe", "".join(seq), s, ((0, 151),)))
    fa, bam = write_contig(tmp_path, contig, reads, "c")
    (pv,) = _rs.prepare_variants([_rs.Variant("1", p - 1, "GA", "G", "DEL")], fa, 5, False, 1, True)
    (counts,) = count_checked(bam, [pv.variant])
    assert (counts.rd, counts.ad) == (21, 20)
    assert counts.n_count == 0, "the error read as unread is no N of the read's"


def test_one_read_of_another_insertion_of_the_alts_length_is_not_alt(tmp_path):
    """+CA in (CA)6, beside 12 REF and 12 ALT reads, one read of another 2-bp insertion
    (its seventh unit reads CC), every base Q30, written with its gap: its one mismatch from
    the ALT lies inside the event, where only an error both readings share is tolerated
    (2026-10-15). It is another allele (partial evidence), never the ALT, however few
    fragments carry it: one error forgiven is one past the windows (2026-10-21)."""
    rng = __import__("random").Random(11)
    left = "".join(rng.choice("ACGT") for _ in range(300))[:-1] + "G"
    right = "T" + "".join(rng.choice("ACGT") for _ in range(300))[1:]
    contig = left + "CA" * 6 + right
    p = len(left)
    alt = contig[:p] + "CA" + contig[p:]
    other = contig[: p + 12] + "CC" + contig[p + 12 :]
    reads = []
    for i in range(12):
        s = p - 100 + 3 * i
        reads.append(make_read(f"ref{i}", contig[s : s + 151], s, ((0, 151),)))
        reads.append(
            make_read(f"alt{i}", alt[s : s + 151], s, ((0, p - s), (1, 2), (0, 149 - (p - s))))
        )
    s = p - 60
    reads.append(make_read("other", other[s : s + 151], s, ((0, 72), (1, 2), (0, 77))))
    fa, bam = write_contig(tmp_path, contig, reads, "c")
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", p - 1, "G", "GCA", "INS")], fa, 5, False, 1, True
    )
    (counts,) = count_checked(bam, [pv.variant])
    assert (counts.rd, counts.ad, counts.partial_alt) == (12, 12, 1)


def test_a_ref_read_with_its_own_deletion_past_the_run_is_not_alt(tmp_path):
    """A 1-bp deletion in an A10 run, beside 20 REF and 20 ALT reads; a REF read whose
    last A is misread (both readings miss it, so it is read on) and which carries its own
    1-bp deletion of the second base past the run, written without its gap. Past the
    windows its bases are shifted one base, as the ALT's are, and miss the ALT at one base
    alone; but the read writes no length change, and an allele forgives an error only
    where the read writes its length (2026-10-15): it is no ALT."""
    rng = __import__("random").Random(3)
    left = "".join(rng.choice("ACGT") for _ in range(300))[:-1] + "G"
    right = "C" + "".join(rng.choice("ACGT") for _ in range(299))
    contig = left + "A" * 10 + right
    p = len(left)  # the run's first A
    assert contig[p + 10] != contig[p + 11], "the ALT misses the read at one base past the run"
    alt = contig[:p] + contig[p + 1 :]
    reads = []
    for i in range(20):
        s = p - 100 + i
        reads.append(make_read(f"ref{i}", contig[s : s + 151], s, ((0, 151),)))
        k = p - s
        reads.append(make_read(f"alt{i}", alt[s : s + 150], s, ((0, k), (2, 1), (0, 150 - k))))
    s = p - 60
    own = contig[: p + 11] + contig[p + 12 :]  # its own deletion of the second base past the run
    seq = list(own[s : s + 151])
    seq[p + 9 - s] = "T"
    reads.append(make_read("own_del", "".join(seq), s, ((0, 151),)))
    fa, bam = write_contig(tmp_path, contig, reads, "c")
    (pv,) = _rs.prepare_variants([_rs.Variant("1", p - 1, "GA", "G", "DEL")], fa, 5, False, 1, True)
    (counts,) = count_checked(bam, [pv.variant])
    assert counts.ad == 20


logging.getLogger(ENGINE_LOGGER)  # the logger the per-read trace comes through
