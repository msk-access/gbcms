"""`TIED_WITH` and listed alleles named at the site (RJ-34): name what held-back reads fit.

Under the sequence judge some reads are held back: they count as depth only, or as another
allele (`partial_alt`), because their bases fit two candidate sequences equally well. This
diagnostic names, in `gbcms_diagnostic`, what those reads fit, and changes no count.

Shapes (made-up bases, allele shapes only):
- An insertion of one base and a listed SNV at the next base whose ALT base is that same
  inserted base. A read that ends right after the base fits both: the judge cannot tell the
  insertion from the SNV. Each row names the other: `TIED_WITH(chrom:pos:REF>ALT:n/m)`, with
  n the reads that fit both and m the row's own alt_count. Reads that run on past the base
  tell the two apart and are each row's own carriers.
- The same pair with the insertion not listed. The tie is with an allele nobody listed (a
  "found" allele), which is never named.
- Some tied reads write a deletion of 25 bases of their own in the CIGAR. They name nothing,
  so n drops by exactly their number.
- A one-base deletion and a listed delins (3 bases to 2) written from the next base. Reads
  that carry the delins are partial at the deletion row, which already says so. The row now
  names the listed delins: `OBSERVED_ALLELE(...:n/0)` when the deletion has no ALT reads,
  `COEXISTING_ALLELE(...:n/m)` when a few reads carry the deletion itself.

Every row in every test is checked against the four counting invariants, and the counts
(ref, alt, partial, depth, fragments) are pinned to what the build gave before the diagnostic
existed: it is count-neutral.
"""

import glob
import random

import pysam
from helpers import make_read, read_maf_output
from typer.testing import CliRunner

from gbcms.cli import app

runner = CliRunner()
READ_LEN = 100
SITE = 700  # 0-based: the SNV base, or the first base after the deleted one's anchor
BASE = "A"  # the inserted base, and the SNV's ALT base
OWN_DEL = 25  # a deletion a read writes itself: long enough to be its own allele
COUNT_COLUMNS = (
    "total_count",
    "ref_count",
    "alt_count",
    "partial_alt",
    "total_count_fragment",
    "ref_count_fragment",
    "alt_count_fragment",
)


def _contig(rng, length=1500):
    return "".join(rng.choice("ACGT") for _ in range(length))


def _run(tmp_path, ref, rows, reads):
    fa = tmp_path / "ref.fasta"
    fa.write_text(">1\n" + ref + "\n")
    pysam.faidx(str(fa))
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "1", "LN": len(ref)}]}
    raw, bam = tmp_path / "raw.bam", tmp_path / "reads.bam"
    with pysam.AlignmentFile(raw, "wb", header=header) as fh:
        for a in reads:
            fh.write(a)
    pysam.sort("-o", str(bam), str(raw))
    pysam.index(str(bam))
    vcf = tmp_path / "rows.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "".join(f"1\t{p}\t.\t{r}\t{a}\t.\t.\t.\n" for p, r, a in rows)
    )
    out = tmp_path / "out"
    args = ["dna", "-v", str(vcf), "-b", f"S:{bam}", "-f", str(fa), "-o", str(out)]
    args += ["--format", "maf", "--read-judge", "sequence"]
    result = runner.invoke(app, args, env={"COLUMNS": "400"})
    assert result.exit_code == 0, result.output
    got = list(read_maf_output(glob.glob(str(out / "*.maf"))[0]))
    for r in got:
        dp, rd, ad = int(r["total_count"]), int(r["ref_count"]), int(r["alt_count"])
        assert dp >= rd + ad
        dpf, rdf, adf = (
            int(r["total_count_fragment"]),
            int(r["ref_count_fragment"]),
            int(r["alt_count_fragment"]),
        )
        assert dpf >= rdf + adf
        assert rd == int(r["ref_count_forward"]) + int(r["ref_count_reverse"])
        assert ad == int(r["alt_count_forward"]) + int(r["alt_count_reverse"])
    return got


def _row(rows, ref_allele, alt_allele):
    (row,) = [
        r
        for r in rows
        if r["Reference_Allele"] == ref_allele and r["Tumor_Seq_Allele2"] == alt_allele
    ]
    return row


def _counts(row):
    return tuple(int(row[c]) for c in COUNT_COLUMNS)


def _tokens(row, name):
    return [t for t in row["gbcms_diagnostic"].split(";") if t.startswith(name)]


# ── An insertion of one base beside a listed SNV ──────────────────────────────────────


def _ins_snv_scene(seed=21):
    """A contig where the inserted base differs from the bases around it (so the insertion
    is in its one canonical form), the insertion row, and the SNV row at the next base."""
    rng = random.Random(seed)
    p = SITE
    while True:
        ref = _contig(rng)
        if all(ref[p + k] != BASE for k in (-2, -1, 0, 1)):
            break
    ins = (p, ref[p - 1], ref[p - 1] + BASE)  # VCF: anchor base, then the inserted base
    snv = (p + 1, ref[p], BASE)
    return ref, ins, snv


def _tied_reads(ref, n, tag, own_deletion=0):
    """Reads that end right after the base: the inserted base after ref[:SITE] and the SNV's
    ALT base at SITE read the same, so the read fits both. With `own_deletion`, the read
    writes a deletion of that many bases (away from the site) in its CIGAR."""
    p, reads = SITE, []
    for i in range(n):
        if own_deletion:
            s = p - 90 + i
            a, b = p - 55, p - 55 + own_deletion
            seq = ref[s:a] + ref[b:p] + BASE
            cigar = ((0, a - s), (2, own_deletion), (0, p - b + 1))
        else:
            s = p - 60 + i
            seq = ref[s:p] + BASE
            cigar = ((0, len(seq)),)
        reads.append(make_read(f"{tag}{i}", seq, s, cigar))
    return reads


def _carriers(ref, n_ins=10, n_snv=12, n_ref=10):
    """Reads that run well past the base, so they tell the insertion, the SNV and the
    reference apart."""
    p, reads = SITE, []
    ins_hap = ref[:p] + BASE + ref[p:]
    for i in range(n_ins):
        s = p - 60 + i
        cigar = ((0, p - s), (1, 1), (0, READ_LEN - (p - s) - 1))
        reads.append(make_read(f"i{i}", ins_hap[s : s + READ_LEN], s, cigar))
    snv_hap = ref[:p] + BASE + ref[p + 1 :]
    for i in range(n_snv):
        s = p - 40 + i
        reads.append(make_read(f"s{i}", snv_hap[s : s + READ_LEN], s, ((0, READ_LEN),)))
    for i in range(n_ref):
        s = p - 50 + i
        reads.append(make_read(f"r{i}", ref[s : s + READ_LEN], s, ((0, READ_LEN),)))
    return reads


def _ins_and_snv_rows(rows, ref):
    return _row(rows, "-", BASE), _row(rows, ref[SITE], BASE)


def test_an_insertion_and_a_listed_snv_beside_it_name_each_other(tmp_path):
    """Both rows listed. 8 reads end right after the base, so they fit the insertion and the
    SNV alike; 10 reads carry the insertion, 12 the SNV, 10 are reference. The SNV row names
    the insertion (its VCF form), and the insertion row names the SNV, 8 reads each. The m
    is each row's own alt_count: 12 for the SNV row, 10 for the insertion row. The 8 are held
    back at both rows, and no count is moved."""
    ref, ins, snv = _ins_snv_scene()
    reads = _tied_reads(ref, 8, "t") + _carriers(ref)
    rows = _run(tmp_path, ref, [ins, snv], reads)
    ins_row, snv_row = _ins_and_snv_rows(rows, ref)
    # Count neutrality: (dp, ref, alt, partial, dp fragments, ref fragments, alt fragments).
    # The 8 tied reads are in dp and in neither ref nor alt (40 - 22 - 10, and 40 - 20 - 12).
    assert _counts(ins_row) == (40, 22, 10, 0, 40, 22, 10)
    assert _counts(snv_row) == (40, 20, 12, 0, 40, 20, 12)
    assert _tokens(snv_row, "TIED_WITH") == [f"TIED_WITH(1:{ins[0]}:{ins[1]}>{ins[2]}:8/12)"]
    assert _tokens(ins_row, "TIED_WITH") == [f"TIED_WITH(1:{snv[0]}:{snv[1]}>{snv[2]}:8/10)"]


def test_a_tie_with_an_allele_nobody_listed_is_not_named(tmp_path):
    """Guard. The same reads, but only the SNV is listed: the insertion is an allele the
    input does not hold (a found allele), so the SNV row names no tie. This passes today."""
    ref, _, snv = _ins_snv_scene()
    reads = _tied_reads(ref, 8, "t") + _carriers(ref)
    rows = _run(tmp_path, ref, [snv], reads)
    (snv_row,) = rows
    assert _counts(snv_row) == (40, 20, 12, 0, 40, 20, 12)
    assert _tokens(snv_row, "TIED_WITH") == []


def test_reads_with_a_deletion_of_their_own_name_nothing(tmp_path):
    """The first shape again, but 3 of the 8 tied reads write a 25-base deletion of their
    own in the CIGAR (away from the site). They are still held back, but they name no row,
    so n is 5, not 8. The counts equal the first shape's."""
    ref, ins, snv = _ins_snv_scene()
    reads = _tied_reads(ref, 5, "t") + _tied_reads(ref, 3, "u", own_deletion=OWN_DEL)
    reads += _carriers(ref)
    rows = _run(tmp_path, ref, [ins, snv], reads)
    ins_row, snv_row = _ins_and_snv_rows(rows, ref)
    assert _counts(ins_row) == (40, 22, 10, 0, 40, 22, 10)
    assert _counts(snv_row) == (40, 20, 12, 0, 40, 20, 12)
    assert _tokens(snv_row, "TIED_WITH") == [f"TIED_WITH(1:{ins[0]}:{ins[1]}>{ins[2]}:5/12)"]
    assert _tokens(ins_row, "TIED_WITH") == [f"TIED_WITH(1:{snv[0]}:{snv[1]}>{snv[2]}:5/10)"]


# ── A one-base deletion beside a listed delins ────────────────────────────────────────


def _del_delins_scene(seed=5):
    """A one-base deletion (anchor at SITE-1, deleted base at SITE) and a listed delins of
    3 bases to 2 written from the next base. The deleted base differs from both neighbours
    (one canonical form), and the delins shares no first or last base with what it writes."""
    rng = random.Random(seed)
    d = SITE
    while True:
        ref = _contig(rng)
        if len({ref[d - 1], ref[d], ref[d + 1]}) == 3:
            break
    x = next(c for c in "ACGT" if c not in (ref[d + 1], ref[d + 2]))
    y = next(c for c in "ACGT" if c not in (ref[d + 3], ref[d + 2]))
    deletion = (d, ref[d - 1 : d + 1], ref[d - 1])
    delins = (d + 2, ref[d + 1 : d + 4], x + y)
    return ref, deletion, delins


def _delins_reads(ref, delins, n):
    """Reads that carry the delins, written as 2 inserted bases and 3 deleted ones."""
    d = SITE
    hap = ref[: d + 1] + delins[2] + ref[d + 4 :]
    reads = []
    for i in range(n):
        s = d - 60 + i
        left = d + 1 - s
        cigar = ((0, left), (1, 2), (2, 3), (0, READ_LEN - left - 2))
        reads.append(make_read(f"x{i}", hap[s : s + READ_LEN], s, cigar))
    return reads


def _deletion_reads(ref, n):
    d = SITE
    hap = ref[:d] + ref[d + 1 :]
    reads = []
    for i in range(n):
        s = d - 50 + i
        cigar = ((0, d - s), (2, 1), (0, READ_LEN - (d - s)))
        reads.append(make_read(f"d{i}", hap[s : s + READ_LEN], s, cigar))
    return reads


def _reference_reads(ref, n):
    d = SITE
    return [
        make_read(f"r{i}", ref[d - 50 + i : d + 50 + i], d - 50 + i, ((0, READ_LEN),))
        for i in range(n)
    ]


def test_a_deletion_names_a_listed_delins_its_reads_carry_with_no_alt_reads(tmp_path):
    """8 reads carry the delins, 10 are reference, none carry the deletion. The deletion row
    has no ALT reads and holds the 8 as partial (counted today). It now names the delins,
    which another input row lists, as observed: n is 8 and m is 0. The delins row counts its
    8 as ALT and names nothing."""
    ref, deletion, delins = _del_delins_scene()
    reads = _delins_reads(ref, delins, 8) + _reference_reads(ref, 10)
    rows = _run(tmp_path, ref, [deletion, delins], reads)
    del_row = _row(rows, ref[SITE], "-")
    delins_row = _row(rows, delins[1], delins[2])
    assert _counts(del_row) == (18, 10, 0, 8, 18, 10, 0)
    assert _counts(delins_row) == (18, 10, 8, 0, 18, 10, 8)
    want = f"OBSERVED_ALLELE(1:{delins[0]}:{delins[1]}>{delins[2]}:8/0)"
    assert _tokens(del_row, "OBSERVED_ALLELE") == [want]
    assert _tokens(del_row, "COEXISTING_ALLELE") == []
    assert _tokens(delins_row, "OBSERVED_ALLELE") == []


def test_a_deletion_with_a_few_alt_reads_names_the_listed_delins_as_coexisting(tmp_path):
    """The same, with 3 reads that carry the deletion itself. The deletion row has AD 3, and
    8 reads carry the listed delins (n above AD and at least 3), so the name is
    COEXISTING_ALLELE with m 3. The 8 are partial; the 3 are the deletion's ALT."""
    ref, deletion, delins = _del_delins_scene()
    reads = _delins_reads(ref, delins, 8) + _deletion_reads(ref, 3) + _reference_reads(ref, 10)
    rows = _run(tmp_path, ref, [deletion, delins], reads)
    del_row = _row(rows, ref[SITE], "-")
    delins_row = _row(rows, delins[1], delins[2])
    assert _counts(del_row) == (21, 10, 3, 8, 21, 10, 3)
    assert _counts(delins_row) == (21, 10, 8, 3, 21, 10, 8)
    want = f"COEXISTING_ALLELE(1:{delins[0]}:{delins[1]}>{delins[2]}:8/3)"
    assert _tokens(del_row, "COEXISTING_ALLELE") == [want]
    assert _tokens(del_row, "OBSERVED_ALLELE") == []


def test_a_few_delins_reads_are_not_named_at_the_deletion(tmp_path):
    """Gate. Only 2 reads carry the listed delins (below the 3 a name needs): the deletion
    row, which does not count them as partial either, and names no OBSERVED_ALLELE."""
    ref, deletion, delins = _del_delins_scene()
    reads = _delins_reads(ref, delins, 2) + _reference_reads(ref, 10)
    rows = _run(tmp_path, ref, [deletion, delins], reads)
    del_row = _row(rows, ref[SITE], "-")
    assert _counts(del_row) == (12, 10, 0, 0, 12, 10, 0)
    assert _tokens(del_row, "OBSERVED_ALLELE") == []
    assert _tokens(del_row, "COEXISTING_ALLELE") == []


def test_a_delins_no_more_common_than_the_deletion_is_not_named_as_coexisting(tmp_path):
    """Gate. 5 reads carry the delins and 5 carry the deletion: n equals the deletion row's
    AD, not above it, so the deletion row names no COEXISTING_ALLELE."""
    ref, deletion, delins = _del_delins_scene()
    reads = _delins_reads(ref, delins, 5) + _deletion_reads(ref, 5) + _reference_reads(ref, 10)
    rows = _run(tmp_path, ref, [deletion, delins], reads)
    del_row = _row(rows, ref[SITE], "-")
    assert _counts(del_row) == (20, 10, 5, 5, 20, 10, 5)
    assert _tokens(del_row, "COEXISTING_ALLELE") == []
    assert _tokens(del_row, "OBSERVED_ALLELE") == []


def test_a_tie_is_named_past_sixty_four_rows_within_reach(tmp_path):
    """The first shape with 70 other listed SNVs (nobody's reads carry them), 4 bases apart
    and ending 12 bases before the insertion, all within 300 bases of the pair. Most of the
    rows in reach come before the pair, so each partner sits past the 64th row in reach. The
    two rows still name each other, 8 reads each, and the counts equal the first shape's."""
    ref, ins, snv = _ins_snv_scene()
    fillers = []
    for i in range(70):
        q = SITE - 12 - 4 * (69 - i)  # 0-based, ascending, the last 12 before the insertion
        fillers.append((q + 1, ref[q], "C" if ref[q] != "C" else "G"))
    reads = _tied_reads(ref, 8, "t") + _carriers(ref)
    rows = _run(tmp_path, ref, fillers + [ins, snv], reads)
    ins_row, snv_row = _ins_and_snv_rows(rows, ref)
    assert _counts(ins_row) == (40, 22, 10, 0, 40, 22, 10)
    assert _counts(snv_row) == (40, 20, 12, 0, 40, 20, 12)
    assert _tokens(snv_row, "TIED_WITH") == [f"TIED_WITH(1:{ins[0]}:{ins[1]}>{ins[2]}:8/12)"]
    assert _tokens(ins_row, "TIED_WITH") == [f"TIED_WITH(1:{snv[0]}:{snv[1]}>{snv[2]}:8/10)"]


# ── A long deletion beside a listed SNV ───────────────────────────────────────────────


def _long_deletion_scene(seed=31):
    """A listed deletion of 25 bases removing [SITE, SITE+25) (VCF: anchored at SITE-1), and a
    listed SNV at SITE whose ALT base is the first base after the deleted span. A read that
    ends right after that base fits both. The deletion is in its one canonical form: its
    anchor differs from its last deleted base."""
    rng = random.Random(seed)
    q, n = SITE, OWN_DEL
    while True:
        ref = _contig(rng)
        if ref[q] != ref[q + n] and ref[q - 1] != ref[q + n - 1]:
            break
    deletion = (q, ref[q - 1 : q + n], ref[q - 1])
    snv = (q + 1, ref[q], ref[q + n])
    return ref, deletion, snv


def test_a_tied_read_that_writes_the_named_deletion_still_names_it(tmp_path):
    """8 reads end right after the SNV's ALT base: 5 written as plain bases with a mismatch
    at the SNV, 3 written with the 25-base deletion in the CIGAR. The deletion IS the row
    they name, so each carries the allele it names and all 8 count. Besides them: 10 deletion
    carriers, 12 SNV carriers and 10 reference reads. All 8 are partial at both rows (the
    other allele), in dp and in neither ref nor alt. Both rows name the other, 8 reads each."""
    ref, deletion, snv = _long_deletion_scene()
    q, n, alt = SITE, OWN_DEL, snv[2]
    reads = []
    for i in range(8):
        s = q - 60 + i
        seq = ref[s:q] + alt
        cigar = ((0, q - s), (2, n), (0, 1)) if i >= 5 else ((0, len(seq)),)
        reads.append(make_read(f"t{i}", seq, s, cigar))
    del_hap = ref[:q] + ref[q + n :]
    for i in range(10):
        s = q - 60 + i
        cigar = ((0, q - s), (2, n), (0, READ_LEN - (q - s)))
        reads.append(make_read(f"l{i}", del_hap[s : s + READ_LEN], s, cigar))
    snv_hap = ref[:q] + alt + ref[q + 1 :]
    for i in range(12):
        s = q - 40 + i
        reads.append(make_read(f"s{i}", snv_hap[s : s + READ_LEN], s, ((0, READ_LEN),)))
    for i in range(10):
        s = q - 50 + i
        reads.append(make_read(f"r{i}", ref[s : s + READ_LEN], s, ((0, READ_LEN),)))
    rows = _run(tmp_path, ref, [deletion, snv], reads)
    del_row, snv_row = _row(rows, ref[q : q + n], "-"), _row(rows, ref[q], alt)
    assert _counts(del_row) == (40, 10, 10, 20, 40, 10, 10)
    assert _counts(snv_row) == (40, 10, 12, 18, 40, 10, 12)
    want_del = f"TIED_WITH(1:{deletion[0]}:{deletion[1]}>{deletion[2]}:8/12)"
    want_snv = f"TIED_WITH(1:{snv[0]}:{snv[1]}>{snv[2]}:8/10)"
    assert _tokens(snv_row, "TIED_WITH") == [want_del]
    assert _tokens(del_row, "TIED_WITH") == [want_snv]
