"""`rna_antisense_depth` counts antisense reads at RNA defaults (#114).

With strandedness enforced (the RNA default), antisense reads were dropped
before the sense/antisense tally, so `rna_antisense_depth` was always 0 although
IGV shows the reads (0.5% of reads at 18% of the FORTE truth loci). An antisense
read is now classified as a sense read would be and tallied there when it is a
first-class REF or ALT read over the anchor, then dropped: REF, ALT, depth,
fragments and every other count are unchanged. The column means the same thing
with `--no-strandedness`, where those reads are also counted.

`STRAND_DISCORDANT` stays a `--no-strandedness` diagnostic: under enforcement no
antisense read reaches the junction tally, so it cannot fire (documented).

Committed red (xfail-strict) before the implementation; flipped green with it.
"""

from helpers import make_read
from rna_fixtures import (
    E1,
    E2,
    E3,
    READ_LEN,
    SENSE,
    mk_ref,
    run_rna,
    spliced,
    write_bam,
    write_fasta,
    write_gtf,
    write_vcf,
)

ANTISENSE = 0  # forward single-end: the opposite transcript strand under dUTP ('reverse')
SNV = 200  # 0-based, mid-E1, ~100bp from any edge
_REF = mk_ref()
_ALT = "A" if _REF[SNV] != "A" else "C"
_THIRD = next(b for b in "ACGT" if b not in (_REF[SNV], _ALT))


def _through(base, n, prefix, flag, mapq=60, start=SNV - 50, tags=()):
    """Unspliced reads covering SNV (from `start`), carrying `base` there."""
    out = []
    for i in range(n):
        s = start + (i % 5)
        seq = _REF[s:SNV] + base + _REF[SNV + 1 : s + READ_LEN]
        r = make_read(f"{prefix}{i}", seq, s, ((0, READ_LEN),), flag=flag, mapq=mapq)
        for tag, value in tags:
            r.set_tag(tag, value)
        out.append(r)
    return out


def _run(tmp_path, reads, extra=(), name="out"):
    (row,) = run_rna(
        tmp_path,
        write_vcf(tmp_path, [(SNV + 1, _REF[SNV], _ALT)]),
        write_bam(tmp_path, _REF, reads, name=f"{name}.bam"),
        write_fasta(tmp_path, _REF),
        write_gtf(tmp_path),
        outname=name,
        extra=extra,
    )
    return row


def _counts(row):
    return tuple(int(row[c]) for c in ("ref_count", "alt_count", "total_count"))


def _depths(row):
    return tuple(
        int(row[c]) for c in ("rna_sense_depth", "rna_antisense_depth", "rna_alt_sense_count")
    )


_MIXED = (
    _through(_REF[SNV], 20, "sref", SENSE)
    + _through(_ALT, 10, "salt", SENSE)
    + _through(_REF[SNV], 3, "aref", ANTISENSE)
    + _through(_ALT, 2, "aalt", ANTISENSE)
)


def test_antisense_reads_are_tallied_at_defaults_without_changing_counts(tmp_path):
    row = _run(tmp_path, _MIXED)
    assert _counts(row) == (20, 10, 30), "antisense reads stay out of REF, ALT and depth"
    assert _depths(row) == (30, 5, 10)


def test_the_column_means_the_same_with_and_without_enforcement(tmp_path):
    enforced = _run(tmp_path, _MIXED, name="enforced")
    free = _run(tmp_path, _MIXED, extra=("--no-strandedness",), name="free")
    assert _counts(free) == (23, 12, 35), "without enforcement antisense reads are counted"
    assert _depths(enforced)[:2] == _depths(free)[:2] == (30, 5)


def test_only_antisense_reads_a_sense_read_would_count_are_tallied(tmp_path):
    """Two antisense REF reads are tallied. Not tallied, as a sense read would
    not be counted: a third allele (neither), a multi-mapper (MAPQ 0, NH:i:3),
    and a read that ends before the SNV (no anchor overlap)."""
    reads = (
        _through(_REF[SNV], 20, "sref", SENSE)
        + _through(_REF[SNV], 2, "aref", ANTISENSE)
        + _through(_THIRD, 1, "athird", ANTISENSE)
        + _through(_REF[SNV], 1, "amulti", ANTISENSE, mapq=0, tags=(("NH", 3),))
        + [
            make_read(
                "aearly", _REF[SNV - 120 : SNV - 20], SNV - 120, ((0, READ_LEN),), flag=ANTISENSE
            )
        ]
    )
    row = _run(tmp_path, reads)
    assert _counts(row) == (20, 0, 20)
    assert _depths(row)[:2] == (20, 2)


# ── STRAND_DISCORDANT is a --no-strandedness diagnostic ──────────────────────
def test_strand_discordant_needs_no_strandedness(tmp_path):
    """Guard (green before and after). REF reads splice E1->E2; ALT reads skip E2
    (E1->E3), six sense and four antisense (a 40% minority strand). Under
    enforcement the antisense junction reads are filtered, so the flag cannot
    fire; with --no-strandedness it does."""
    snv = E1[1] - 20  # covered by the spliced reads' first block
    alt = "A" if _REF[snv] != "A" else "C"
    mutant = _REF[:snv] + alt + _REF[snv + 1 :]
    reads = spliced(_REF, E1[1], E2[0], 30, "ref") + spliced(mutant, E1[1], E3[0], 6, "skip_s")
    for r in spliced(mutant, E1[1], E3[0], 4, "skip_a"):
        r.flag = ANTISENSE
        reads.append(r)
    rows = {}
    for tag, extra in (("enforced", ()), ("free", ("--no-strandedness",))):
        (rows[tag],) = run_rna(
            tmp_path,
            write_vcf(tmp_path, [(snv + 1, _REF[snv], alt)]),
            write_bam(tmp_path, _REF, reads, name=f"{tag}.bam"),
            write_fasta(tmp_path, _REF),
            write_gtf(tmp_path),
            outname=tag,
            extra=extra,
        )
    assert "STRAND_DISCORDANT" not in rows["enforced"]["asjd_diagnostic"]
    assert "STRAND_DISCORDANT" in rows["free"]["asjd_diagnostic"], rows["free"]["asjd_diagnostic"]
