"""A GTF transcript ID that another chromosome (or strand) reuses is another transcript.

The parser grouped exons by transcript ID alone, so an ID used on two chromosomes
(UCSC/RefSeq PAR copies on X and Y, alternate loci, versions that differ only by suffix)
merged into one transcript: its exons from both places were sorted together, the
consecutive pairs made introns that exist nowhere and missed the real ones, and the
whole transcript sat on the first chromosome seen. A spliced read over a real intron then
matched none of the merged introns, so it was left out of the transcript's counts.

Contracts:
- the transcript is read per ID, chromosome and strand: the per-transcript counts of a
  variant on either chromosome include the reads that splice over that place's own intron;
- the parser warns once when a GTF reuses an ID, and not otherwise.
"""

import pysam
from helpers import make_read
from rna_fixtures import E1, E2, E3, READ_LEN, SENSE, mk_ref, run_rna
from typer.testing import CliRunner

from gbcms.cli import app

runner = CliRunner()

# chr1 reads the shared three-exon gene model (introns [300,500) and [600,800)); chr2 gives
# the same ID two exons of its own, at coordinates that interleave chr1's, so the merged
# exons make introns where neither gene has one.
CHR2_EXONS = ((200, 400), (700, 900))  # one intron, [400,700)
ID = "TD"
DONOR1, ACCEPTOR1 = E1[1], E2[0]  # chr1's first intron
DONOR2, ACCEPTOR2 = CHR2_EXONS[0][1], CHR2_EXONS[1][0]
POS1 = ACCEPTOR1 + 30  # 0-based, inside E2, 30 bases past the acceptor
POS2 = ACCEPTOR2 + 30


def _gtf(tmp_path, chr2_id=ID, name="gene.gtf"):
    lines = [("chr1", s, e, ID) for s, e in (E1, E2, E3)]
    lines += [("chr2", s, e, chr2_id) for s, e in CHR2_EXONS]
    path = tmp_path / name
    path.write_text(
        "".join(
            f'{c}\tTEST\texon\t{s + 1}\t{e}\t.\t+\t.\tgene_id "G_{c}"; transcript_id "{t}";\n'
            for c, s, e, t in lines
        )
    )
    return path


def _spliced(ref, donor, acceptor, pos, base, n, prefix, tid):
    """Reads whose first block ends at `donor`, then skip to `acceptor`, and carry `base`
    at `pos` (30 bases into the second block)."""
    out = []
    for i in range(n):
        s = donor - 40 - (i % 5)
        left = donor - s
        right = READ_LEN - left
        seq = list(ref[s:donor] + ref[acceptor : acceptor + right])
        seq[left + (pos - acceptor)] = base
        cigar = ((0, left), (3, acceptor - donor), (0, right))
        r = make_read(f"{prefix}{i}", "".join(seq), s, cigar, flag=SENSE)
        r.reference_id = tid
        out.append(r)
    return out


def _files(tmp_path, ref):
    fa = tmp_path / "ref.fasta"
    fa.write_text(f">chr1\n{ref}\n>chr2\n{ref}\n")
    pysam.faidx(str(fa))
    return fa


def _bam(tmp_path, ref, reads):
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "1", "LN": len(ref)}, {"SN": "2", "LN": len(ref)}],
    }
    raw = tmp_path / "raw.bam"
    with pysam.AlignmentFile(raw, "wb", header=header) as fh:
        for a in sorted(reads, key=lambda a: (a.reference_id, a.reference_start)):
            fh.write(a)
    out = tmp_path / "s.bam"
    pysam.sort("-o", str(out), str(raw))
    pysam.index(str(out))
    return out


def _vcf(tmp_path, rows):
    vcf = tmp_path / "variants.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "".join(f"{c}\t{p + 1}\t.\t{r}\t{a}\t.\t.\t.\n" for c, p, r, a in rows)
    )
    return vcf


def _other(base):
    return next(b for b in "ACGT" if b != base)


def _setup(tmp_path, chr2_id=ID):
    """A variant on each chromosome: chr1 has 4 ALT and 6 REF reads spliced over chr1's
    intron, chr2 has 3 ALT and 5 REF reads spliced over chr2's."""
    ref = mk_ref()
    a1, a2 = _other(ref[POS1]), _other(ref[POS2])
    reads = (
        _spliced(ref, DONOR1, ACCEPTOR1, POS1, a1, 4, "a1_", 0)
        + _spliced(ref, DONOR1, ACCEPTOR1, POS1, ref[POS1], 6, "r1_", 0)
        + _spliced(ref, DONOR2, ACCEPTOR2, POS2, a2, 3, "a2_", 1)
        + _spliced(ref, DONOR2, ACCEPTOR2, POS2, ref[POS2], 5, "r2_", 1)
    )
    variants = _vcf(tmp_path, [("chr1", POS1, ref[POS1], a1), ("chr2", POS2, ref[POS2], a2)])
    return variants, _bam(tmp_path, ref, reads), _files(tmp_path, ref), _gtf(tmp_path, chr2_id)


def test_each_chromosomes_variant_counts_the_reads_spliced_over_its_own_intron(tmp_path):
    """ID TD is a three-exon transcript on chr1 and a two-exon transcript on chr2. A
    read spliced over chr1's intron is compatible with the chr1 transcript (it was
    left out: the merged introns held none of chr1's), and so on chr2. The main counts
    agree with the per-transcript ones, as each place has one transcript."""
    variants, bam, fa, gtf = _setup(tmp_path)
    r1, r2 = run_rna(tmp_path, variants, bam, fa, gtf)
    assert (r1["alt_count"], r1["ref_count"], r1["total_count"]) == ("4", "6", "10")
    assert (r2["alt_count"], r2["ref_count"], r2["total_count"]) == ("3", "5", "8")
    assert r1["transcript_read_counts"] == f"{ID}:4,6,10"
    assert r2["transcript_read_counts"] == f"{ID}:3,5,8"
    assert r1["transcript_fragment_counts"] == f"{ID}:4,6,10"
    assert r2["transcript_fragment_counts"] == f"{ID}:3,5,8"


def test_distinct_ids_give_the_same_counts(tmp_path):
    """Guard: the same data with a different ID on chr2 (no reuse) counts the same, so
    keying by place reads each chromosome as if its ID were unique."""
    variants, bam, fa, gtf = _setup(tmp_path, chr2_id="TE")
    r1, r2 = run_rna(tmp_path, variants, bam, fa, gtf)
    assert r1["transcript_read_counts"] == f"{ID}:4,6,10"
    assert r2["transcript_read_counts"] == "TE:3,5,8"


def _log(tmp_path, chr2_id=ID):
    """The log text of a run over `_setup`'s files (unwrapped: Rich wraps log lines at the
    terminal width)."""
    variants, bam, fa, gtf = _setup(tmp_path, chr2_id)
    out = tmp_path / "again"
    out.mkdir()
    args = ["rna", "-v", str(variants), "-b", f"S:{bam}", "-f", str(fa), "-o", str(out)]
    result = runner.invoke(app, [*args, "--gtf", str(gtf)], env={"COLUMNS": "3000"})
    assert result.exit_code == 0, result.output
    return result.output


def test_a_gtf_that_reuses_an_id_is_warned_about_once(tmp_path):
    hits = [m for m in _log(tmp_path).splitlines() if "more than one chromosome" in m]
    assert len(hits) == 1, hits
    assert "1 transcript ID" in hits[0] and ID in hits[0], hits


def test_a_gtf_without_reuse_is_not_warned_about(tmp_path):
    assert "more than one chromosome" not in _log(tmp_path, chr2_id="TE")
