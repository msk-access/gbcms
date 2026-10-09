"""The read-judgment cases: the executable table behind
`docs/reference/read-judgment.md`.

Each case is one read shape at one variant: four reads of that shape, counted
alone, and the call they get as (ref_count, alt_count, partial_alt). Every case
carries the decision it follows. A decided case states the decided call; an open
case pins today's call until its decision lands, so a change to either shows in
review (`tests/test_read_judgment_spec.py`). Changing an expectation is a
decision: it needs the operator's say first (AGENTS.md).

Synthetic and PHI-free. `python tests/read_judgment_cases.py` prints the table,
with each read's own verdict where the read census can judge it (pure indels).
"""

from __future__ import annotations

import random
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from census import census, tract  # noqa: E402
from helpers import make_read, write_contig  # noqa: E402

from gbcms import _rs  # noqa: E402

A, L, N_READS = 400, 900, 4
ARGS = {
    "min_mapq": 20,
    "min_baseq": 20,
    "filter_duplicates": True,
    "filter_secondary": True,
    "filter_supplementary": True,
    "filter_qc_failed": False,
    "filter_improper_pair": False,
    "filter_indel": False,
    "threads": 1,
}

# Pure indels at A (motif written from the anchor), as in the read census.
PURE = {
    "hp-A": ("GAAAAAT", "GA", "G"),
    "hp-AA": ("GAAAAAAT", "GAA", "G"),
    "hp+A": ("GAAAAAT", "G", "GA"),
    "u-8": ("GATCGGATAC", "GATCGGATA", "G"),
    "u+10": ("GT", "G", "GACGTTGCATC"),
    "ca-CA": ("G" + "CA" * 5 + "T", "GCA", "G"),
    "dup-8": ("G" + "ACGTTGCA" * 2 + "T", "GACGTTGCA", "G"),
    "dup+8": ("G" + "ACGTTGCA" * 2 + "T", "G", "GACGTTGCA"),
}
# Anchor-changing events before an A run: a whole-window event (10 A's) and a
# long one judged by junction windows (60 A's).
COMPLEX = {
    "C>TA run10": (10, "C", "TA"),
    "C>TA run60": (60, "C", "TA"),
    "CA>T run60": (60, "CA", "T"),
}

# Delins in unique sequence, read by one-sided junction windows (RJ-23): one whose
# windows are whole (28 reference bases replaced by 6) and a long one (60 by 2).
DELINS = {
    "delins 28>6": (28, "ATATGA"),
    "delins 60>2": (60, "AT"),
}
# How far past the event a recurring change sits on the ALT reads of the guard case.
RECUR_AT = 7


# What a read contributes, at an SNV (C>A at A + 1) and an unprepared Del+SNV
# (GCT>A at A, which the previous complex classifier judges): read-through bases
# past the fragment end, absent base qualities, hard clips.
# The Illumina TruSeq adapter a read runs into past a short insert.
ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" * 2

READ_INPUTS = {
    "read-through: the base on the SNV is adapter, showing the ALT",
    "read-through: the base on the SNV is adapter, showing REF",
    "read-through: the read aligns on past its TLEN end (R2's 5' end clipped), showing the ALT",
    "absent base qualities (QUAL '*'), showing the ALT",
    "previous complex classifier: REF read, no clip",
    "previous complex classifier: REF read with a leading hard clip",
}


@dataclass(frozen=True)
class Case:
    id: str
    group: str
    variant: str
    shape: str


def _contig(motif: str, seed: int = 189) -> str:
    rng = random.Random(seed)
    seq = [rng.choice("CGT") for _ in range(L)]
    seq[A : A + len(motif)] = list(motif)
    return "".join(seq)


def _complex_contig(run: int) -> str:
    rng = random.Random(61)
    left = "".join(rng.choice("CGT") for _ in range(A))
    right = "".join(rng.choice("CGT") for _ in range(L - A - run - 2))
    return left + "C" + "A" * run + "G" + right


def _build(contig: str, start: int, length: int, events) -> tuple[str, tuple]:
    """A read of `length` bases from `start`, with events (pos, "I", bases) or
    (pos, "D", n) placed before their reference positions."""
    seq, cigar, p = "", [], start
    for pos, op, x in sorted(events, key=lambda e: e[0]):
        if pos > p:
            seq += contig[p:pos]
            cigar.append([0, pos - p])
            p = pos
        if op == "I":
            seq += x
            cigar.append([1, len(x)])
        elif op == "X":  # a substitution: one base
            seq += x
            cigar.append([0, 1])
            p += 1
        else:
            cigar.append([2, x])
            p += x
    seq += contig[p : p + length]
    cigar.append([0, length])
    out, q = [], 0
    for op, n in cigar:
        if q >= length:
            break
        if op == 2:
            out.append([2, n])
            continue
        take = min(n, length - q)
        out.append([op, take])
        q += take
    merged: list[list[int]] = []
    for op, n in out:
        if merged and merged[-1][0] == op == 0:
            merged[-1][1] += n
        else:
            merged.append([op, n])
    return seq[:length], tuple(map(tuple, merged))


def _own(ref: str, alt: str):
    return (A + 1, "D", len(ref) - 1) if len(ref) > len(alt) else (A + 1, "I", alt[1:])


def _pure_shapes(contig: str, ref: str, alt: str):
    """Read shapes at a pure indel: (events, read start offset, read length)."""
    lo, hi = tract(A, ref, alt, contig)
    own = _own(ref, alt)
    other = "C" if contig[hi + 4] != "C" else "G"
    in_d1 = (lo + 1, "D", 1) if hi - lo > 1 else (A + 2, "D", 1)
    # The discrimination window: the tract plus one base each side.
    w_lo, w_hi = lo - 1, hi + 1

    def where(event) -> str:
        pos, op, x = event
        inside = (w_lo < pos < w_hi) if op == "I" else (pos < w_hi and pos + x > w_lo)
        return "inside" if inside else "outside"

    def is_the_alt(event) -> bool:
        """The variant itself at another placement in its tract."""
        pos, op, x = event
        if op == "D":
            return len(ref) > len(alt) and x == len(ref) - 1 and lo <= pos and pos + x <= hi
        ins = alt[1:]
        return len(alt) > len(ref) and len(x) == len(ins) and x in ins + ins and lo <= pos <= hi

    shapes = {
        "REF read spanning the tract": ([], 0, 100),
        "REF read ending inside the tract": ([], 0, (lo + 1) - (A - 50)),
        "exact carrier": ([own], 0, 100),
        "carrier with another indel outside the window": ([own, (hi + 6, "D", 1)], 0, 100),
    }
    for label, ev in (
        ("D1 near the anchor", in_d1),
        ("I1 near the anchor", (lo + 1, "I", contig[lo + 1])),
        ("D2 after the anchor", (A + 2, "D", 2)),
        ("D1 past the tract", (hi + 4, "D", 1)),
        ("D5 past the tract", (hi + 4, "D", 5)),
        ("I1 past the tract", (hi + 4, "I", other)),
    ):
        if is_the_alt(ev):
            shapes[f"the ALT at another placement: {label}"] = ([ev], 0, 100)
        else:
            shapes[f"another indel {where(ev)} the window: {label}"] = ([ev], 0, 100)
    if len(ref) > 2:  # a deletion of two or more bases written as two ops
        k = len(ref) - 1
        shapes["the ALT written as two deletions"] = (
            [(A + 1, "D", 1), (A + 2, "D", k - 1)],
            0,
            100,
        )
    run_base = contig[A + 1]
    in_run = contig[A + 1 : A + 6] == run_base * 5
    if in_run and len(ref) == 2 and len(alt) == 1:  # -1 in a run as a shifted pair
        shapes["the ALT written as D2 + I1 at the anchor"] = (
            [(A + 1, "D", 2), (A + 3, "I", run_base)],
            0,
            100,
        )
        shapes["the ALT written as D2 + I1 inside the run"] = (
            [(A + 2, "D", 2), (A + 4, "I", run_base)],
            0,
            100,
        )
    if in_run and len(ref) == 1 and len(alt) == 2:  # +1 in a run as a shifted pair
        shapes["the ALT written as I2 + D1 at the anchor"] = (
            [(A + 1, "I", run_base * 2), (A + 1, "D", 1)],
            0,
            100,
        )
    if len(ref) > 3:  # the anchor deleted
        shapes["anchor deleted, then an insertion"] = ([(A, "D", 3), (A + 5, "I", "T")], 0, 100)
        shapes["anchor deleted"] = ([(A, "D", 3)], 0, 100)
    if len(alt) > 2:  # an insertion row whose read deletes the anchor
        shapes["anchor deleted"] = ([(A, "D", 2)], 0, 100)
    return shapes


def _unreadable_shapes(contig: str, ref: str, alt: str):
    """Carriers whose inserted bases are masked: (events, (kind, first, count)), the
    read's inserted bases [first, first + count) masked as N at Q2 ("N") or kept at
    Q5 ("low")."""
    ins = alt[1:]
    if len(ins) == 1:  # +1 in a run: the aligner writes the masked base as the insertion
        return {
            "an unreadable insert at the anchor (N)": ([(A + 1, "I", "N")], ("N", 0, 1)),
            "an unreadable insert inside the run (N)": ([(A + 3, "I", "N")], ("N", 0, 1)),
            "an insert below min BQ inside the run": ([(A + 3, "I", ins)], ("low", 0, 1)),
        }
    k = len(ins)
    return {
        "an unreadable insert (N)": ([(A + 1, "I", ins)], ("N", 0, k)),
        "an insert below min BQ": ([(A + 1, "I", ins)], ("low", 0, k)),
        "a partly masked insert": ([(A + 1, "I", ins)], ("N", 0, k // 2)),
        "anchor deleted, the insert unreadable": (
            [(A, "D", 1), (A + 1, "I", contig[A] + ins)],
            ("N", 1, k),
        ),
    }


def _other_bases_shapes(contig: str, ref: str, alt: str):
    """Same-length inserts of other bases near a pure insertion, and the ALT written
    one junction off with a compensating mismatch: read-judgment events."""
    ins = alt[1:]
    other = "CTTAGCCTAG"[: len(ins)]
    if len(ins) == 8:  # dup+8: inside its duplicated tract
        return {"a same-length insert of other bases inside the tract": [(A + 3, "I", other)]}
    return {
        "a same-length insert of other bases one junction left": [(A, "I", other)],
        "a same-length insert of other bases four junctions right": [(A + 5, "I", other)],
        "the ALT written one junction off with a compensating mismatch": [
            (A + 1, "X", ins[0]),
            (A + 2, "I", ins[1:] + contig[A + 1]),
        ],
        # Its first inserted base misread as the reference base after the anchor:
        # the aligner writes it one junction right with no mismatch; it slides back.
        "the ALT with an error at its first inserted base, written one junction right": [
            (A + 2, "I", ins[1:] + contig[A + 1]),
        ],
    }


def _behind_substitution_shapes(contig: str, ref: str, alt: str):
    """The ALT with a substitution next to its insert: at the junction with a
    substitution in the run (the strict path keeps it), and with the anchor
    substituted, written one junction left (no slide carries it back): events."""
    ins = alt[1:]
    out = {
        # The aligner's tie: one mismatch either way, at the anchor; the first
        # inserted base (the read's anchor, then the ALT's first bases) goes left.
        "the ALT with its anchor substituted, written one junction left": [
            (A, "I", "T" + ins[:-1]),
            (A, "X", ins[-1]),
        ],
    }
    if len(set(ins)) == 1:  # in a run: a substitution the slide would absorb
        out["the ALT at the junction with a substitution in the run"] = [
            (A + 1, "I", ins),
            (A + 3, "X", "C"),
        ]
    return out


def _mask_inserted(read, mask) -> None:
    """Mask the read's inserted bases [first, first + count), in CIGAR order."""
    kind, first, count = mask
    seq, quals = list(read.query_sequence), list(read.query_qualities)
    q, k = 0, 0
    for op, n in read.cigartuples:
        if op == 1:
            for i in range(q, q + n):
                if first <= k < first + count:
                    seq[i], quals[i] = ("N", 2) if kind == "N" else (seq[i], 5)
                k += 1
        if op in (0, 1, 4):
            q += n
    read.query_sequence = "".join(seq)
    read.query_qualities = quals


def _complex_reads(run: int, ref_allele: str, alt_allele: str):
    """Read haplotypes at an anchor-changing event before an A run, ending inside
    the run, on the base after it, or past it: (events, last reference base)."""
    end_run = A + 1 + run  # the G after the run
    haps = {
        "REF": [],
        "substitution only": [(A, "X", "T")],
        "anchor kept, one A more": [(A + 1, "I", "A")],
        "anchor kept, one A fewer": [(A + 1, "D", 1)],
        "exact carrier": [(A, "X", "T")]
        + ([(A + 1, "I", "A")] if alt_allele == "TA" else [(A + 1, "D", 1)]),
    }
    ends = {
        "ends inside the run": A + run // 2,
        "ends on the base after the run": end_run,
        "ends past the run": end_run + 20,
    }
    return haps, ends


def _delins_reads(n: int, alt: str, shape: str):
    """Reads (start, query length, events) of a delins shape: the delins replaces
    reference [A, A + n) by `alt`; starts and ends vary per read so no window length
    is special."""
    carrier = [(A, "I", alt), (A, "D", n)]
    m = len(alt)
    out = []
    for i in range(N_READS):
        far = 60 + 3 * i
        if shape == "exact carrier, eight flank bases right":
            # eight: the 28-to-6 event grows three bases right through a TCG repeat
            out.append((A - far, far + m + 8, carrier))
        elif shape == "exact carrier, four flank bases left":
            out.append((A - 4, 4 + m + far, carrier))
        elif shape == "exact carrier, long flanks":
            out.append((A - far, far + m + 60, carrier))
        elif shape == "REF read, four flank bases left":
            out.append((A - 4, 4 + n + far, []))
        elif shape == "REF read, four flank bases right":
            out.append((A - far, far + n + 4, []))
        elif shape == "REF read starting inside the event":
            out.append((A + n // 2, n - n // 2 + far, []))
        elif shape == "carriers ending before a change that recurs on ALT reads only":
            # four carriers end four bases past the event (inside its right flank, as
            # it grows three bases through a TCG repeat: depth only); four more read
            # on past a change seven bases out, which every one of them shows
            out.append((A - far, far + m + 4, carrier))
            sub = (A + n + RECUR_AT - 1, "X", "A")
            out.append((A - far, far + m + 40, carrier + [sub]))
        else:
            raise KeyError(shape)
    return out


DELINS_SHAPES = {
    "delins 28>6": (
        "exact carrier, eight flank bases right",
        "exact carrier, four flank bases left",
        "exact carrier, long flanks",
        "REF read, four flank bases left",
        "REF read, four flank bases right",
        "REF read starting inside the event",
        "carriers ending before a change that recurs on ALT reads only",
    ),
    "delins 60>2": (
        "exact carrier, eight flank bases right",
        "exact carrier, long flanks",
        "REF read, four flank bases left",
        "REF read starting inside the event",
    ),
}


def cases() -> list[Case]:
    out: list[Case] = []
    for v, (motif, ref, alt) in PURE.items():
        for shape in _pure_shapes(_contig(motif), ref, alt):
            group = _pure_group(shape)
            out.append(Case(f"{v} | {shape}", group, v, shape))
    for v, (run, r, a) in COMPLEX.items():
        haps, ends = _complex_reads(run, r, a)
        for h in haps:
            for e in ends:
                group = "C25 long events" if run > 40 else "C1 exact carriers"
                out.append(Case(f"{v} | {h}, {e}", group, v, f"{h}, {e}"))
    for v, shapes in DELINS_SHAPES.items():
        for shape in shapes:
            out.append(Case(f"{v} | {shape}", "C39 one-sided windows", v, shape))
    out.append(
        Case(
            "u-8 | sibling SNV inside the span", "C2 siblings", "u-8", "sibling SNV inside the span"
        )
    )
    out.append(
        Case(
            "u-8 | sibling SNV outside the window",
            "C2 siblings",
            "u-8",
            "sibling SNV outside the window",
        )
    )
    for shape in sorted(READ_INPUTS):
        out.append(Case(f"read inputs | {shape}", "read inputs", "read inputs", shape))
    for v in ("hp+A", "u+10"):
        motif, ref, alt = PURE[v]
        for shape in _unreadable_shapes(_contig(motif), ref, alt):
            out.append(Case(f"{v} | {shape}", "C35 unreadable inserts", v, shape))
    for v in ("u+10", "dup+8"):
        motif, ref, alt = PURE[v]
        for shape in _other_bases_shapes(_contig(motif), ref, alt):
            out.append(Case(f"{v} | {shape}", "C36 inserts of other bases", v, shape))
    for v in ("hp+A", "u+10"):
        motif, ref, alt = PURE[v]
        for shape in _behind_substitution_shapes(_contig(motif), ref, alt):
            out.append(Case(f"{v} | {shape}", "C37 inserts behind a substitution", v, shape))
    return out


def _read_input_bam(d: Path, shape: str):
    """The reads of a read-input case, the contig, and the variant they cover."""
    contig = _contig("GCT")
    reads = []
    if shape.startswith("read-through"):
        # R1's molecule ends just before the SNV and it reads on into adapter: the
        # adapter's first base is aligned on the SNV by chance (showing the ALT or
        # REF), the rest soft-clipped. With "aligns on", its bases past TLEN's end
        # are the genome's instead (R2's 5' end clipped): the molecule goes on.
        base = "A" if "ALT" in shape else contig[A + 1]
        on = "aligns on" in shape
        for i in range(N_READS):
            s, end = A - 60 + i, (A - 9 if on else A + 1)
            if on:
                r1 = list(contig[s : s + 100])
                r1[A + 1 - s] = base
                seq, cig = "".join(r1), ((0, 100),)
                r2, r2_cig = contig[end - 80 : end + 20], ((0, 80), (4, 20))
            else:
                tail = ADAPTER[1 : 100 - (A + 2 - s) + 1]
                seq, cig = contig[s : A + 1] + base + tail, ((0, A + 2 - s), (4, len(tail)))
                r2, r2_cig = contig[end - 100 : end], ((0, 100),)
            a = make_read(f"f{i}", seq, s, cig, flag=99)
            b = make_read(f"f{i}", r2, end - len(r2) + (20 if on else 0), r2_cig, flag=147)
            a.template_length, b.template_length = end - s, s - end
            a.next_reference_id = b.next_reference_id = 0
            a.next_reference_start, b.next_reference_start = b.reference_start, s
            reads += [a, b]
        return contig, reads, _rs.Variant("1", A + 1, contig[A + 1], "A", "SNP")
    if shape.startswith("absent"):
        for i in range(N_READS):
            s = A - 50 + i
            seq = list(contig[s : s + 100])
            seq[A + 1 - s] = "A"
            a = make_read(f"q{i}", "".join(seq), s, ((0, 100),))
            a.query_qualities = None  # QUAL '*': the BAM stores 0xFF
            reads.append(a)
        return contig, reads, _rs.Variant("1", A + 1, contig[A + 1], "A", "SNP")
    clip = "hard clip" in shape
    for i in range(N_READS):
        s = A - 60 + i
        cig = ((5, 20), (0, 70)) if clip else ((0, 70),)
        reads.append(make_read(f"h{i}", contig[s : s + 70], s, cig))
    # Unprepared (no reference context): the exact-carrier rule cannot judge it.
    return contig, reads, _rs.Variant("1", A, contig[A : A + 3], "A", "COMPLEX")


def _pure_group(shape: str) -> str:
    if shape.startswith("REF read"):
        return "C10 REF reads"
    if shape == "exact carrier" or shape.startswith("the ALT at another placement"):
        return "invariant 7 carriers"
    if "inside the window" in shape:
        return "C26 inside the window"
    if "outside the window" in shape:
        return "C26 outside the window"
    if shape.startswith("the ALT written as"):
        return "C27 the ALT across ops"
    return "C28 anchor deleted"


def run(case: Case, workdir: Path) -> tuple[tuple[int, int, int], str]:
    """Count the case's reads; returns ((rd, ad, partial), census verdicts)."""
    d = Path(tempfile.mkdtemp(dir=workdir))
    if case.variant in PURE:
        motif, ref, alt = PURE[case.variant]
        contig = _contig(motif)
        siblings: list = []
        mask = None
        if case.group == "C2 siblings":
            snv_at = A + 3 if "inside" in case.shape else tract(A, ref, alt, contig)[1] + 6
            snv_alt = "C" if contig[snv_at] != "C" else "G"
            events, s0, length = [(snv_at, "X", snv_alt)], 0, 100
        elif case.group == "C35 unreadable inserts":
            (events, mask), s0, length = _unreadable_shapes(contig, ref, alt)[case.shape], 0, 100
        elif case.group == "C36 inserts of other bases":
            events, s0, length = _other_bases_shapes(contig, ref, alt)[case.shape], 0, 100
        elif case.group == "C37 inserts behind a substitution":
            events, s0, length = _behind_substitution_shapes(contig, ref, alt)[case.shape], 0, 100
        else:
            events, s0, length = _pure_shapes(contig, ref, alt)[case.shape]
        reads = []
        for i in range(N_READS):
            s = A - 50 + i + s0
            seq, cig = _build(contig, s, length - i, events)
            reads.append(make_read(f"r{i}", seq, s, cig))
            if mask:
                _mask_inserted(reads[-1], mask)
        fa, bam = write_contig(d, contig, reads, "c")
        rows = [_rs.Variant("1", A, ref, alt, "X")]
        if case.group == "C2 siblings":
            rows.append(_rs.Variant("1", snv_at, contig[snv_at], snv_alt, "SNP"))
        pvs = _rs.prepare_variants(rows, fa, 5, False, 1, True)
        vs = [p.variant for p in pvs]
        if case.group == "C2 siblings":
            siblings = [[vs[1]], [vs[0]]]
        counts = _rs.count_bam_binned(
            bam, vs, [None] * len(vs), sibling_variants=siblings or [[] for _ in vs], **ARGS
        )
        c = counts[0]
        verdicts = Counter(x.value for x in census(bam, contig, vs[0]).counts.elements())
        return (c.rd, c.ad, c.partial_alt), ", ".join(
            f"{k} {n}" for k, n in sorted(verdicts.items())
        )
    if case.variant == "read inputs":
        contig, reads, v = _read_input_bam(d, case.shape)
        _, bam = write_contig(d, contig, reads, "c")
        (c,) = _rs.count_bam_binned(bam, [v], [None], **ARGS)
        return (c.rd, c.ad, c.partial_alt), ""
    if case.variant in DELINS:
        n, alt = DELINS[case.variant]
        contig = _contig("")
        reads = []
        for i, (s, length, events) in enumerate(_delins_reads(n, alt, case.shape)):
            seq, cig = _build(contig, s, length, events)
            reads.append(make_read(f"r{i}", seq, s, cig))
        fa, bam = write_contig(d, contig, reads, "c")
        (pv,) = _rs.prepare_variants(
            [_rs.Variant("1", A, contig[A : A + n], alt, "X")], fa, 5, False, 1, True
        )
        (c,) = _rs.count_bam_binned(bam, [pv.variant], [None], **ARGS)
        return (c.rd, c.ad, c.partial_alt), ""
    run_len, r, a = COMPLEX[case.variant]
    contig = _complex_contig(run_len)
    haps, ends = _complex_reads(run_len, r, a)
    h, e = case.shape.rsplit(", ", 1)
    events, last = haps[h], ends[e]
    reads = []
    for i in range(N_READS):
        s = A - 40 - 3 * i  # same last base, different starts
        gap = sum((len(x) if op == "I" else -x) for _, op, x in events if op in "ID")
        seq, cig = _build(contig, s, last + 1 - s + gap, events)
        reads.append(make_read(f"r{i}", seq, s, cig))
    fa, bam = write_contig(d, contig, reads, "c")
    (pv,) = _rs.prepare_variants([_rs.Variant("1", A, r, a, "X")], fa, 5, False, 1, True)
    (c,) = _rs.count_bam_binned(bam, [pv.variant], [None], **ARGS)
    return (c.rd, c.ad, c.partial_alt), ""


# The decision each group follows: decided (the call is the rule) or open (the
# call is today's, pinned until the decision lands).
DECISIONS = {
    "C10 REF reads": "decided: REF only from informative reads (C10 #157)",
    "invariant 7 carriers": "decided: ALT only where the read's bases hold it (#161, #188, #191, #192, #204)",
    "C1 exact carriers": "decided: complex variants count exact carriers (C1 #141)",
    "C2 siblings": "decided: a sibling's event inside the window withdraws REF (C2 #119)",
    "C26 inside the window": "decided: another indel inside the window, not REF (C26 #200)",
    "C26 outside the window": "decided: another indel outside the window is a separate event (C26 #200)",
    "C27 the ALT across ops": "decided: the ALT across several ops counts ALT, by its bases (C27 #201)",
    "C28 anchor deleted": "decided: a read deleting the anchor is judged by its bases (C28 #202)",
    "C25 long events": "decided: junction windows read on through the read (C25 #199)",
    "C39 one-sided windows": "decided: one pair of equal-length junction windows at one flank (C39 #253)",
    "read inputs": "decided: a read contributes its molecule's bases, with qualities (C17 #176, C19 #182; C29 #207 a fix)",
    "C35 unreadable inserts": "decided: an insertion's ALT needs one of the read's own inserted bases read (C35 #240)",
    "C36 inserts of other bases": "decided: a same-length insert of other bases is judged by its bases, never Phase 3 (C36 #243)",
    "C37 inserts behind a substitution": "decided: the junction ALT stays; a carrier behind a substitution is another allele, not REF (C37 #245)",
}

# Decided calls the engine does not give yet: the red-first cases of a rule being
# built, strict xfails in test_read_judgment_spec.py until it lands.
PENDING: set[str] = set()

# (ref_count, alt_count, partial_alt) for the case's four reads.
EXPECT = {
    "hp-A | REF read spanning the tract": (4, 0, 0),
    "hp-A | REF read ending inside the tract": (0, 0, 0),
    "hp-A | exact carrier": (0, 4, 0),
    "hp-A | carrier with another indel outside the window": (0, 4, 0),
    "hp-A | the ALT at another placement: D1 near the anchor": (0, 4, 0),
    "hp-A | another indel inside the window: I1 near the anchor": (0, 0, 4),
    "hp-A | another indel inside the window: D2 after the anchor": (0, 0, 4),
    "hp-A | another indel outside the window: D1 past the tract": (4, 0, 0),
    "hp-A | another indel outside the window: D5 past the tract": (4, 0, 0),
    "hp-A | another indel outside the window: I1 past the tract": (4, 0, 0),
    "hp-AA | REF read spanning the tract": (4, 0, 0),
    "hp-AA | REF read ending inside the tract": (0, 0, 0),
    "hp-AA | exact carrier": (0, 4, 0),
    "hp-AA | carrier with another indel outside the window": (0, 4, 0),
    "hp-AA | another indel inside the window: D1 near the anchor": (0, 0, 4),
    "hp-AA | another indel inside the window: I1 near the anchor": (0, 0, 4),
    "hp-AA | the ALT at another placement: D2 after the anchor": (0, 4, 0),
    "hp-AA | another indel outside the window: D1 past the tract": (4, 0, 0),
    "hp-AA | another indel outside the window: D5 past the tract": (4, 0, 0),
    "hp-AA | another indel outside the window: I1 past the tract": (4, 0, 0),
    "hp-AA | the ALT written as two deletions": (0, 4, 0),
    "hp+A | REF read spanning the tract": (4, 0, 0),
    "hp+A | REF read ending inside the tract": (0, 0, 0),
    "hp+A | exact carrier": (0, 4, 0),
    "hp+A | carrier with another indel outside the window": (0, 4, 0),
    "hp+A | another indel inside the window: D1 near the anchor": (0, 0, 4),
    "hp+A | the ALT at another placement: I1 near the anchor": (0, 4, 0),
    "hp+A | another indel inside the window: D2 after the anchor": (0, 0, 4),
    "hp+A | another indel outside the window: D1 past the tract": (4, 0, 0),
    "hp+A | another indel outside the window: D5 past the tract": (4, 0, 0),
    "hp+A | another indel outside the window: I1 past the tract": (4, 0, 0),
    "u-8 | REF read spanning the tract": (4, 0, 0),
    "u-8 | REF read ending inside the tract": (0, 0, 0),
    "u-8 | exact carrier": (0, 4, 0),
    "u-8 | carrier with another indel outside the window": (0, 4, 0),
    "u-8 | another indel inside the window: D1 near the anchor": (0, 0, 4),
    "u-8 | another indel inside the window: I1 near the anchor": (0, 0, 4),
    "u-8 | another indel inside the window: D2 after the anchor": (0, 0, 4),
    "u-8 | another indel outside the window: D1 past the tract": (4, 0, 0),
    "u-8 | another indel outside the window: D5 past the tract": (4, 0, 0),
    "u-8 | another indel outside the window: I1 past the tract": (4, 0, 0),
    "u-8 | the ALT written as two deletions": (0, 4, 0),
    "u-8 | anchor deleted, then an insertion": (0, 0, 4),
    "u-8 | anchor deleted": (0, 0, 4),
    "u+10 | anchor deleted": (0, 0, 4),
    "dup+8 | anchor deleted": (0, 0, 4),
    "hp-A | the ALT written as D2 + I1 at the anchor": (0, 4, 0),
    "hp-A | the ALT written as D2 + I1 inside the run": (0, 4, 0),
    "hp+A | the ALT written as I2 + D1 at the anchor": (0, 4, 0),
    "u+10 | REF read spanning the tract": (4, 0, 0),
    "u+10 | REF read ending inside the tract": (4, 0, 0),
    "u+10 | exact carrier": (0, 4, 0),
    "u+10 | carrier with another indel outside the window": (0, 4, 0),
    "u+10 | another indel outside the window: D1 near the anchor": (4, 0, 0),
    "u+10 | another indel outside the window: I1 near the anchor": (4, 0, 0),
    "u+10 | another indel outside the window: D2 after the anchor": (4, 0, 0),
    "u+10 | another indel outside the window: D1 past the tract": (4, 0, 0),
    "u+10 | another indel outside the window: D5 past the tract": (4, 0, 0),
    "u+10 | another indel outside the window: I1 past the tract": (4, 0, 0),
    "ca-CA | REF read spanning the tract": (4, 0, 0),
    "ca-CA | REF read ending inside the tract": (0, 0, 0),
    "ca-CA | exact carrier": (0, 4, 0),
    "ca-CA | carrier with another indel outside the window": (0, 4, 0),
    "ca-CA | another indel inside the window: D1 near the anchor": (0, 0, 4),
    "ca-CA | another indel inside the window: I1 near the anchor": (0, 0, 4),
    "ca-CA | the ALT at another placement: D2 after the anchor": (0, 4, 0),
    "ca-CA | another indel outside the window: D1 past the tract": (4, 0, 0),
    "ca-CA | another indel outside the window: D5 past the tract": (4, 0, 0),
    "ca-CA | another indel outside the window: I1 past the tract": (4, 0, 0),
    "ca-CA | the ALT written as two deletions": (0, 4, 0),
    "dup-8 | REF read spanning the tract": (4, 0, 0),
    "dup-8 | REF read ending inside the tract": (0, 0, 0),
    "dup-8 | exact carrier": (0, 4, 0),
    "dup-8 | carrier with another indel outside the window": (0, 4, 0),
    "dup-8 | another indel inside the window: D1 near the anchor": (0, 0, 4),
    "dup-8 | another indel inside the window: I1 near the anchor": (0, 0, 4),
    "dup-8 | another indel inside the window: D2 after the anchor": (0, 0, 4),
    "dup-8 | another indel outside the window: D1 past the tract": (4, 0, 0),
    "dup-8 | another indel outside the window: D5 past the tract": (4, 0, 0),
    "dup-8 | another indel outside the window: I1 past the tract": (4, 0, 0),
    "dup-8 | the ALT written as two deletions": (0, 4, 0),
    "dup-8 | anchor deleted, then an insertion": (0, 0, 4),
    "dup-8 | anchor deleted": (0, 0, 4),
    "dup+8 | REF read spanning the tract": (4, 0, 0),
    "dup+8 | REF read ending inside the tract": (0, 0, 0),
    "dup+8 | exact carrier": (0, 4, 0),
    "dup+8 | carrier with another indel outside the window": (0, 4, 0),
    "dup+8 | another indel inside the window: D1 near the anchor": (0, 0, 4),
    "dup+8 | another indel inside the window: I1 near the anchor": (0, 0, 4),
    "dup+8 | another indel inside the window: D2 after the anchor": (0, 0, 4),
    "dup+8 | another indel outside the window: D1 past the tract": (4, 0, 0),
    "dup+8 | another indel outside the window: D5 past the tract": (4, 0, 0),
    "dup+8 | another indel outside the window: I1 past the tract": (4, 0, 0),
    "C>TA run10 | REF, ends inside the run": (0, 0, 0),
    "C>TA run10 | REF, ends on the base after the run": (0, 0, 0),
    "C>TA run10 | REF, ends past the run": (4, 0, 0),
    "C>TA run10 | substitution only, ends inside the run": (0, 0, 0),
    "C>TA run10 | substitution only, ends on the base after the run": (0, 0, 0),
    "C>TA run10 | substitution only, ends past the run": (0, 0, 0),
    "C>TA run10 | anchor kept, one A more, ends inside the run": (0, 0, 0),
    "C>TA run10 | anchor kept, one A more, ends on the base after the run": (0, 0, 0),
    # one base from each allele (an A more than REF, C for ALT's T): neither (RJ-23)
    "C>TA run10 | anchor kept, one A more, ends past the run": (0, 0, 0),
    "C>TA run10 | anchor kept, one A fewer, ends inside the run": (0, 0, 0),
    "C>TA run10 | anchor kept, one A fewer, ends on the base after the run": (0, 0, 0),
    "C>TA run10 | anchor kept, one A fewer, ends past the run": (0, 0, 0),
    "C>TA run10 | exact carrier, ends inside the run": (0, 0, 0),
    # holds the REF window's 15 bases: the T, all eleven A's and the G after (RJ-23)
    "C>TA run10 | exact carrier, ends on the base after the run": (0, 4, 0),
    "C>TA run10 | exact carrier, ends past the run": (0, 4, 0),
    "C>TA run60 | REF, ends inside the run": (4, 0, 0),
    "C>TA run60 | REF, ends on the base after the run": (4, 0, 0),
    "C>TA run60 | REF, ends past the run": (4, 0, 0),
    "C>TA run60 | substitution only, ends inside the run": (0, 4, 0),
    "C>TA run60 | substitution only, ends on the base after the run": (0, 0, 0),
    "C>TA run60 | substitution only, ends past the run": (0, 0, 0),
    "C>TA run60 | anchor kept, one A more, ends inside the run": (4, 0, 0),
    "C>TA run60 | anchor kept, one A more, ends on the base after the run": (0, 0, 0),
    "C>TA run60 | anchor kept, one A more, ends past the run": (0, 0, 0),
    "C>TA run60 | anchor kept, one A fewer, ends inside the run": (4, 0, 0),
    "C>TA run60 | anchor kept, one A fewer, ends on the base after the run": (0, 0, 0),
    "C>TA run60 | anchor kept, one A fewer, ends past the run": (0, 0, 0),
    "C>TA run60 | exact carrier, ends inside the run": (0, 4, 0),
    "C>TA run60 | exact carrier, ends on the base after the run": (0, 4, 0),
    "C>TA run60 | exact carrier, ends past the run": (0, 4, 0),
    "CA>T run60 | REF, ends inside the run": (4, 0, 0),
    "CA>T run60 | REF, ends on the base after the run": (4, 0, 0),
    "CA>T run60 | REF, ends past the run": (4, 0, 0),
    "CA>T run60 | substitution only, ends inside the run": (0, 4, 0),
    "CA>T run60 | substitution only, ends on the base after the run": (0, 0, 0),
    "CA>T run60 | substitution only, ends past the run": (0, 0, 0),
    "CA>T run60 | anchor kept, one A more, ends inside the run": (4, 0, 0),
    "CA>T run60 | anchor kept, one A more, ends on the base after the run": (0, 0, 0),
    "CA>T run60 | anchor kept, one A more, ends past the run": (0, 0, 0),
    "CA>T run60 | anchor kept, one A fewer, ends inside the run": (4, 0, 0),
    "CA>T run60 | anchor kept, one A fewer, ends on the base after the run": (0, 0, 0),
    "CA>T run60 | anchor kept, one A fewer, ends past the run": (0, 0, 4),
    "CA>T run60 | exact carrier, ends inside the run": (0, 4, 0),
    "CA>T run60 | exact carrier, ends on the base after the run": (0, 4, 0),
    "CA>T run60 | exact carrier, ends past the run": (0, 4, 0),
    # C39 #253: one-sided, equal-length junction windows (RJ-23).
    "delins 28>6 | exact carrier, eight flank bases right": (0, 4, 0),
    "delins 28>6 | exact carrier, four flank bases left": (0, 4, 0),
    "delins 28>6 | exact carrier, long flanks": (0, 4, 0),
    "delins 28>6 | REF read, four flank bases left": (4, 0, 0),
    "delins 28>6 | REF read, four flank bases right": (4, 0, 0),
    "delins 28>6 | REF read starting inside the event": (0, 0, 0),
    # the change is a separate event, masked: the carriers past it count ALT (RJ-23)
    "delins 28>6 | carriers ending before a change that recurs on ALT reads only": (0, 4, 0),
    "delins 60>2 | exact carrier, eight flank bases right": (0, 4, 0),
    "delins 60>2 | exact carrier, long flanks": (0, 4, 0),
    "delins 60>2 | REF read, four flank bases left": (4, 0, 0),
    "delins 60>2 | REF read starting inside the event": (0, 0, 0),
    "u-8 | sibling SNV inside the span": (0, 0, 4),
    "u-8 | sibling SNV outside the window": (4, 0, 0),
    "read inputs | absent base qualities (QUAL '*'), showing the ALT": (0, 0, 0),
    "read inputs | previous complex classifier: REF read with a leading hard clip": (4, 0, 0),
    "read inputs | previous complex classifier: REF read, no clip": (4, 0, 0),
    "read inputs | read-through: the base on the SNV is adapter, showing REF": (0, 0, 0),
    "read inputs | read-through: the base on the SNV is adapter, showing the ALT": (0, 0, 0),
    "read inputs | read-through: the read aligns on past its TLEN end (R2's 5' end clipped), showing the ALT": (
        0,
        4,
        0,
    ),
    # C35 #240: the read's own inserted bases unreadable (N at Q2, or letters at Q5).
    "hp+A | an unreadable insert at the anchor (N)": (0, 0, 4),
    "hp+A | an unreadable insert inside the run (N)": (0, 0, 4),
    "hp+A | an insert below min BQ inside the run": (0, 0, 4),
    "u+10 | an unreadable insert (N)": (0, 0, 4),
    "u+10 | an insert below min BQ": (0, 0, 4),
    "u+10 | a partly masked insert": (0, 4, 0),
    "u+10 | anchor deleted, the insert unreadable": (0, 0, 4),
    # C36 #243: same-length inserts of other bases, judged by their bases.
    "u+10 | a same-length insert of other bases one junction left": (4, 0, 0),
    "u+10 | a same-length insert of other bases four junctions right": (4, 0, 0),
    "u+10 | the ALT written one junction off with a compensating mismatch": (0, 4, 0),
    "u+10 | the ALT with an error at its first inserted base, written one junction right": (
        0,
        0,
        4,
    ),
    "dup+8 | a same-length insert of other bases inside the tract": (0, 0, 4),
    # C37 #245: the strict path keeps its junction ALT; a readable insert whose
    # fewest-mismatch placement sits inside the window is another allele (RJ-22).
    "hp+A | the ALT at the junction with a substitution in the run": (0, 4, 0),
    "hp+A | the ALT with its anchor substituted, written one junction left": (0, 0, 4),
    "u+10 | the ALT with its anchor substituted, written one junction left": (0, 0, 4),
}


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as tmp:
        for case in cases():
            got, verdicts = run(case, Path(tmp))
            print(f"{case.group}\t{case.id}\t{got}\t{verdicts}", flush=True)
