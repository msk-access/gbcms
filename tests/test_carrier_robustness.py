"""Reads judged by their bases never crash the count (#278).

Pure indels and delins are judged by the exact-carrier rule, which places a read's own
indels by its bases over the reference prepared around the variant. Real reads run far
past that reference, carry indels well beside the event, soft clips, N and low-quality
bases: a seeded sweep of such reads at random pure indels and delins must count without
error and keep the counting invariants.
"""

import random

from helpers import count_checked, make_read, write_contig

from gbcms import _rs as gbcms_rs

READ = 250


def _contig(rng: random.Random, n: int = 3000) -> str:
    motifs = ["A", "CA", "GGC", "ACGT"]
    seq = []
    while len(seq) < n:
        if rng.random() < 0.15:
            seq += list(rng.choice(motifs) * rng.randint(2, 12))
        else:
            seq.append(rng.choice("ACGT"))
    return "".join(seq[:n])


def _variant(rng: random.Random, contig: str):
    """A pure insertion, a pure deletion, or a delins."""
    p = rng.randint(800, 2200)
    kind = rng.randrange(3)
    if kind == 0:
        ins = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 12)))
        return p, contig[p], contig[p] + ins
    if kind == 1:
        k = rng.randint(1, 20)
        return p, contig[p : p + 1 + k], contig[p]
    ref = contig[p : p + rng.randint(2, 6)]
    alt = ref
    while alt == ref or alt[0] == ref[0] and alt[-1] == ref[-1]:
        alt = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 5)))
    return p, ref, alt


def _read(rng: random.Random, contig: str, p: int, ref: str, alt: str, name: str):
    """A read over the variant: the ALT or REF, then up to two more indels anywhere,
    substitutions, N, low-quality bases and a soft clip at either end."""
    hap, events = list(contig), []
    if rng.random() < 0.5:
        hap[p : p + len(ref)] = list(alt)
        events.append((p + 1, len(alt) - len(ref)))
    shift = sum(n for _, n in events)
    s = rng.randint(p - READ + 10, p - 5)
    for _ in range(rng.randint(0, 2)):
        at = rng.randint(max(s + 5, p - 70), p + 70)
        n = rng.choice([-3, -1, 1, 2, 5])
        if any(abs(at - e) < 4 for e, _ in events):
            continue
        h = at + (shift if at > p else 0)
        if n > 0:
            hap[h:h] = list("".join(rng.choice("ACGT") for _ in range(n)))
        else:
            del hap[h : h - n]
        events.append((at, n))
    seq = hap[s : s + READ]
    quals = [30] * READ
    for _ in range(rng.randint(0, 6)):
        i = rng.randrange(READ)
        r = rng.random()
        if r < 0.4:
            seq[i] = rng.choice("ACGT")
        elif r < 0.7:
            seq[i] = "N"
            quals[i] = 2
        else:
            quals[i] = 5
    # the CIGAR from the read's own events, in genome order
    ops, gpos, qpos = [], s, 0
    for at, n in sorted(events):
        if at <= gpos or at >= s + READ - 20:
            continue
        ops.append((0, at - gpos))
        qpos += at - gpos
        if n > 0:
            ops.append((1, n))
            qpos += n
            gpos = at
        else:
            ops.append((2, -n))
            gpos = at - n
    ops.append((0, READ - qpos))
    clip_l, clip_r = rng.choice([0, 0, 0, 7]), rng.choice([0, 0, 0, 9])
    if clip_l and ops[0][0] == 0 and ops[0][1] > clip_l:
        ops[0] = (0, ops[0][1] - clip_l)
        ops.insert(0, (4, clip_l))
        s += clip_l
    if clip_r and ops[-1][0] == 0 and ops[-1][1] > clip_r:
        ops[-1] = (0, ops[-1][1] - clip_r)
        ops.append((4, clip_r))
    return make_read(name, "".join(seq), s, tuple(ops), quals=quals)


def test_random_reads_at_random_indels_count_without_error(tmp_path):
    rng = random.Random(278)
    for case in range(40):
        contig = _contig(rng)
        p, ref, alt = _variant(rng, contig)
        reads = [_read(rng, contig, p, ref, alt, f"q{i}") for i in range(60)]
        d = tmp_path / f"c{case}"
        d.mkdir()
        fa, bam = write_contig(d, contig, reads, "c")
        (pv,) = gbcms_rs.prepare_variants(
            [gbcms_rs.Variant("1", p, ref, alt, "X")], fa, 5, False, 1, True
        )
        (c,) = count_checked(bam, [pv.variant])
        assert c.dp >= c.rd + c.ad
        assert c.dpf >= c.rdf + c.adf
        assert c.rd == c.rd_fwd + c.rd_rev
        assert c.ad == c.ad_fwd + c.ad_rev
