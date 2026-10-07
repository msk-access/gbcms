"""Compare two builds' regression-panel outputs, and each against the sign-out counts
(docs/development/regression-panel.md; the reports name runs and loci, so they stay
where the BAMs are).

usage: compare_panel.py TIER_DIR OUT_ROOT BASE_BUILD NEW_BUILD [REPORT_DIR]

Reads TIER_DIR/runs.tsv, the builds' outputs and their times.<shard>.tsv (a run's
latest attempt is its outcome; failed runs are counted and left out of the
comparison). The panel MAFs' signout_* and panel_strata columns pass through gbcms.
Rows are matched by locus (and occurrence, for a locus listed twice), never by
position. Writes to REPORT_DIR (default OUT_ROOT/compare_BASE_vs_NEW):
  gate.txt                    the summary the release gate reads
  header_diff.tsv             columns (MAF) or INFO/FORMAT IDs (VCF) in one version only
  version_summary.tsv         per arm: runs, rows, identical rows, count-changed rows,
                              rows (records) only in one version
  version_cells.tsv           changed cells per arm and column
  version_rows.tsv            every changed cell (run, arm, locus, column, base, new);
                              VCF output field by field (FILTER, INFO:<id>, FORMAT:<id>)
  concordance_by_stratum.tsv  BASE and NEW against the sign-out ALT count, per stratum
                              and level (read, fragment, and matched: the level the
                              sign-out counted, reads for IMPACT and duplex + simplex
                              fragments for ACCESS): n, median delta, fraction within
                              max(2, 10%), and its change. The mq0 arm's rows (PMS2 at
                              --min-mapq 0) are strata prefixed "mq0:".
  discordant.tsv              the largest NEW vs sign-out disagreements at the matched
                              level, and the top 20 per stratum: to adjudicate read by
                              read (the BAM is the truth; the sign-out is a comparison)
  normals.tsv                 matched-normal fillouts against the sign-out normal counts
  times.tsv                   wall time and peak memory per arm, BASE vs NEW
  parquet.tsv                 the mfsd arm's Parquet tables, BASE vs NEW (needs pyarrow)
Exits 1 when the gate needs attention: failed or missing runs, a stratum whose matched
concordance fell by more than 5 points, or an arm 1.5x slower or larger.
"""

import collections
import csv
import glob
import os
import statistics
import sys

COUNTS = (
    "ref_count",
    "alt_count",
    "partial_alt",
    "total_count",
    "ref_count_fragment",
    "alt_count_fragment",
)
KEY = ("Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2")
FALL = 0.05  # a stratum whose matched within-fraction falls more than this is adjudicated
SLOWER = 1.5  # an arm whose median wall time or peak memory grows past this is explained


def read_runs(tier):
    with open(os.path.join(tier, "runs.tsv")) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def output_file(out, build, rid):
    for ext in ("maf", "vcf"):
        f = glob.glob(os.path.join(out, build, rid, f"*.{ext}"))
        if f:
            return f[0], ext
    return None, None


def maf_rows(path):
    lines = [x for x in open(path) if not x.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def by_locus(rows):
    """Rows by locus, with the occurrence number for a locus listed twice."""
    seen, out = collections.Counter(), {}
    for y in rows:
        k = tuple(y.get(c, "") for c in KEY)
        out[k + (seen[k],)] = y
        seen[k] += 1
    return out


def vcf_records(path):
    """(the INFO and FORMAT IDs the header defines, {record key: {field: value}}): a
    record keyed by CHROM, POS, REF, ALT (and its occurrence), its fields FILTER,
    INFO:<id> and FORMAT:<id> (the first sample)."""
    ids, recs, seen = set(), {}, collections.Counter()
    for x in open(path):
        if x.startswith("##INFO=<ID=") or x.startswith("##FORMAT=<ID="):
            ids.add(x[2:].split("=", 1)[0] + ":" + x.split("ID=", 1)[1].split(",", 1)[0])
            continue
        if x.startswith("#"):
            continue
        f = x.rstrip("\n").split("\t")
        k = (f[0], f[1], f[3], f[4])
        fields = {"FILTER": f[6]}
        for kv in f[7].split(";"):
            key, _, val = kv.partition("=")
            fields["INFO:" + key] = val
        if len(f) > 9:
            for key, val in zip(f[8].split(":"), f[9].split(":"), strict=False):
                fields["FORMAT:" + key] = val
        recs[k + (seen[k],)] = fields
        seen[k] += 1
    return ids, recs


def as_int(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


def as_num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def within(d, so):
    return abs(d) <= max(2, 0.1 * so)


def attempts(out, build):
    """Each run's latest attempt (by end time, else file order): {run_id: row}."""
    latest = {}
    for p in sorted(glob.glob(os.path.join(out, build, "times.*.tsv"))):
        with open(p) as fh:
            for n, x in enumerate(csv.DictReader(fh, delimiter="\t")):
                order = (as_num(x.get("end_epoch")) or 0, p, n)
                if x["run_id"] not in latest or order > latest[x["run_id"]][0]:
                    latest[x["run_id"]] = (order, x)
    return {rid: x for rid, (_o, x) in latest.items()}


def main():
    tier, out, base, new = sys.argv[1:5]
    rep = sys.argv[5] if len(sys.argv) > 5 else os.path.join(out, f"compare_{base}_vs_{new}")
    os.makedirs(rep, exist_ok=True)
    runs = read_runs(tier)
    att = {b: attempts(out, b) for b in (base, new)}
    failed = {
        b: sorted(r for r, x in att[b].items() if str(x.get("exit")) != "0") for b in (base, new)
    }
    hdr = collections.defaultdict(lambda: [set(), set()])
    summ = collections.defaultdict(collections.Counter)
    cells = collections.Counter()
    changed = []
    conc_in = {b: {} for b in (base, new)}  # (tag, locus key) -> sums over the tag's runs
    flavours = collections.defaultdict(set)  # ACCESS tag -> flavours compared
    normals = []
    missing = collections.Counter()
    for r in runs:
        rid, arm = r["run_id"], r["arm"]
        if rid in failed[base] or rid in failed[new]:
            continue
        fb, eb = output_file(out, base, rid)
        fn, en = output_file(out, new, rid)
        if fn is None or fb is None:
            missing[(arm, "base" if fb is None else "new")] += 1
            continue
        summ[arm]["runs"] += 1
        if en == "vcf":
            ib, vb = vcf_records(fb)
            in_, vn = vcf_records(fn)
            hdr[(arm, "vcf")][0].update(ib - in_)
            hdr[(arm, "vcf")][1].update(in_ - ib)
            one_side = ib ^ in_  # declared in one version only: a header difference
            summ[arm]["records only in base"] += len(set(vb) - set(vn))
            for k, y in vn.items():
                x = vb.get(k)
                summ[arm]["records"] += 1
                if x is None:
                    summ[arm]["records only in new"] += 1
                    continue
                fields = (set(x) | set(y)) - one_side
                diff = sorted(c for c in fields if x.get(c) != y.get(c))
                if not diff:
                    summ[arm]["records identical"] += 1
                for c in diff:
                    cells[(arm, c)] += 1
                    changed.append(
                        (rid, arm, k[0], k[1], "", k[2], k[3], c, x.get(c, ""), y.get(c, ""))
                    )
            continue
        lb, ln = by_locus(maf_rows(fb)), by_locus(maf_rows(fn))
        if lb and ln:
            any_b, any_n = next(iter(lb.values())), next(iter(ln.values()))
            hdr[(arm, "maf")][0].update(set(any_b) - set(any_n))
            hdr[(arm, "maf")][1].update(set(any_n) - set(any_b))
        if len(lb) != len(ln):
            summ[arm]["runs with a different row count"] += 1
        summ[arm]["rows only in base"] += len(set(lb) - set(ln))
        summ[arm]["rows only in new"] += len(set(ln) - set(lb))
        for k, y in ln.items():
            x = lb.get(k)
            if x is None:
                continue
            summ[arm]["rows"] += 1
            diff = [c for c in y if c in x and x[c] != y[c]]
            if not diff:
                summ[arm]["rows identical"] += 1
            if any(x.get(c) != y.get(c) for c in COUNTS):
                summ[arm]["rows count-changed"] += 1
            for c in diff:
                cells[(arm, c)] += 1
                changed.append((rid, arm, *k[:5], c, x[c], y[c]))
        if arm == "normal":
            for k, y in ln.items():
                normals.append(
                    (
                        rid,
                        *k[:5],
                        y.get("signout_n_alt_count", ""),
                        y.get("alt_count", ""),
                        y.get("alt_count_fragment", ""),
                        y.get("total_count", ""),
                    )
                )
            continue
        if arm not in ("dna", "mq0"):
            continue
        access = rid.endswith(("_duplex", "_simplex"))
        if access:
            flavours[(arm, r["tag"])].add(rid.rsplit("_", 1)[1])
        for b, rows_ in ((base, lb), (new, ln)):
            for k, y in rows_.items():
                v = conc_in[b].setdefault(
                    (arm, r["tag"], k),
                    {"alt": 0, "altf": 0, "so": None, "strata": "", "access": access},
                )
                v["alt"] += as_int(y.get("alt_count")) or 0
                v["altf"] += as_int(y.get("alt_count_fragment")) or 0
                v["so"] = as_int(y.get("signout_t_alt_count"))
                v["strata"] = y.get("panel_strata", "")

    # Concordance by stratum; an ACCESS sample missing a flavour is left out.
    incomplete = {k for k, f in flavours.items() if f != {"duplex", "simplex"}}
    conc = {b: collections.defaultdict(list) for b in (base, new)}
    disc = []
    for b in (base, new):
        for (arm, tag, k), v in conc_in[b].items():
            if v["so"] is None or (v["access"] and (arm, tag) in incomplete):
                continue
            m = v["altf"] if v["access"] else v["alt"]
            prefix = "mq0:" if arm == "mq0" else ""
            for s in [x for x in v["strata"].split(";") if x] + ["all"]:
                s = prefix + s
                conc[b][(s, "read")].append((v["alt"] - v["so"], v["so"]))
                conc[b][(s, "fragment")].append((v["altf"] - v["so"], v["so"]))
                conc[b][(s, "matched")].append((m - v["so"], v["so"]))
            if b == new:
                disc.append(
                    (abs(m - v["so"]), tag, *k[:5], v["so"], v["alt"], v["altf"], v["strata"])
                )
    worst = []
    with open(os.path.join(rep, "concordance_by_stratum.tsv"), "w") as fh:
        fh.write(
            "stratum\tlevel\tn\tmedian_delta_base\tmedian_delta_new\twithin_base\twithin_new\t"
            "within_change\n"
        )
        for k in sorted(set(conc[base]) | set(conc[new])):
            db, dn = conc[base].get(k, []), conc[new].get(k, [])
            if not dn:
                continue
            wb = sum(within(d, so) for d, so in db) / len(db) if db else float("nan")
            wn = sum(within(d, so) for d, so in dn) / len(dn)
            mb = statistics.median(d for d, _ in db) if db else float("nan")
            mn = statistics.median(d for d, _ in dn)
            fh.write(
                f"{k[0]}\t{k[1]}\t{len(dn)}\t{mb:.1f}\t{mn:.1f}\t{wb:.3f}\t{wn:.3f}\t{wn - wb:+.3f}\n"
            )
            if k[1] == "matched" and db and wn - wb < 0:
                worst.append((wn - wb, k[0], len(dn)))
    with open(os.path.join(rep, "discordant.tsv"), "w") as fh:
        fh.write(
            "abs_delta_matched\ttag\t"
            + "\t".join(KEY)
            + "\tsignout_alt\talt_count\talt_count_fragment\tstrata\n"
        )
        top = sorted(disc, reverse=True)
        picked = set()
        for r in top[:500]:
            picked.add(r[1:7])
            fh.write("\t".join(map(str, r)) + "\n")
        per = collections.Counter()
        for r in top:
            for s in [x for x in r[-1].split(";") if x]:
                if per[s] < 20 and r[1:7] not in picked:
                    per[s] += 1
                    picked.add(r[1:7])
                    fh.write("\t".join(map(str, r)) + "\n")

    with open(os.path.join(rep, "header_diff.tsv"), "w") as fh:
        fh.write("arm\tformat\tonly_in_base\tonly_in_new\n")
        for (arm, fmt), (ob, on) in sorted(hdr.items()):
            fh.write(f"{arm}\t{fmt}\t{','.join(sorted(ob))}\t{','.join(sorted(on))}\n")
    with open(os.path.join(rep, "version_summary.tsv"), "w") as fh:
        keys = sorted({k for c in summ.values() for k in c})
        fh.write("arm\t" + "\t".join(keys) + "\n")
        for arm, c in sorted(summ.items()):
            fh.write(arm + "\t" + "\t".join(str(c.get(k, 0)) for k in keys) + "\n")
    with open(os.path.join(rep, "version_cells.tsv"), "w") as fh:
        fh.write("arm\tcolumn\tchanged_cells\n")
        for (arm, c), n in cells.most_common():
            fh.write(f"{arm}\t{c}\t{n}\n")
    with open(os.path.join(rep, "version_rows.tsv"), "w") as fh:
        fh.write("run_id\tarm\t" + "\t".join(KEY) + "\tcolumn\tbase\tnew\n")
        for r in changed:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(os.path.join(rep, "normals.tsv"), "w") as fh:
        fh.write(
            "run_id\t"
            + "\t".join(KEY)
            + "\tsignout_n_alt\talt_count\talt_count_fragment\ttotal_count\n"
        )
        for r in normals:
            fh.write("\t".join(map(str, r)) + "\n")

    # Times: each run's latest successful attempt.
    slow = []
    with open(os.path.join(rep, "times.tsv"), "w") as fh:
        fh.write(
            "arm\truns\tmedian_wall_base\tmedian_wall_new\twall_ratio\tmax_rss_base_kb\tmax_rss_new_kb\n"
        )
        per_arm = {b: collections.defaultdict(list) for b in (base, new)}
        for b in (base, new):
            for x in att[b].values():
                if str(x.get("exit")) == "0":
                    per_arm[b][x["arm"]].append(
                        (as_num(x.get("wall_s")), as_num(x.get("max_rss_kb")))
                    )
        for arm in sorted(set(per_arm[base]) | set(per_arm[new])):
            wb = [w for w, _ in per_arm[base].get(arm, []) if w is not None]
            wn = [w for w, _ in per_arm[new].get(arm, []) if w is not None]
            rb = [m for _, m in per_arm[base].get(arm, []) if m is not None]
            rn = [m for _, m in per_arm[new].get(arm, []) if m is not None]
            mwb = statistics.median(wb) if wb else None
            mwn = statistics.median(wn) if wn else None
            ratio = mwn / mwb if mwb and mwn is not None else None
            fh.write(
                f"{arm}\t{len(wn)}\t{mwb if mwb is not None else 'NA'}\t{mwn if mwn is not None else 'NA'}\t"
                f"{f'{ratio:.2f}' if ratio is not None else 'NA'}\t{max(rb) if rb else 'NA'}\t"
                f"{max(rn) if rn else 'NA'}\n"
            )
            if ratio is not None and ratio > SLOWER:
                slow.append(f"{arm} wall time x{ratio:.2f}")
            if rb and rn and max(rn) > SLOWER * max(rb):
                slow.append(f"{arm} peak memory x{max(rn) / max(rb):.2f}")

    # Parquet (mfsd arm), when pyarrow is available
    pq_lines = []
    try:
        import pyarrow.parquet as pq

        for r in runs:
            if r["arm"] != "mfsd":
                continue
            for fb in glob.glob(os.path.join(out, base, r["run_id"], "*.parquet")):
                fn = os.path.join(out, new, r["run_id"], os.path.basename(fb))
                if not os.path.exists(fn):
                    pq_lines.append((r["run_id"], os.path.basename(fb), "missing in new", "", ""))
                    continue
                a, b2 = pq.read_table(fb), pq.read_table(fn)
                same = a.schema.equals(b2.schema) and a.equals(b2)
                pq_lines.append(
                    (
                        r["run_id"],
                        os.path.basename(fb),
                        "identical" if same else "differs",
                        a.num_rows,
                        b2.num_rows,
                    )
                )
    except ImportError:
        pq_lines.append(("", "", "pyarrow not available: Parquet not compared", "", ""))
    with open(os.path.join(rep, "parquet.tsv"), "w") as fh:
        fh.write("run_id\ttable\tresult\trows_base\trows_new\n")
        for r in pq_lines:
            fh.write("\t".join(map(str, r)) + "\n")

    worst.sort()
    flagged = [(d, s, n) for d, s, n in worst if d < -FALL]
    attention = []
    with open(os.path.join(rep, "gate.txt"), "w") as fh:
        fh.write(f"runs in list: {len(runs)}; missing outputs: {dict(missing) or 0}\n")
        fh.write(f"failed runs: {base} {len(failed[base])}, {new} {len(failed[new])}\n")
        for b in (base, new):
            if failed[b]:
                fh.write(
                    f"  failed in {b}: {', '.join(failed[b][:20])}{' ...' if len(failed[b]) > 20 else ''}\n"
                )
        for arm, c in sorted(summ.items()):
            fh.write(f"{arm}: {dict(c)}\n")
        fh.write(
            f"changed cells: {sum(cells.values())} across {len({c for _, c in cells})} columns\n"
        )
        for (arm, fmt), (ob, on) in sorted(hdr.items()):
            if ob or on:
                fh.write(
                    f"columns ({arm}, {fmt}): only in {base}: {sorted(ob)}; only in {new}: {sorted(on)}\n"
                )
        if incomplete:
            fh.write(
                f"ACCESS samples missing a flavour (left out of concordance): {len(incomplete)}\n"
            )
        fh.write(
            "strata whose concordance with the sign-out fell (matched level: reads for IMPACT,"
            " fragments for ACCESS; within-fraction change):\n"
        )
        for d, s, n in worst[:10]:
            fh.write(f"  {s} ({n}): {d:+.3f}{'  ADJUDICATE' if d < -FALL else ''}\n")
        if not worst:
            fh.write("  none\n")
        for s in slow:
            fh.write(f"  EXPLAIN: {s}\n")
        if failed[base] or failed[new]:
            attention.append("failed runs")
        if missing:
            attention.append("missing outputs")
        if flagged:
            attention.append(f"{len(flagged)} strata fell more than {FALL:.0%}")
        if slow:
            attention.append("run time or memory")
        fh.write(f"GATE: {'ATTENTION: ' + '; '.join(attention) if attention else 'clean'}\n")
    print(open(os.path.join(rep, "gate.txt")).read())
    print(f"reports in {rep}")
    if attention:
        sys.exit(1)


if __name__ == "__main__":
    main()
