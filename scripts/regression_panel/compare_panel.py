"""Compare two builds' regression-panel outputs, and each against the sign-out counts
(docs/development/regression-panel.md; the reports name runs and loci, so they stay
where the BAMs are).

usage: compare_panel.py TIER_DIR OUT_ROOT BASE_BUILD NEW_BUILD [REPORT_DIR]

Reads TIER_DIR/runs.tsv and the panel MAFs (their signout_* and panel_strata columns
pass through gbcms unchanged). Writes to REPORT_DIR (default OUT_ROOT/compare_BASE_vs_NEW):
  gate.txt                    the summary the release gate reads
  header_diff.tsv             columns only in BASE or only in NEW, per arm and format
  version_summary.tsv         per arm: runs, rows, identical rows, count-changed rows
  version_cells.tsv           changed cells per arm and column
  version_rows.tsv            every changed cell (run, arm, locus, column, base, new)
  concordance_by_stratum.tsv  BASE and NEW against the sign-out ALT count, per stratum
                              and level (read, fragment): n, median delta, fraction
                              within max(2 reads, 10%), and the change in that fraction
  discordant.tsv              the largest NEW vs sign-out fragment disagreements, and
                              the top 20 per stratum: to adjudicate read by read (the
                              BAM is the truth; the sign-out is a comparison)
  normals.tsv                 matched-normal fillouts against the sign-out normal counts
  times.tsv                   wall time and peak memory per arm, BASE vs NEW
  parquet.tsv                 the mfsd arm's Parquet tables, BASE vs NEW (needs pyarrow)
ACCESS rows are compared with the sign-out after summing duplex and simplex.
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


def vcf_records(path):
    head, recs = [], []
    for x in open(path):
        if x.startswith("##"):
            continue
        (head if x.startswith("#") else recs).append(x.rstrip("\n"))
    return head, recs


def as_int(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


def within(d, so):
    return abs(d) <= max(2, 0.1 * so)


def main():
    tier, out, base, new = sys.argv[1:5]
    rep = sys.argv[5] if len(sys.argv) > 5 else os.path.join(out, f"compare_{base}_vs_{new}")
    os.makedirs(rep, exist_ok=True)
    runs = read_runs(tier)
    hdr = collections.defaultdict(lambda: [set(), set()])
    summ = collections.defaultdict(collections.Counter)
    cells = collections.Counter()
    changed = []
    # per (build, tag, key-index) sums for concordance: ACCESS duplex + simplex summed
    conc_in = {
        b: collections.defaultdict(lambda: {"alt": 0, "altf": 0, "so": None, "strata": ""})
        for b in (base, new)
    }
    normals = []
    missing = collections.Counter()
    for r in runs:
        rid, arm = r["run_id"], r["arm"]
        fb, eb = output_file(out, base, rid)
        fn, en = output_file(out, new, rid)
        if fn is None or fb is None:
            missing[(arm, "base" if fb is None else "new")] += 1
            continue
        summ[arm]["runs"] += 1
        if en == "vcf":
            hb, rb = vcf_records(fb)
            hn, rn = vcf_records(fn)
            summ[arm]["records"] += len(rn)
            summ[arm]["records changed"] += sum(
                1 for x, y in zip(rb, rn, strict=False) if x != y
            ) + abs(len(rb) - len(rn))
            if hb != hn:
                hdr[(arm, "vcf")][0].add("column header line differs")
            continue
        rb, rn = maf_rows(fb), maf_rows(fn)
        if rb and rn:
            hdr[(arm, "maf")][0].update(set(rb[0]) - set(rn[0]))
            hdr[(arm, "maf")][1].update(set(rn[0]) - set(rb[0]))
        if len(rb) != len(rn):
            summ[arm]["runs with a different row count"] += 1
        for x, y in zip(rb, rn, strict=False):  # a different row count is counted above
            key = tuple(y.get(k, "") for k in KEY)
            summ[arm]["rows"] += 1
            common = [c for c in y if c in x]
            diff = [c for c in common if x[c] != y[c]]
            if not diff:
                summ[arm]["rows identical"] += 1
            if any(x.get(c) != y.get(c) for c in COUNTS):
                summ[arm]["rows count-changed"] += 1
            for c in diff:
                cells[(arm, c)] += 1
                changed.append((rid, arm, *key, c, x[c], y[c]))
        if arm == "normal":
            for y in rn:
                normals.append(
                    (
                        rid,
                        *(y.get(k, "") for k in KEY),
                        y.get("signout_n_alt_count", ""),
                        y.get("alt_count", ""),
                        y.get("alt_count_fragment", ""),
                        y.get("total_count", ""),
                    )
                )
            continue
        if arm != "dna":
            continue
        for b, rows_ in ((base, rb), (new, rn)):
            for i, y in enumerate(rows_):
                v = conc_in[b][(r["tag"], i)]
                v["alt"] += as_int(y.get("alt_count")) or 0
                v["altf"] += as_int(y.get("alt_count_fragment")) or 0
                v["so"] = as_int(y.get("signout_t_alt_count"))
                v["strata"] = y.get("panel_strata", "")
                v["key"] = tuple(y.get(k, "") for k in KEY)

    # Concordance by stratum
    conc = {b: collections.defaultdict(list) for b in (base, new)}
    disc = []
    for b in (base, new):
        for (tag, _i), v in conc_in[b].items():
            if v["so"] is None:
                continue
            for s in [x for x in v["strata"].split(";") if x] + ["all"]:
                conc[b][(s, "read")].append((v["alt"] - v["so"], v["so"]))
                conc[b][(s, "fragment")].append((v["altf"] - v["so"], v["so"]))
            if b == new:
                disc.append(
                    (
                        abs(v["altf"] - v["so"]),
                        tag,
                        *v["key"],
                        v["so"],
                        v["alt"],
                        v["altf"],
                        v["strata"],
                    )
                )
    worst = []
    with open(os.path.join(rep, "concordance_by_stratum.tsv"), "w") as fh:
        fh.write(
            "stratum\tlevel\tn\tmedian_delta_base\tmedian_delta_new\twithin_base\twithin_new\twithin_change\n"
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
            if k[1] == "fragment" and db:
                worst.append((wn - wb, k[0], len(dn)))
    with open(os.path.join(rep, "discordant.tsv"), "w") as fh:
        fh.write(
            "abs_delta_fragment\ttag\t"
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

    # Times
    def times(b):
        t = collections.defaultdict(list)
        for p in glob.glob(os.path.join(out, b, "times.*.tsv")):
            for x in csv.DictReader(open(p), delimiter="\t"):
                t[x["arm"]].append(
                    (
                        float(x["wall_s"]) if x["wall_s"] not in ("NA", "") else None,
                        int(x["max_rss_kb"]) if x["max_rss_kb"] not in ("NA", "") else None,
                    )
                )
        return t

    tb, tn = times(base), times(new)
    with open(os.path.join(rep, "times.tsv"), "w") as fh:
        fh.write(
            "arm\truns\tmedian_wall_base\tmedian_wall_new\twall_ratio\tmax_rss_base_kb\tmax_rss_new_kb\n"
        )
        for arm in sorted(set(tb) | set(tn)):
            wb = [w for w, _ in tb.get(arm, []) if w is not None]
            wn = [w for w, _ in tn.get(arm, []) if w is not None]
            rb = [m for _, m in tb.get(arm, []) if m is not None]
            rn = [m for _, m in tn.get(arm, []) if m is not None]
            mwb = statistics.median(wb) if wb else float("nan")
            mwn = statistics.median(wn) if wn else float("nan")
            fh.write(
                f"{arm}\t{len(wn)}\t{mwb:.1f}\t{mwn:.1f}\t{(mwn / mwb if wb and mwb else float('nan')):.2f}\t"
                f"{max(rb) if rb else 'NA'}\t{max(rn) if rn else 'NA'}\n"
            )

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

    failed = {
        b: [
            x.split("\t")[0]
            for p in glob.glob(os.path.join(out, b, "failed.*.txt"))
            for x in open(p).read().splitlines()
            if x
        ]
        for b in (base, new)
    }
    with open(os.path.join(rep, "gate.txt"), "w") as fh:
        fh.write(f"runs in list: {len(runs)}; missing outputs: {dict(missing) or 0}\n")
        fh.write(f"failed runs: {base} {len(failed[base])}, {new} {len(failed[new])}\n")
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
        worst = sorted(worst)[:10]
        fh.write(
            "strata whose concordance with the sign-out fell most (fragment, within-fraction change):\n"
        )
        for d, s, n in worst:
            fh.write(f"  {s} ({n}): {d:+.3f}\n")
    print(open(os.path.join(rep, "gate.txt")).read())
    print(f"reports in {rep}")


if __name__ == "__main__":
    main()
