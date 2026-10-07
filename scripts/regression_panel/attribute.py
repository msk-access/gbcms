"""Attribute every changed cell of a version comparison to the merge(s) where it changed
(docs/development/regression-panel.md; the outputs name runs and loci).

  attribute.py prepare TIER OUT BASE NEW ATTRIB [--pad 100]
      For every run whose MAF output changed between BASE and NEW: a reduced variant
      file holding the changed rows and every row within PAD bp of one (siblings
      classify together), and ATTRIB/runs.tsv to run each checkpoint build on.
  attribute.py check TIER OUT BASE NEW ATTRIB
      The reduction must not change a count: BASE and NEW run on the reduced files
      (run_panel.sh with ATTRIB/runs.tsv into ATTRIB/out) must give
      the full runs' values for every changed row.
  attribute.py attribute TIER OUT BASE NEW ATTRIB CHECKPOINTS
      Trace each changed cell through BASE, then each checkpoint build in order (named
      cp_<sha> under ATTRIB/out), and name each interval where its value changes.
      Writes ATTRIB/attribution_cells.tsv, attribution_summary.tsv, unattributed.tsv.
The MAF arms are attributed (the vcf arm's counts are the dna arm's, written as VCF).
"""

import collections
import csv
import glob
import os
import sys

KEY = ("Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2")


def maf(path):
    lines = [x for x in open(path) if not x.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def out_maf(root, build, rid):
    f = glob.glob(os.path.join(root, build, rid, "*.maf"))
    return f[0] if f else None


def runs_of(tier):
    with open(os.path.join(tier, "runs.tsv")) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def changed_rows(rb, rn):
    """Indexes of rows whose cells differ (columns both builds write)."""
    out = []
    for i, (x, y) in enumerate(zip(rb, rn, strict=True)):
        if any(x[c] != y[c] for c in y if c in x):
            out.append(i)
    return out


def prepare(tier, out, base, new, attrib, pad):
    os.makedirs(os.path.join(attrib, "mafs"), exist_ok=True)
    keep_runs, n_rows, n_changed = [], 0, 0
    for r in runs_of(tier):
        fb, fn = out_maf(out, base, r["run_id"]), out_maf(out, new, r["run_id"])
        if not fb or not fn:
            continue
        rb, rn = maf(fb), maf(fn)
        idx = changed_rows(rb, rn)
        if not idx:
            continue
        src = maf(os.path.join(tier, r["variants"]))
        assert len(src) == len(rn), f"{r['run_id']}: output rows do not follow input rows"
        sites = {(rn[i]["Chromosome"], int(rn[i]["Start_Position"])) for i in idx}
        keep = [
            v
            for v in src
            if any(
                v["Chromosome"] == c and abs(int(v["Start_Position"]) - p) <= pad for c, p in sites
            )
        ]
        path = os.path.join("mafs", f"{r['run_id']}.maf")
        with open(os.path.join(attrib, path), "w") as fh:
            w = csv.DictWriter(
                fh, fieldnames=list(src[0].keys()), delimiter="\t", lineterminator="\n"
            )
            w.writeheader()
            w.writerows(keep)
        keep_runs.append({**r, "variants": os.path.join(os.path.relpath(attrib, tier), path)})
        n_rows += len(keep)
        n_changed += len(idx)
    with open(os.path.join(attrib, "runs.tsv"), "w") as fh:
        cols = ["run_id", "tag", "arm", "mode", "root", "bam_relpath", "variants", "extra_args"]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows([{c: r[c] for c in cols} for r in keep_runs if r["arm"] != "vcf"])
    print(
        f"runs with changes: {len(keep_runs)}; changed rows: {n_changed}; rows to re-run: {n_rows}"
    )


def by_key(rows):
    """Rows by locus, with the occurrence number for loci listed twice (cohort MAFs)."""
    seen, out = collections.Counter(), {}
    for y in rows:
        k = tuple(y[c] for c in KEY)
        out[k + (seen[k],)] = y
        seen[k] += 1
    return out


def check(tier, out, base, new, attrib):
    bad = collections.Counter()
    for r in csv.DictReader(open(os.path.join(attrib, "runs.tsv")), delimiter="\t"):
        for build in (base, new):
            full = out_maf(out, build, r["run_id"])
            red = out_maf(os.path.join(attrib, "out"), build, r["run_id"])
            if not full or not red:
                bad[f"{build}: missing"] += 1
                continue
            f, d = by_key(maf(full)), by_key(maf(red))
            for k, y in d.items():
                x = f.get(k)
                if x is None or any(
                    x[c] != y[c] for c in y if c in x and c not in ("gbcms_status_reason",)
                ):
                    bad[f"{build}: a reduced row differs from the full run"] += 1
    print(dict(bad) or "reduction exact: every reduced row equals the full run's")


def attribute(tier, out, base, new, attrib, cps_path):
    cps = list(csv.DictReader((x for x in open(cps_path) if not x.startswith("# ")), delimiter="	"))
    builds = [base] + [f"cp_{c['sha']}" for c in cps]
    label = {f"cp_{c['sha']}": f"{c['label']} [{c['interval_merges']}]" for c in cps}
    aroot = os.path.join(attrib, "out")
    cells_out, summary, unattr = [], collections.Counter(), []
    for r in csv.DictReader(open(os.path.join(attrib, "runs.tsv")), delimiter="\t"):
        fb, fn = out_maf(out, base, r["run_id"]), out_maf(out, new, r["run_id"])
        rb, rn = by_key(maf(fb)), by_key(maf(fn))
        trail = {}
        for b in builds:
            p = out_maf(aroot, b, r["run_id"])
            trail[b] = by_key(maf(p)) if p else None
        for k, y in rn.items():
            x = rb.get(k)
            if x is None:
                continue
            for c in y:
                if c not in x or x[c] == y[c]:
                    continue
                vals = []
                for b in builds:
                    t = trail[b]
                    vals.append(t.get(k, {}).get(c) if t is not None else None)
                steps = [
                    (builds[i], vals[i - 1], vals[i])
                    for i in range(1, len(builds))
                    if vals[i] is not None and vals[i - 1] is not None and vals[i] != vals[i - 1]
                ]
                if vals[-1] != y[c] or not steps:
                    unattr.append(
                        (
                            r["run_id"],
                            r["arm"],
                            *k[:5],
                            c,
                            x[c],
                            y[c],
                            "|".join(str(v) for v in vals),
                        )
                    )
                    continue
                for b, _v0, _v1 in steps:
                    summary[(label[b], r["arm"], c)] += 1
                cells_out.append(
                    (
                        r["run_id"],
                        r["arm"],
                        *k[:5],
                        c,
                        x[c],
                        y[c],
                        " ; ".join(f"{label[b]}: {v0} -> {v1}" for b, v0, v1 in steps),
                    )
                )
    with open(os.path.join(attrib, "attribution_cells.tsv"), "w") as fh:
        fh.write("run_id\tarm\t" + "\t".join(KEY) + "\tcolumn\tbase\tnew\tchanged_at\n")
        for row in cells_out:
            fh.write("\t".join(map(str, row)) + "\n")
    with open(os.path.join(attrib, "attribution_summary.tsv"), "w") as fh:
        fh.write("checkpoint\tarm\tcolumn\tcells\n")
        for (lab, arm, c), n in sorted(summary.items()):
            fh.write(f"{lab}\t{arm}\t{c}\t{n}\n")
    with open(os.path.join(attrib, "unattributed.tsv"), "w") as fh:
        fh.write("run_id\tarm\t" + "\t".join(KEY) + "\tcolumn\tbase\tnew\tvalues_through_builds\n")
        for row in unattr:
            fh.write("\t".join(map(str, row)) + "\n")
    per_cp = collections.Counter()
    for (lab, _a, _c), n in summary.items():
        per_cp[lab] += n
    print(f"changed cells attributed: {len(cells_out)}; unattributed: {len(unattr)}")
    for lab, n in per_cp.most_common():
        print(f"  {n:7d}  {lab}")


def main():
    cmd = sys.argv[1]
    a = sys.argv[2:]
    if cmd == "prepare":
        pad = int(a[a.index("--pad") + 1]) if "--pad" in a else 100
        prepare(*a[:5], pad)
    elif cmd == "check":
        check(*a[:5])
    elif cmd == "attribute":
        attribute(*a[:6])
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
