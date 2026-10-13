"""Attribute every changed cell of a version comparison to the merge(s) where it changed
(docs/development/regression-panel.md; the outputs name runs and loci).

  attribute.py prepare TIER OUT BASE NEW ATTRIB [--pad 100]
      For every run whose MAF output changed between BASE and NEW (a changed cell, or a
      row in one version only): a variant file holding the changed rows and the rows
      near them, ATTRIB/runs.tsv to run each build on, and ATTRIB/changed.tsv. Near =
      within PAD bp of a kept row's span, closed transitively (siblings classify
      together). The arms whose outputs hold Benjamini-Hochberg q-values computed across
      the run's rows (mfsd, rna) keep their full variant files: a reduced run would
      change every q-value.
  attribute.py check TIER OUT BASE NEW ATTRIB
      The reduction must not change a cell: BASE and NEW run on the reduced files
      (run_panel.sh with ATTRIB/runs.tsv into ATTRIB/out) must give the full runs'
      values. gbcms_status_reason is exempt on the context rows only (it names the
      co-annotation group, which the reduction may cut at its edge), never on a changed
      row. Exits 1 when a row differs.
  attribute.py attribute TIER OUT BASE NEW ATTRIB CHECKPOINTS
      Trace each changed cell through BASE and then each checkpoint build in order
      (cp_<sha> under ATTRIB/out), and name each interval where its value changes; a
      row in one version only is traced as its presence, column "(row)". A cell whose
      trail misses a checkpoint's output, or does not end at NEW's value, is
      unattributed. Writes ATTRIB/attribution_cells.tsv, attribution_summary.tsv and
      unattributed.tsv; exits 1 when any cell is unattributed.
The MAF arms are attributed (the vcf arm's counts are the dna arm's, written as VCF).
Every step reads the time records as compare_panel.py does: a run whose latest attempt
failed (it may have left a partial MAF) is left out of every step when it failed in BASE
or NEW, and is named; a reduced run whose latest attempt failed is not evidence (check
fails, attribute treats its output as missing) and is named to re-submit.
"""

import collections
import csv
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_panel import failed_runs  # noqa: E402  (a sibling script, not a package)

KEY = ("Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2")
FULL_ARMS = {"mfsd", "rna"}  # outputs hold BH q-values across the run's rows
MISSING = "<no output>"


def maf(path):
    lines = [x for x in open(path) if not x.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def out_maf(root, build, rid):
    f = glob.glob(os.path.join(root, build, rid, "*.maf"))
    return f[0] if f else None


def runs_of(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def by_key(rows):
    """Rows by locus, with the occurrence number for a locus listed twice."""
    seen, out = collections.Counter(), {}
    for y in rows:
        k = tuple(y.get(c, "") for c in KEY)
        out[k + (seen[k],)] = y
        seen[k] += 1
    return out


def changed_keys(rb, rn):
    """Loci whose cells differ (columns both versions write), or in one version only."""
    out = set(rb) ^ set(rn)
    for k in set(rb) & set(rn):
        x, y = rb[k], rn[k]
        if any(x[c] != y[c] for c in y if c in x):
            out.add(k)
    return out


def span(v):
    c, s = v["Chromosome"], int(v["Start_Position"])
    try:
        e = int(v["End_Position"])
    except (KeyError, ValueError):
        e = s
    return c, min(s, e), max(s, e)


def near(a, b, pad):
    return a[0] == b[0] and max(b[1] - a[2], a[1] - b[2], 0) <= pad


def failed_in(root, builds, what):
    """Runs whose latest attempt failed in any of the builds under ROOT, named on stdout
    (WHAT says what is done with them): {build: set of run ids}."""
    failed = {b: set(failed_runs(root, b)) for b in builds}
    for b, rids in failed.items():
        if rids:
            print(f"{what} ({b}, latest attempt failed): {', '.join(sorted(rids))}")
    return failed


def prepare(tier, out, base, new, attrib, pad):
    os.makedirs(os.path.join(attrib, "mafs"), exist_ok=True)
    keep_runs, n_rows, n_changed, changed_rows = [], 0, 0, []
    failed = failed_in(out, (base, new), "skipped")
    for r in runs_of(os.path.join(tier, "runs.tsv")):
        if r["arm"] == "vcf" or any(r["run_id"] in f for f in failed.values()):
            continue
        fb, fn = out_maf(out, base, r["run_id"]), out_maf(out, new, r["run_id"])
        if not fb or not fn:
            continue
        rb, rn = by_key(maf(fb)), by_key(maf(fn))
        ch = changed_keys(rb, rn)
        if not ch:
            continue
        src = maf(os.path.join(tier, r["variants"]))
        if r["arm"] in FULL_ARMS:
            keep = src
        else:
            seeds = [span(rn.get(k) or rb[k]) for k in ch]
            kept, rest = [], [(span(v), v) for v in src]
            frontier = seeds
            while frontier:
                take = [(s, v) for s, v in rest if any(near(s, f, pad) for f in frontier)]
                rest = [(s, v) for s, v in rest if not any(near(s, f, pad) for f in frontier)]
                kept += take
                frontier = [s for s, _ in take]
            keep = [v for v in src if any(v is w for _, w in kept)]
        path = os.path.join("mafs", f"{r['run_id']}.maf")
        with open(os.path.join(attrib, path), "w") as fh:
            w = csv.DictWriter(
                fh, fieldnames=list(src[0].keys()), delimiter="\t", lineterminator="\n"
            )
            w.writeheader()
            w.writerows(keep)
        keep_runs.append({**r, "variants": os.path.join(os.path.relpath(attrib, tier), path)})
        changed_rows += [(r["run_id"], *k) for k in sorted(ch)]
        n_rows += len(keep)
        n_changed += len(ch)
    cols = ["run_id", "tag", "arm", "mode", "root", "bam_relpath", "variants", "extra_args"]
    with open(os.path.join(attrib, "runs.tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows([{c: r[c] for c in cols} for r in keep_runs])
    with open(os.path.join(attrib, "changed.tsv"), "w") as fh:
        fh.write("run_id\t" + "\t".join(KEY) + "\toccurrence\n")
        for row in changed_rows:
            fh.write("\t".join(map(str, row)) + "\n")
    print(
        f"runs with changes: {len(keep_runs)}; changed rows: {n_changed}; rows to re-run: {n_rows}"
    )


def check(tier, out, base, new, attrib):
    bad = collections.Counter()
    failed = failed_in(out, (base, new), "skipped")
    reduced_failed = failed_in(os.path.join(attrib, "out"), (base, new), "reduced run to re-submit")
    for r in runs_of(os.path.join(attrib, "runs.tsv")):
        if any(r["run_id"] in f for f in failed.values()):
            continue
        fb, fn = out_maf(out, base, r["run_id"]), out_maf(out, new, r["run_id"])
        ch = changed_keys(by_key(maf(fb)), by_key(maf(fn))) if fb and fn else set()
        for build in (base, new):
            if r["run_id"] in reduced_failed[build]:
                bad[f"{build}: a reduced run's latest attempt failed"] += 1
                continue
            full = out_maf(out, build, r["run_id"])
            red = out_maf(os.path.join(attrib, "out"), build, r["run_id"])
            if not full or not red:
                bad[f"{build}: missing output"] += 1
                continue
            f, d = by_key(maf(full)), by_key(maf(red))
            for k, y in d.items():
                x = f.get(k)
                if x is None:
                    bad[f"{build}: a reduced row is not in the full run"] += 1
                    continue
                exempt = {"gbcms_status_reason"} if k not in ch else set()
                if any(x[c] != y[c] for c in y if c in x and c not in exempt):
                    bad[f"{build}: a reduced row differs from the full run"] += 1
    print(dict(bad) or "reduction exact: every reduced row equals the full run's")
    if bad:
        sys.exit(1)


def attribute(tier, out, base, new, attrib, cps_path):
    cps = list(
        csv.DictReader((x for x in open(cps_path) if not x.startswith("# ")), delimiter="\t")
    )
    builds = [base] + [f"cp_{c['sha']}" for c in cps]
    label = {f"cp_{c['sha']}": f"{c['label']} [{c['interval_merges']}]" for c in cps}
    aroot = os.path.join(attrib, "out")
    cells_out, summary, unattr = [], collections.Counter(), []
    failed = failed_in(out, (base, new), "skipped")
    # a failed reduced run's output is missing evidence, whatever it holds
    reduced_failed = failed_in(aroot, builds, "reduced run to re-submit")
    for r in runs_of(os.path.join(attrib, "runs.tsv")):
        if any(r["run_id"] in f for f in failed.values()):
            continue
        fb, fn = out_maf(out, base, r["run_id"]), out_maf(out, new, r["run_id"])
        rb, rn = by_key(maf(fb)), by_key(maf(fn))
        trail = {}
        for b in builds:
            p = out_maf(aroot, b, r["run_id"])
            trail[b] = by_key(maf(p)) if p and r["run_id"] not in reduced_failed[b] else None
        for k in sorted(changed_keys(rb, rn)):
            x, y = rb.get(k), rn.get(k)
            if x is None or y is None:  # a row in one version only: trace its presence
                cols = {"(row)": ("present" if x else "absent", "present" if y else "absent")}
            else:
                cols = {c: (x[c], y[c]) for c in y if c in x and x[c] != y[c]}
            for c, (v_base, v_new) in cols.items():
                vals = []
                for b in builds:
                    t = trail[b]
                    if t is None:
                        vals.append(MISSING)
                    elif c == "(row)":
                        vals.append("present" if k in t else "absent")
                    else:
                        vals.append(t[k].get(c, MISSING) if k in t else MISSING)
                steps = [
                    (builds[i], vals[i - 1], vals[i])
                    for i in range(1, len(builds))
                    if vals[i] != vals[i - 1]
                ]
                if MISSING in vals or vals[-1] != v_new or vals[0] != v_base or not steps:
                    unattr.append(
                        (r["run_id"], r["arm"], *k[:5], c, v_base, v_new, " | ".join(vals))
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
                        v_base,
                        v_new,
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
    if unattr:
        sys.exit(1)


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
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
