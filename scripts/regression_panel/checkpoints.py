"""The attribution checkpoints for a version comparison (docs/development/
regression-panel.md): one build per
count-affecting merge into develop since BASE_TAG (first-parent history), each
interval naming every merge it covers, so a merge skipped as output-neutral is still
named if a change lands in its interval. Decided 2026-10-07: skip only merges whose
own acceptance was byte-identical (or that change no code path that counts).
usage: checkpoints.py REPO BASE_TAG NEW_REF > checkpoints.tsv"""

import subprocess
import sys

# Merges whose own acceptance recorded byte-identical output, or that touch no
# counting code path (version bumps, docs, logging, packaging): not built.
NEUTRAL = {
    "7f1273b4": "back-merge of the 6.5.0 release (version)",
    "9270b4ac": "docs (#158)",
    "e7cfa47e": "trace read names: logging only (#162)",
    "bc090548": "H3 PR B: byte-identical to PR A (#209)",
    "9b4b7c6b": "release infrastructure: version sources (#229)",
    "61f7558c": "GTF loading: byte-identical (#230)",
    "f4b81b3b": "Rust dependency majors: byte-identical (#238)",
    "f73f9219": "bio 4.2.1: byte-identical (#248)",
}


def git(repo, *a):
    return subprocess.run(
        ["git", "-C", repo, *a], capture_output=True, text=True, check=True
    ).stdout


def main():
    repo, base, new = sys.argv[1:4]
    merges = git(
        repo, "log", "--first-parent", "--merges", "--reverse", "--format=%h\t%s", f"{base}..{new}"
    ).splitlines()
    print("order\tlabel\tsha\tinterval_merges")
    pending, n = [], 0
    for line in merges:
        sha, subj = line.split("\t", 1)
        label = subj.replace("Merge pull request ", "").replace(" from msk-access/", " ")
        touched = git(
            repo, "diff", "--name-only", f"{sha}^1", sha, "--", "rust/src", "src/gbcms"
        ).split()
        pending.append(label.split()[0] if label.startswith("#") else sha)
        if not touched or sha in NEUTRAL:
            continue
        n += 1
        print(f"{n}\t{label}\t{sha}\t{','.join(pending)}")
        pending = []
    if pending:
        print(f"# merges after the last checkpoint (neutral, not built): {','.join(pending)}")


if __name__ == "__main__":
    main()
