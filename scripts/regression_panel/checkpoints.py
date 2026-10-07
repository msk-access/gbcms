"""The attribution checkpoints for a version comparison (docs/development/
regression-panel.md): one build per count-affecting commit on develop's first-parent
history since BASE_TAG (in practice each merged PR), each interval naming every commit
it covers, so a commit skipped as output-neutral is still named if a change lands in
its interval. Decided 2026-10-07: skip only merges whose own acceptance was
byte-identical (or that change no code path that counts).

Two rows frame the trail:
- row 0 is BASE_TAG itself, built as a wheel like every checkpoint: a cell that
  changes between the BASE image's output and row 0 is a build difference (image
  against wheel), not a merge;
- the last row is NEW_REF, the commit the release ships, even when the commits after
  the last count-affecting one touch nothing that counts.
usage: checkpoints.py REPO BASE_TAG NEW_REF > checkpoints.tsv"""

import subprocess
import sys

# Commits whose own acceptance recorded byte-identical output, or that touch no
# counting code path (version bumps, docs, logging, packaging): not built. Matched by
# prefix, so a clone's abbreviation does not matter.
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


def neutral(sha):
    return any(sha.startswith(k) for k in NEUTRAL)


def main():
    repo, base, new = sys.argv[1:4]
    base_sha = git(repo, "rev-parse", f"{base}^{{commit}}").strip()
    new_sha = git(repo, "rev-parse", f"{new}^{{commit}}").strip()
    commits = git(
        repo, "log", "--first-parent", "--reverse", "--format=%H\t%s", f"{base}..{new}"
    ).splitlines()
    print("order\tlabel\tsha\tinterval_merges")
    print(f"0\t{base} (as a wheel: build environment)\t{base_sha}\t{base}")
    pending, n, last = [], 0, None
    for line in commits:
        sha, subj = line.split("\t", 1)
        label = subj.replace("Merge pull request ", "").replace(" from msk-access/", " ")
        pending.append(label.split()[0] if label.startswith("#") else sha[:8])
        touched = git(
            repo, "diff", "--name-only", f"{sha}^1", sha, "--", "rust/src", "src/gbcms"
        ).split()
        if not touched or neutral(sha):
            continue
        n += 1
        print(f"{n}\t{label}\t{sha}\t{','.join(pending)}")
        pending, last = [], sha
    if last != new_sha:
        n += 1
        print(
            f"{n}\t{new} (the candidate; no counting change since the previous row)\t{new_sha}\t"
            f"{','.join(pending) or new}"
        )


if __name__ == "__main__":
    main()
