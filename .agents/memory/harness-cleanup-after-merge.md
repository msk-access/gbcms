---
name: harness-cleanup-after-merge
description: After a ticket's PR merges, delete its isolated harness builds (~1 GB each, rebuildable from BUILT_FROM), gzip raw --trace logs, keep outputs; cargo clean when rust/target bloats.
metadata:
  type: feedback
---

The operator asked (2026-09-30) to make disk cleanup part of the routine: the disk
reached 97% with 20+ isolated acceptance builds (`~/test/gbcms/harness/*/src_*`,
~1 GB each: a git archive + its own .venv + cargo target), 8.6 GB of raw `--trace`
logs, and a 12.5 GiB `rust/target` (the debug profile bloats past 6 GB).

**Why:** acceptance builds pile up one or two per ticket and are only needed until
the PR merges; everything they hold is rebuildable from the commit in their
`BUILT_FROM`.

**How to apply:** once a ticket's PR is merged (and its issue closed):
- delete that ticket's `src_*` trees with `rm -r <explicit path>` (the repo hook
  blocks `rm -rf`; guard on `<tree>/rust/Cargo.toml`). If a tree has local source
  edits (a measurement prototype), save `diff -u` of the edited files as
  `<harness>/<tree>.patch` and check it re-applies to `git archive <sha>` first;
  add the tree to `~/test/gbcms/harness/REMOVED_BUILDS.md`;
- `gzip` raw `--trace` logs (≈400× smaller); keep outputs, scripts and logs, the
  acceptance evidence;
- `cargo clean` in `rust/` when `rust/target` exceeds a few GB.
Leave global caches (pip/uv/cargo registry) and the app's worktree pool alone. See
[[commit-before-review-workflows]].
