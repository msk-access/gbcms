---
name: merge-with-gh-when-asked
description: "msk-access/gbcms has GitHub auto-merge disabled; merge a PR with `gh pr merge --merge` when the operator asks, after checking its CI is green; don't wait on or poll CI to merge on your own."
metadata:
  type: feedback
---

On 2026-10-26, auto-merge could not be enabled on #329 and #330 ("Auto merge is not allowed
for this repository"). The operator said "Use gh to merge", then "When i ask".

**Why:** the repository disallows GitHub auto-merge, and CI must not be polled from the
session. The operator decides when a PR lands.

**How to apply:** keep Auto-fix on, so failures, conflicts and review comments wake the
session. When the operator asks to merge, check the PR's checks once with
`gh pr checks <n> --repo msk-access/gbcms`. If they're green, run
`gh pr merge <n> --repo msk-access/gbcms --merge`, a merge commit as for earlier PRs. If
they're still pending, say so and wait to be asked again. Related:
[[commit-before-review-workflows]].
