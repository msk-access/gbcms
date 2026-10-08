---
name: github-sub-issues
description: Group related GitHub work as sub-issues (parent + children), not combined issues or checklists
metadata:
  type: feedback
---

When filing related GitHub issues (a cycle's tickets, an umbrella issue's items, a
cluster of small tickets), use sub-issues wherever possible: a parent issue with one
sub-issue per ticket, nested as needed, instead of one combined issue or a
checklist of items. Each cycle gets a tracking issue with every ticket under it.

**Why:** operator preference (2026-09-25): sub-issues make each item trackable and
closable on its own while the parent shows progress.

**How to apply:** `gh issue create --parent N` for new issues; `gh issue edit N
--parent P` / `--add-sub-issue` for existing ones (gh ≥ 2.9x). The 6.6.0 tree is
tracker #140 → #92, #112, #133, #134, #135 with their own sub-issues.
