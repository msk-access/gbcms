---
description: gbcms security posture — destructive ops, secrets, outward-facing actions. Enforced by hooks, not the model's judgment.
alwaysApply: true
---

# Security Posture

A rule states intent; a **hook enforces it** (a rule in a config competes with the
model's habits and loses under load). The deterministic floor lives in
`.claude/hooks/` (AI layer) and `.githooks/` (tool-agnostic). See
`docs/harness/00-plan.md` and `.agents/skills/` for wiring.

## Destructive ops
- No `rm -rf` — narrow the target or use the session scratchpad. (Enforced:
  `.claude/hooks/block-rm-rf.py`.)
- No direct push to `main`; Git Flow is `feature/* → develop → main` via PR.
  (Enforced: `.claude/hooks/block-push-main.py` + `.githooks/pre-push`.)
- Ask before bulk/irreversible operations.

## Secrets
- Never commit secrets, API keys, credentials, or PHI/patient identifiers — this
  is a clinical-genomics tool; treat sample IDs and paths as sensitive.
- Keep keys in a secret manager, not in `.env` files or hardcoded in scripts.
  (Enforced: `.claude/hooks/block-env-edits.py`.)
- Every push is secret-scanned before it leaves the machine. (Enforced:
  `.githooks/pre-push` runs gitleaks.)
- Scan tracked **code and fixtures**, not only data files — secrets hide in test
  fixtures and hardcoded example lists.

## Outward-facing actions
- **Ask before public.** The default destination for any artifact is private. A
  public host, a deploy, a published file — each needs explicit, per-artifact
  approval that names the destination. (Reflex: making a private repo public
  leaves pre-scrub secrets in git *history* — build a clean export, don't flip.)
- **No comms without approval** — no email/message/PR-to-public goes out without
  showing the draft and getting an explicit "send it."

## Injection is the real threat, not fat-fingering
With broad Bash + live web access, the realistic threat is an *injected*
instruction (web page, repo file, MCP response) chaining an allowed command into
exfiltration. Keep `curl`/`wget` off any blanket allowlist. (Defense-in-depth:
`.claude/hooks/injection-guard.py` asks on the anomalous remote-pipe / secret-exfil shapes.)
