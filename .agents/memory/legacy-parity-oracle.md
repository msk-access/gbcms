---
name: legacy-parity-oracle
description: count_bam is the feature-gated parity oracle; mirror binned changes in both paths or parity breaks
metadata:
  type: project
---

`count_bam` (per-variant) + its helper `count_single_variant` are the **parity oracle**
for production `count_bam_binned`, gated behind the default `legacy-parity` Cargo
feature. Dev/test builds include it (the parity tests run against it); the **shipped
wheel** builds `--no-default-features` (Dockerfile, release.yml) to omit it. Production
never calls `count_bam`.

**Going-forward rule:** any change to read classification / filtering / fragment
consensus / fetch-window logic in the binned path (`count_bin_shared` /
`count_variant_from_cache`) MUST be mirrored in `count_single_variant`, or the
binned↔legacy parity tests (`count_both`, `test_filters`, `test_parity_large_deletion`)
fail. RNA / mFSD / ASJD / strandedness are exempt (binned-only; not in `PARITY_FIELDS`)
— do not add them to the legacy path. Remove `count_bam` only after parity sign-off;
never let the two paths silently diverge.

Full contract: `.agents/rules/architecture.md` §"Legacy count_bam parity oracle".
Parity holds only *without* siblings — see [[siblings-break-binned-legacy-parity]].
