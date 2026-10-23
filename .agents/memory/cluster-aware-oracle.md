---
name: cluster-aware-oracle
description: "Score a change at co-annotated rows against every allele of the cluster (REF, the row's ALT, each sibling's ALT); a REF-vs-ALT-only read check calls a removed sibling molecule \"lost\"."
metadata:
  node_type: memory
  type: reference
  originSessionId: a8529157-ff34-4400-8321-c28f0357971d
  modified: 2026-10-23T09:46:12.339Z
---

The plan670 read check `m294/oracle_changed.py` aligns each changed read only to the
row's REF and ALT. At co-annotated rows that misleads: a sibling allele's molecule
often fits this row's ALT better than REF, so removing it from AD shows as "lost".
On #300 (sibling guard by bases) it reported 157 panel and 204 gate reads lost;
`m294/oracle_sibling.py`, which also aligns to each sibling's ALT (rows of the
same sample within 30 bp), found 556 of the removals were the sibling's molecules
or ties, and only 17 gained calls were wrong.

Use `oracle_sibling.py <arm> <label A> <label B>` for any change that can move a read
between co-annotated rows; `oracle_changed.py` stays right for isolated rows.
(`oracle_sibling.py` imports `oracle_changed.py`, which runs its own check at
import, so a run takes twice as long.)

Related: [[pysam-validation-oracle]], [[test-across-data-types]], [[decisions-worked-examples]].
