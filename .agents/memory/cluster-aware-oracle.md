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

**Containment exception (operator, 2026-10-25, #304):** a read whose sibling allele *contains*
this row's change (an MNV listed beside its component SNV) is this row's ALT, not the
sibling's: it carries X's base with a phased neighbour. The oracle as written scores it SIB
whenever the sibling also places it, which encodes the rejected reading (it made "4,611 ->
165 wrongly credited" look like a judge gain). Score containment by the decided rule; at the
MNV's row, a read carrying only one component is still another allele.

Related: [[pysam-validation-oracle]], [[test-across-data-types]], [[decisions-worked-examples]].
