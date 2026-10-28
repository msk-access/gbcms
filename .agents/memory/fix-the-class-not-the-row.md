---
name: fix-the-class-not-the-row
description: "A count-rule fix starts from the mechanism, stated without the row that showed it; the row is one instance. Enumerate the class as tests, derive cutoffs from a mechanism, justify any narrowing, judge on data the bug wasn't found in."
metadata:
  type: feedback
---

When a diagnosis hands over failing rows, fix the mechanism that makes them fail, not the rows'
shape. On 2026-10-27 the operator asked that "when we fix bugs, we should think of bigger
picture of universal solution then just changing it for the data".

**Why:** our fix rounds start from rows. An agent writes a rule matching the row's shape,
measures it, and keeps or reverts it. That loop drifts into fitting the data. The judge's
second fix round (#304) showed every sign of it:
- 3 of 6 fixes were reverted;
- a tie rule was narrowed to substitution rows to pass its bars, then reverted anyway;
- two kept fixes are narrower than their principle (the guard-indel fix covers indels only,
  one direction, at substitution rows; the tie fix covers REF-vs-other ties only);
- the third-base switch's 20 bp cut sat where the measurements split, and the patient
  oracle uses the same cut.

The read judge exists to replace per-shape rules with one principle
([[universal-haplotype-judge]]); fixing it row by row regrows the rules path inside it.

**How to apply:**
1. **State the mechanism first**: what about the read, the alleles or the judge makes the
   call wrong, with no reference to the row that showed it. Then derive the rule from the
   principles (count the given allele, judge bases not placement, parsimony). The failing row
   is one instance.
2. **Enumerate the class.** Write the other shapes of the same mechanism as red-first
   synthetic tests (SNV, MNV, indel, delins; repeat and unique; read ends; strands;
   qualities), not only the shape that was seen.
3. **Derive cutoffs from a mechanism**, such as an error rate per telling base, a quality,
   or a read length. Never take them from where a measurement splits or from an oracle's own
   cut. Sweep any constant users can change ([[parameters-are-user-settable]]).
4. **Justify narrowing.** A rule narrower than its mechanism says why. Otherwise record it as
   a known partial, and ticket and measure the general version.
5. **Judge on data the bug wasn't found in**: other arms, other data types, and a truth that
   doesn't share the rule (HG002, the bases themselves, the mates). Tuning until an oracle
   agrees converges on the oracle's definitions.
6. **Treat a revert or a narrowing-to-pass as a signal**: ask whether the mechanism was
   wrong, rather than narrowing further.
7. **Report each fix** as mechanism → general rule → this instance → what else it touches
   ([[holistic-effects-map]]).

Promoted to `docs/reference/read-judgment.md` ("How a change to read judgment is made"; each
register row states its mechanism) and the counting-engine skill. Related:
[[best-not-least-bad]], [[genotyper-not-caller]], [[decisions-worked-examples]].
