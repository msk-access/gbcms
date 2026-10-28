---
name: read-rule-register-to-the-end
description: "Before stating what a decided read rule says, read its whole row in docs/reference/read-judgment.md, where amendments are appended at the end, and the contract test's docstring; the reference pages (allele-classification.md, glossary.md) can lag an amendment."
metadata:
  type: feedback
---

On 2026-10-26 I told the operator that the read judge "breaks the insertion readability
rule (RJ-20)": an insert of unreadable bases should be neither + partial, not ALT. I
queued a fix to "port" it, recorded it on #304, and wrote it into #329's Read Filters
text. All of it was wrong.

RJ-20's register row states the old rule first and ends with "Amended 2026-10-15 (RJ-27):
an unreadable insert of the ALT's length inside the discrimination window is the ALT". That
was the operator's decision on #278. The judge and the rules path both follow it, and
`tests/test_unreadable_insert_contract.py` pins it. I had read the stale
`allele-classification.md` §Insertion and the start of the row. The fix round found the
mistake, and the contract tests would have blocked a change.

**Why:** the register appends amendments at the end of a row, and the reference pages were
not updated when the rule changed. So the first text found is the superseded rule.

**How to apply:**
- Before asserting a decided rule, read its register row to the end, then the contract
  test's docstring, which names the current decision.
- Treat `allele-classification.md` and `glossary.md` as possibly behind the register.
- When a rule is amended, put the current rule at the start of its register row and
  update the reference pages in the same change.

Related: [[decisions-worked-examples]], [[test-changes-need-operator-notice]].
