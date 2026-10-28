---
name: parameters-are-user-settable
description: "Thresholds like --min-baseq (default 20) are user parameters: truth and oracles must take the run's own value, never hard-code the default, and a decided rule must be checked across the values users choose, not only at the default."
metadata:
  type: feedback
---

On 2026-10-26 a scoring rule for the read judge's bars was written as "reads decidable by
bases ≥ Q20". The operator corrected it: "remember the q20 can be changed by the end user as
it a parameter". `--min-baseq` (and `--min-mapq`, gap penalties, margins exposed as flags) is a
user setting; 20 is only the default.

**Why:** a truth or oracle that hard-codes the default silently measures a different gate from
the run when a user changes it, and a rule tuned only at the default can misbehave at the
values users actually pick (the patient oracle and truth v2's `v2q20` both hard-coded 20).

**How to apply:** truth and oracle scripts take the threshold from the run (its command line or
config), never a constant; write bars as "decidable at the run's `--min-baseq`"; and any
count-changing rule whose behaviour depends on a user parameter gets a sweep over plausible
values (e.g. `--min-baseq` 10/20/30) in its acceptance. MSK runs at 20 (operator, 2026-10-26: "We use
20 at MSK so that is ok"), so the default stays; the sweep checks that a rule holds at the other
values a run may need, and is not a search for a new default. The harness reads the default from `src/gbcms/cli.py` (`truthv2_lib.cli_default_min_baseq`, `xoracle_v2.cli_default_min_baseq`), so its pattern must allow the option to span lines: #329 broke a pattern that didn't (fixed 2026-10-26). Related: [[decisions-worked-examples]],
[[test-across-data-types]], [[giab-truth-arm]].
