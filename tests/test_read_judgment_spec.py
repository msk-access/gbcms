"""The read-judgment spec, executed (docs/reference/read-judgment.md).

Every case in `read_judgment_cases.py` must get its expected call. A decided case
fails when the engine breaks its rule; an open case fails when today's call
changes before its decision lands. Either way the change is a decision, made with
the operator before an expectation is edited (AGENTS.md).
"""

import pytest
from read_judgment_cases import DECISIONS, EXPECT, cases, run

CASES = cases()


def test_every_case_has_an_expectation_and_a_decision():
    assert {c.id for c in CASES} == set(EXPECT)
    assert {c.group for c in CASES} == set(DECISIONS)


@pytest.mark.parametrize("case", CASES, ids=lambda c: c.id)
def test_the_case_gets_its_call(tmp_path, case):
    got, verdicts = run(case, tmp_path)
    want = EXPECT[case.id]
    assert got == want, f"{case.id} [{DECISIONS[case.group]}]: got {got}, want {want} ({verdicts})"
