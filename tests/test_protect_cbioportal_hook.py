"""The cBioPortal guard's Bash rule (`.claude/hooks/protect-cbioportal.py`).

Tool output is sent to the model provider, so a command that names protected study data passes
only when both of its output streams, stdout and stderr, end in files under a protected root
(no pipes or chains). The hook runs as it does in Claude Code: a process reading the tool call
as JSON on stdin, with the protected roots listed in ~/.claude/protected-paths.txt (here a
temporary home). The file names are made up.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

HOOK = Path(__file__).resolve().parents[1] / ".claude" / "hooks" / "protect-cbioportal.py"
DATA = "data_clinical_patient.txt"


@pytest.fixture
def run(tmp_path):
    home = tmp_path / "home"
    (home / ".claude").mkdir(parents=True)
    root = tmp_path / "protected"
    root.mkdir()
    (home / ".claude" / "protected-paths.txt").write_text(f"{root}\n")

    def denied(command: str) -> bool:
        p = subprocess.run(
            [sys.executable, "-I", str(HOOK)],
            input=json.dumps(
                {"tool_name": "Bash", "tool_input": {"command": command.format(p=root)}}
            ),
            capture_output=True,
            text=True,
            env={**os.environ, "HOME": str(home)},
            check=True,
        )
        out = p.stdout.strip()
        return (
            bool(out)
            and json.loads(out).get("hookSpecificOutput", {}).get("permissionDecision") == "deny"
        )

    return denied


@pytest.mark.parametrize(
    "command",
    [
        f"python3 -I s.py {DATA} > {{p}}/out.txt 2> {{p}}/err.txt",
        f"python3 -I s.py {DATA} > {{p}}/out.txt 2>&1",
        f"python3 -I s.py {DATA} >> {{p}}/out.txt 2>> {{p}}/err.txt",
        f"python3 -I s.py {DATA} 1> {{p}}/out.txt 2> {{p}}/err.txt",
        f"python3 -I s.py {DATA} &> {{p}}/all.txt",
    ],
)
def test_both_streams_into_a_protected_file_pass(run, command):
    assert not run(command)


@pytest.mark.parametrize(
    "command",
    [
        # stderr alone into the protected root: stdout prints the rows (the 2> is not stdout)
        f"cat {DATA} 2> {{p}}/err.txt",
        f"cat {DATA} 2>> {{p}}/err.txt",
        # stdout alone: stderr reaches the console
        f"python3 -I s.py {DATA} > {{p}}/out.txt",
        # stderr joined to stdout before stdout is redirected: stderr reaches the console
        f"python3 -I s.py {DATA} 2>&1 > {{p}}/out.txt",
        # one stream outside the protected root
        f"python3 -I s.py {DATA} > {{p}}/out.txt 2> /tmp/err.txt",
        f"python3 -I s.py {DATA} > /tmp/out.txt 2> {{p}}/err.txt",
        f"python3 -I s.py {DATA} &> /tmp/all.txt",
        f"python3 -I s.py {DATA} > {{p}}/out.txt 2> {{p}}/err.txt > /tmp/again.txt",
        # stdout sent to stderr, or no redirect at all
        f"cat {DATA} >&2 2> {{p}}/err.txt",
        f"cat {DATA}",
        # a pipe or a chain
        f"cat {DATA} | head > {{p}}/out.txt 2> {{p}}/err.txt",
        f"cat {DATA} > {{p}}/out.txt 2> {{p}}/err.txt; cat {{p}}/out.txt",
    ],
)
def test_a_stream_that_can_reach_the_console_is_denied(run, command):
    assert run(command)


def test_a_command_that_names_no_protected_data_is_not_stopped(run):
    assert not run("ls -la /tmp")
