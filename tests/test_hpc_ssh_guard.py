"""The HPC guard hook (.claude/hooks/hpc-ssh-guard.py) decides as the HPC rules say.

Each case is a command the model might run and the decision the hook must give: "pass"
(no objection), "ask" (the operator approves it) or "deny". The settings are made-up
examples given through $HPC_GUARD_CONFIG, with HOME pointed elsewhere, so an operator's own
settings never reach the tests.
"""

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

HOOK = Path(__file__).resolve().parents[1] / ".claude" / "hooks" / "hpc-ssh-guard.py"
HOST = "hpc-test"
WORK = "/work/tester/project"
TMP = "/scratch/tester"
TEST = "/work/tester"  # the read-only root
LOCAL = "/local/harness"
SHORT, LONG, PREFIX = "short_q", "long_q", "proj-"
SSH = f"ssh -a -x {HOST}"
RSYNC = 'rsync -a -e "ssh -a -x"'
CONFIG = {
    "host": HOST,
    "hpc_names": [HOST, "login01.example.org"],
    "remote_roots": [WORK, TMP],
    "read_roots": [TEST],
    "local_roots": [LOCAL],
    "short_partition": SHORT,
    "long_partition": LONG,
    "short_limit_minutes": 180,
    "job_prefix": PREFIX,
    "max_array_tasks": 20,
    "deny_patterns": ["cbio"],
}


@pytest.fixture(scope="module")
def env(tmp_path_factory):
    home = tmp_path_factory.mktemp("home")
    config = home / "guard.json"
    config.write_text(json.dumps(CONFIG))
    return {**os.environ, "HOME": str(home), "HPC_GUARD_CONFIG": str(config)}


def run(hook, command, env):
    return subprocess.run(
        [sys.executable, "-I", str(hook)],
        input=json.dumps({"tool_input": {"command": command}}),
        capture_output=True,
        text=True,
        env=env,
    )


def decision(command, env):
    done = run(HOOK, command, env)
    if not done.stdout.strip():
        assert done.returncode == 0, done.stderr
        return "pass"
    return json.loads(done.stdout)["hookSpecificOutput"]["permissionDecision"]


def sbatch(opts, script=f"{WORK}/kit/a.sbatch"):
    return f"{SSH} 'sbatch {opts} {script}'"


PASS = [
    f"{SSH} squeue -u tester -n {PREFIX}pilot",
    f"{SSH} 'hostname; whoami'",
    f"{SSH} 'sinfo -s' 2>&1 | head -5",
    f"{SSH} 'sacct -n {PREFIX}pilot --format=JobID,State'",
    f"ssh -O check {HOST}",
    "ls -la",
    "git push origin feature/x",
    "dig +short example.org SRV",
    # Local tools that only mention the host in their text.
    f"grep -n {HOST} tests/test_hpc_ssh_guard.py",
    f"sed -n '/{HOST}/p' docs/harness/measurement-harness.md",
]

ASK = [
    f"{SSH} 'df -h {WORK} {TMP}'",
    f"{SSH} 'ls {WORK}'",
    f"{SSH} 'mkdir -p {WORK}/runs'",
    f"{SSH} 'cd {WORK} && ls runs'",
    sbatch(f"--job-name={PREFIX}pilot --partition={SHORT} --time=02:00:00 --array=1-133%20"),
    sbatch(f"--job-name={PREFIX}long --partition={LONG} --time=2-00:00:00"),
    sbatch(f"-J {PREFIX}t -p {SHORT} -t 179"),
    f"{SSH} 'scancel --name={PREFIX}pilot'",
    f"{RSYNC} {LOCAL}/kit/ {HOST}:{WORK}/kit/",
    f"rsync -a -e 'ssh -a -x' {HOST}:{WORK}/runs/x/ {LOCAL}/out/x/",
    f"scp {LOCAL}/a.txt {HOST}:{TMP}/a.txt",
    f"{SSH} 'ls {TEST}'",  # reads may reach the read-only root
    f"{SSH} 'du -sh {TEST}/kit'",
    f"{RSYNC} {HOST}:{TEST}/kit/README.md {LOCAL}/",  # so may a copy from it
]

DENY = [
    f"ssh {HOST} squeue",  # forwarding not turned off
    f"ssh -x {HOST} squeue",
    f"ssh -a -x -A {HOST} squeue",
    "ssh -a -x otherhost squeue",
    f"ssh -a -x {HOST}",  # interactive shell
    f"ssh -a -x -L 8080:localhost:80 {HOST} squeue",
    f"ssh -a -x -J jump {HOST} squeue",
    f"ssh -a -x -o ProxyCommand=nc {HOST} squeue",
    f"ssh -a -x -o ForwardAgent=yes {HOST} squeue",
    f"ssh -O exit {HOST}",
    f"{SSH} 'rm -rf {WORK}/x'",
    f"{SSH} 'ls /data1/clinical'",
    f"{SSH} 'cat /home/tester/.bashrc'",
    f"{SSH} 'ls ~/'",
    f"{SSH} 'echo $(whoami)'",
    f"{SSH} 'squeue > {WORK}/x'",
    f"{SSH} 'python3 x.py'",
    f"{SSH} 'bash {WORK}/kit/run.sh'",
    sbatch(f"--partition={SHORT} --time=02:00:00"),  # no job name
    sbatch(f"--job-name={PREFIX}x --partition={SHORT} --time=02:00:00 --array=1-200"),
    sbatch(f"--job-name={PREFIX}x --partition={SHORT} --time=02:00:00 --array=1-200%50"),
    f"{SSH} 'sbatch --job-name={PREFIX}x --partition={SHORT} --time=01:00:00 --wrap=hostname'",
    sbatch(f"--job-name={PREFIX}x --partition=other --time=02:00:00"),
    sbatch(f"--job-name={PREFIX}x --partition={SHORT} --time=02:00:00", "a.sbatch"),
    sbatch(f"--job-name={PREFIX}x --partition={SHORT}"),  # no --time
    sbatch(f"--job-name={PREFIX}x --partition={SHORT} --time=04:00:00"),
    sbatch(f"--job-name={PREFIX}x --partition={LONG} --time=01:00:00"),
    sbatch(f"--job-name={PREFIX}x --partition={SHORT} --time=03:00:01"),
    sbatch(f"--job-name={PREFIX}x --partition={SHORT} --time=02:00:00", f"{TEST}/a.sbatch"),
    f"{SSH} 'scancel 12345'",
    f"{SSH} 'find {WORK} -delete'",
    f"{SSH} 'find {WORK} -exec cat {{}} +'",
    f"{SSH} 'ls {WORK}/cbioportal'",
    f"{SSH} 'mkdir -p {TEST}/x'",  # writes stay in the work folders
    f"{SSH} 'ls /work'",  # above the read-only root
    f"rsync -a {LOCAL}/x {HOST}:{WORK}/",  # no -e "ssh -a -x"
    f"{RSYNC} --delete {LOCAL}/x {HOST}:{WORK}/",
    f"{RSYNC} /Users/someone/Documents/x {HOST}:{WORK}/",
    f"{RSYNC} relative/path {HOST}:{WORK}/",
    f"{RSYNC} {LOCAL}/x {HOST}:/data1/clinical/",
    f"{RSYNC} {LOCAL}/x {HOST}:{TEST}/x",  # a copy to HPC stays in the work folders
    f"{RSYNC} {LOCAL}/x otherhost:/tmp/",
    f"scp {LOCAL}/x {HOST}:/home/tester/x",
    f"sftp {HOST}",
    f'bash -c "ssh -a -x {HOST} squeue"',
    f"python3 -c \"import os; os.system('ssh {HOST} ls')\"",
    f"echo 'ssh -a -x {HOST} squeue' | bash",
    f"env X=1 ssh {HOST} squeue",
    f"xargs -I{{}} ssh -a -x {HOST} {{}} < cmds.txt",
    f"{SSH} 'squeue' > {LOCAL}/out.txt",
]


@pytest.mark.parametrize("command", PASS)
def test_status_and_unrelated_commands_pass(command, env):
    assert decision(command, env) == "pass"


@pytest.mark.parametrize("command", ASK)
def test_reads_copies_and_submissions_ask_the_operator(command, env):
    assert decision(command, env) == "ask"


@pytest.mark.parametrize("command", DENY)
def test_everything_outside_the_rules_is_denied(command, env):
    assert decision(command, env) == "deny"


@pytest.mark.parametrize(
    "value,minutes",
    [
        ("90", 90),
        ("02:00:00", 120),
        ("03:00:00", 180),
        ("03:00:01", 181),
        ("2-00:00:00", 2880),
        ("1-12", 2160),
        ("45:30", 46),
        ("bad", None),
        ("", None),
    ],
)
def test_slurm_time_values_read_as_minutes(value, minutes):
    spec = importlib.util.spec_from_file_location("hpc_ssh_guard", HOOK)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.slurm_minutes(value) == minutes


def test_without_settings_the_hook_guards_nothing(tmp_path):
    """No settings on the machine means no HPC to guard: nothing changes."""
    env = {k: v for k, v in os.environ.items() if k != "HPC_GUARD_CONFIG"}
    env["HOME"] = str(tmp_path)
    hooks = tmp_path / ".claude" / "hooks"
    hooks.mkdir(parents=True)
    (hooks / HOOK.name).write_text(HOOK.read_text())
    done = run(hooks / HOOK.name, "ssh otherhost ls /", env)
    assert done.returncode == 0 and not done.stdout.strip()


def test_an_unreadable_settings_file_denies_every_ssh(tmp_path):
    config = tmp_path / "broken.json"
    config.write_text("{ not json")
    env = {**os.environ, "HOME": str(tmp_path), "HPC_GUARD_CONFIG": str(config)}
    done = run(HOOK, f"ssh -a -x {HOST} squeue", env)
    assert json.loads(done.stdout)["hookSpecificOutput"]["permissionDecision"] == "deny"


def test_the_example_settings_file_parses():
    example = HOOK.parent / "hpc-guard.example.json"
    assert set(json.loads(example.read_text())) == set(CONFIG)
