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


@pytest.fixture()
def guard(env, monkeypatch):
    """The hook as a module, loaded with the made-up settings, never the operator's."""
    monkeypatch.setenv("HOME", env["HOME"])
    monkeypatch.setenv("HPC_GUARD_CONFIG", env["HPC_GUARD_CONFIG"])
    spec = importlib.util.spec_from_file_location("hpc_ssh_guard", HOOK)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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
    f"ssh -a -x -o BatchMode=yes -o ConnectTimeout=10 {HOST} squeue",
    f"ssh -T -n -a -x {HOST} 'sinfo -s'",
    # Copies and commands that do not reach the HPC are not the guard's business.
    "rsync --version 2>&1",
    "rsync -a /tmp/a/ /tmp/b/",
    "git commit -m 'docs: describe the ssh guard'",
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
    # Everyday reads: a pattern, a count or a format is not a file name.
    f"{SSH} 'squeue -u tester | grep RUNNING'",
    f"{SSH} 'head -n 5 {WORK}/a.txt'",
    f"{SSH} 'cut -d, -f1 {WORK}/a.csv | sort -u | uniq -c'",
    f"{SSH} \"find {WORK} -name '*.out' -type f -newer {WORK}/start\"",
    f"{SSH} 'ls {WORK}/runs/*/done'",
    f"{SSH} 'cd {WORK} && grep -c ALT runs/a.maf'",
    f"{SSH} 'cd {WORK} && ls'",
    sbatch(f"--job-name={PREFIX}o --partition={SHORT} --time=30 --output={WORK}/logs/%x-%j.out"),
    sbatch(f"--job-name={PREFIX}a --partition={SHORT} --time=30", f"{WORK}/a.sbatch {TEST}/x.bam"),
    f"rsync -avz -e 'ssh -a -x -o BatchMode=yes' --exclude='*.tmp' {LOCAL}/kit/ {HOST}:{WORK}/kit/",
    f"scp -p -o BatchMode=yes {LOCAL}/a.txt {HOST}:{TMP}/a.txt",
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
    # A second command the remote shell runs: a newline, a quoted separator, a substitution.
    f"{SSH} 'hostname\nrm {WORK}/obsolete.txt'",
    f"{SSH} squeue\nssh {HOST} 'rm x'",
    f"{SSH} squeue '' ';' 'rm {WORK}/x'",
    f'echo "$(ssh {HOST} rm x)"',
    f"echo `ssh {HOST} rm x`",
    f"{SSH} 'cd {WORK} || cat .bashrc'",
    f"{SSH} 'ls {WORK} & cat .bashrc'",
    # A script or path in the place of a listed command.
    f"{SSH} {WORK}/hostname",
    f"{SSH} './squeue'",
    # Paths the folder check must see: relative, the home folder, '..', braces, variables.
    f"{SSH} 'cd {WORK} && mkdir -p ../../outside'",
    f"{SSH} 'cat $HOME/.bashrc'",
    f"{SSH} 'cat .bashrc'",
    f"{SSH} 'date -f .bashrc'",
    f"{SSH} 'ls'",
    f"{SSH} 'cd {WORK}/missing; cat .bashrc'",
    f"{SSH} 'cd {WORK} | cat .bashrc'",
    f"{SSH} 'ls {{/etc,{WORK}}}'",
    f"{SSH} 'cat {WORK}/../../../etc/passwd'",
    f"{SSH} 'cat ~/.ssh/id_rsa'",
    f"{SSH} 'grep -f/etc/passwd x {WORK}/a'",
    # Read commands that write or run a program.
    f"{SSH} 'sort --compress-program=sh {WORK}/x'",
    f"{SSH} 'sort -o {WORK}/y {WORK}/x'",
    f"{SSH} 'sort --out={WORK}/y {WORK}/x'",
    f"{SSH} 'uniq {WORK}/a {TEST}/b'",
    f"{SSH} 'find {WORK} -fls {TEST}/x'",
    f"{SSH} 'file -C -m {WORK}/magic'",
    # ssh options off the list, wherever they stand.
    f"ssh -a -x -o Hostname=other.example {HOST} hostname",
    f"ssh -a -x -K {HOST} squeue",
    f"ssh -a -x -o GSSAPIDelegateCredentials=yes {HOST} squeue",
    f"ssh -a -x -l other {HOST} squeue",
    f"ssh -a -x {HOST} -A squeue",
    f"X=1 ssh -a -x {HOST} squeue",
    # Copies: the same ssh list, the rsync and scp lists, plain remote paths.
    f"rsync -a -e 'ssh -a -x -A -J jump' {LOCAL}/x {HOST}:{WORK}/",
    f"rsync -a -e 'ssh -a -x {HOST} rm' {LOCAL}/x {HOST}:{WORK}/",
    f"scp -o ProxyCommand=nc {LOCAL}/a {HOST}:{TMP}/a",
    f"scp -A {LOCAL}/a {HOST}:{TMP}/a",
    f"{RSYNC} --rsync-path='sh {WORK}/run.sh' {LOCAL}/x {HOST}:{WORK}/",
    f"/usr/bin/rsync --delete {LOCAL}/in/ {HOST}:{WORK}/out/",
    f"{RSYNC} -M--log-file=/etc/x {LOCAL}/x {HOST}:{WORK}/",
    f"{RSYNC} --temp-dir=/etc {LOCAL}/x {HOST}:{WORK}/",
    f"{RSYNC} '{HOST}:{WORK}/x /etc/passwd' {LOCAL}/",
    f"{RSYNC} '{HOST}:{WORK}/x;id' {LOCAL}/",
    # Other ways to reach the HPC.
    f"git clone {HOST}:{WORK}/repo {LOCAL}/r",
    f"sed '1e ssh {HOST} rm x' {LOCAL}/notes.txt",
    'bash -c "ssh 10.0.0.5 rm x"',
    # sbatch options read as SLURM reads them.
    sbatch(f"--job-name={PREFIX}x --partition={SHORT} --time=60 --output {WORK}/x.sh", "a.sbatch"),
    sbatch(f"--job-name={PREFIX}x --partition={SHORT} --time=60 --wr=hostname"),
    sbatch(f"--job-name={PREFIX}x --partition={SHORT} --time=60 --output=/etc/x"),
    sbatch(f"--job-name={PREFIX}x --partition={SHORT} --time=60 --output={TEST}/x"),
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
def test_slurm_time_values_read_as_minutes(value, minutes, guard):
    assert guard.slurm_minutes(value) == minutes


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


def test_the_example_settings_file_guards(tmp_path):
    example = HOOK.parent / "hpc-guard.example.json"
    env = {**os.environ, "HOME": str(tmp_path), "HPC_GUARD_CONFIG": str(example)}
    host = json.loads(example.read_text())["host"]
    assert decision(f"ssh -a -x {host} squeue", env) == "pass"
    assert decision(f"ssh {host} squeue", env) == "deny"


def without(key):
    return {k: v for k, v in CONFIG.items() if k != key}


@pytest.mark.parametrize(
    "settings",
    [
        *(without(k) for k in ("host", "job_prefix", "remote_roots", "local_roots")),
        *(without(k) for k in ("short_partition", "long_partition")),
        {**CONFIG, "remote_root": [WORK]},  # a misspelt key
        {**CONFIG, "remote_roots": ["work/tester"]},  # a relative folder
        {**CONFIG, "remote_roots": ["/"]},
        {**CONFIG, "local_roots": []},
        {**CONFIG, "deny_patterns": ["("]},
        {**CONFIG, "max_array_tasks": "20"},
        {},
        [],
    ],
)
def test_incomplete_settings_deny_every_ssh(settings, tmp_path):
    """A settings file that is present but incomplete guards nothing, so it denies."""
    config = tmp_path / "guard.json"
    config.write_text(json.dumps(settings))
    env = {**os.environ, "HOME": str(tmp_path), "HPC_GUARD_CONFIG": str(config)}
    assert decision(f"ssh -a -x {HOST} squeue", env) == "deny"
    assert decision("ssh -a -x otherhost rm -rf /", env) == "deny"
    assert decision("ls -la", env) == "pass"


def test_a_missing_named_settings_file_denies_every_ssh(tmp_path):
    env = {**os.environ, "HOME": str(tmp_path), "HPC_GUARD_CONFIG": str(tmp_path / "none.json")}
    assert decision(f"ssh -a -x {HOST} squeue", env) == "deny"
