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
    "allow_in_remote_roots": False,
    "ask_patterns": [],
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
    # Text in single quotes is printed, not run.
    f"echo '$(ssh {HOST} rm x)'",
    f"grep -n '`ssh {HOST}' {LOCAL}/notes.md",
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
    f"{SSH} 'head --lines 5 {WORK}/a.txt'",  # a long option's separate value
    f"{SSH} 'tail --lin 5 {WORK}/a.txt'",  # a unique prefix, as GNU tools take
    f"{SSH} 'grep --max-count 3 ALT {WORK}/a.maf'",
    f"{SSH} 'grep --context=3 ALT {WORK}/a.maf'",
    f"{SSH} 'file --magic-file={WORK}/m1:{WORK}/m2 {WORK}/a'",  # a list, every name inside
    f"{SSH} 'cut -d, -f1 {WORK}/a.csv | sort -u | uniq -c'",
    f"{SSH} \"find {WORK} -name '*.out' -type f -newer {WORK}/start\"",
    f"{SSH} 'ls {WORK}/runs/*/done'",
    f"{SSH} 'cd {WORK} && grep -c ALT runs/a.maf'",
    f"{SSH} 'cd {WORK} && ls'",
    f"{SSH} 'mkdir -m 755 {WORK}/x'",
    f"{SSH} 'squeue -u tester | uniq -c -'",
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
    f'echo "it\'s $(ssh {HOST} rm x)"',  # a quote inside double quotes is literal
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
    f"{SSH} 'uniq - {TEST}/out'",  # standard input is the first operand
    f"{SSH} 'uniq {WORK}/a 5'",  # a number left over is a file name
    f"{SSH} 'grep 5 /etc/passwd'",  # the pattern is 5; the file is checked
    f"{SSH} 'grep --file .patterns {WORK}/a'",  # a long option whose value is a file
    f"{SSH} 'grep --file=notes {WORK}/a'",
    f"{SSH} 'file --magic-file /etc/magic {WORK}/a'",
    f"{SSH} 'head --lines 5 .bashrc'",
    f"{SSH} 'ls --hide /etc'",  # a path-like value is checked whatever the option
    # A value that may be optional must not shift a file into the pattern's place.
    f"{SSH} 'grep --context 3 /etc/passwd'",
    f"{SSH} 'grep -C 3 /etc/passwd'",
    f"{SSH} 'grep /etc/passwd'",  # a path-like pattern is checked as a path
    # A list of files: every name is checked.
    f"{SSH} 'file --magic-file={WORK}/magic:/etc/magic {WORK}/a'",
    f"{SSH} 'file -m {WORK}/magic:/etc/magic {WORK}/a'",
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


# The opt-in setting allow_in_remote_roots: what stays in the work folders runs without asking.

IN_WORK = [
    # Reads in the work folders, by path or after a cd into them.
    f"{SSH} 'ls {WORK}'",
    f"{SSH} 'tail -n 50 {TMP}/logs/a.out'",
    f"{SSH} 'cd {WORK} && grep -c ALT runs/a.maf'",
    f"{SSH} 'cd {WORK} && ls'",
    f"{SSH} 'sha256sum {WORK}/kit/a.bam {TMP}/b.bam'",
    f"{SSH} \"find {WORK} -name '*.out' -newer {WORK}/start\"",
    f"{SSH} 'squeue -u tester | grep RUNNING'",  # a read of standard input names no path
    f"{SSH} 'cat {WORK}/a.tsv' | grep ALT | head -5",  # local text tools after it
    f"ssh -O check {HOST} && {SSH} 'ls {WORK}'",
    f"{SSH} 'mkdir -p {WORK}/runs/a {TMP}/b'",
    # sbatch: path options, a dependency, an array, the script's absolute-path arguments.
    sbatch(
        f"--job-name={PREFIX}a --partition={SHORT} --time=60 --array=1-40%20"
        f" --output={WORK}/logs/%x-%A_%a.out --error={WORK}/logs/%x-%A_%a.err"
        f" --chdir={WORK}/runs --dependency=afterok:123:456",
        f"{WORK}/kit/a.sbatch {WORK}/kit/in.tsv {TMP}/out label 5",
    ),
    sbatch(f"-J {PREFIX}t -p {SHORT} -t 179 -o {TMP}/t.out"),
    # Copies up and down.
    f"{RSYNC} {LOCAL}/kit/ {HOST}:{WORK}/kit/",
    f"{RSYNC} {HOST}:{WORK}/runs/x/summary.tsv {LOCAL}/out/",
    f"{RSYNC} {HOST}:{TMP}/logs/ {LOCAL}/logs/",
    f"scp {LOCAL}/a.txt {HOST}:{TMP}/a.txt",
    f"scp {HOST}:{TMP}/a.txt {LOCAL}/a.txt",
]

STILL_ASK = [
    # In the read-only folders only: reads, a cd, a copy back, a job's input.
    f"{SSH} 'ls {TEST}'",
    f"{SSH} 'cat {TEST}/kit/README.md'",
    f"{SSH} 'cd {TEST} && ls'",
    f"{SSH} 'cd {TEST} && cat {WORK}/a'",
    f"{RSYNC} {HOST}:{TEST}/kit/README.md {LOCAL}/",
    sbatch(f"--job-name={PREFIX}a --partition={SHORT} --time=30", f"{WORK}/a.sbatch {TEST}/x.bam"),
    sbatch(f"--job-name={PREFIX}a --partition={SHORT} --time=30 --export=ALL,REF={TEST}/ref.fa"),
    # Snapshot folders: in the read-only folders, and inside a work folder.
    f"{SSH} 'ls {TEST}/.snapshots/daily/project'",
    f"{SSH} 'cat {WORK}/.snapshots/daily/a.tsv'",
    f"{RSYNC} {HOST}:{WORK}/.snapshot/a {LOCAL}/",
    # A mixed command: one path outside the work folders is enough.
    f"{SSH} 'ls {WORK} {TEST}/kit'",
    f"{SSH} 'ls {WORK} && cat {TEST}/x'",
    f"{SSH} 'grep -f {TEST}/patterns {WORK}/a'",
    f"{SSH} 'ls {WORK}' && {SSH} 'ls {TEST}'",
    f"{RSYNC} {HOST}:{WORK}/a {HOST}:{TEST}/b {LOCAL}/out/",
    # A copy that stays inside HPC: between work folders, from a snapshot, one source of several.
    f"scp {HOST}:{WORK}/a {HOST}:{TMP}/a",
    f"scp {HOST}:{TEST}/.snapshots/d/a {HOST}:{WORK}/a",
    f"scp {LOCAL}/a {HOST}:{WORK}/b {HOST}:{TMP}/c",
    # scancel always asks, alone or in a chain.
    f"{SSH} 'scancel --name={PREFIX}pilot'",
    f"{SSH} 'ls {WORK} && scancel --name={PREFIX}pilot'",
    # What the guard cannot place: a relative sbatch path option, files listed in a file,
    # what a substitution prints, a local sed or sort (they can run a program).
    sbatch(f"--job-name={PREFIX}a --partition={SHORT} --time=30 --output=slurm.out"),
    sbatch(f"--job-name={PREFIX}a --partition={SHORT} --time=30 --out=slurm.out"),
    sbatch(f"--job-name={PREFIX}a --partition={SHORT} --time=30 -D runs"),
    f"{SSH} 'cd {WORK}/kit && sha256sum -c SHA256SUMS'",
    f"{SSH} 'sort --files0-from={WORK}/list'",
    f"{SSH} 'file -f {WORK}/list'",
    f"echo \"$({SSH} 'ls {WORK}')\"",
    f"{SSH} 'cat {WORK}/a.tsv' | sed -n 1,5p",
    f"{SSH} 'cat {WORK}/a.tsv' | sort",
]


# Patient data inside the work folders, as the operator's ask_patterns might mark it.
PATIENT = [
    r"/inputs/patient(/|$)",
    r"/inputs/prod(/|$)",
    r"/\.private(/|$)",
    r"\.(maf|bam|bai|cram|crai|sam|vcf|bcf|fastq|fq)(\.gz)?$",
]
NEW_KEYS = ("allow_in_remote_roots", "ask_patterns")


def settings_env(tmp_path_factory, name, settings):
    home = tmp_path_factory.mktemp(name)
    config = home / "guard.json"
    config.write_text(json.dumps(settings))
    return {**os.environ, "HOME": str(home), "HPC_GUARD_CONFIG": str(config)}


@pytest.fixture(scope="module")
def absent_env(tmp_path_factory):
    """Settings without either key: the hook as it was before them."""
    settings = {k: v for k, v in CONFIG.items() if k not in NEW_KEYS}
    return settings_env(tmp_path_factory, "absent", settings)


@pytest.fixture(scope="module")
def off_env(tmp_path_factory):
    """allow_in_remote_roots off, with ask patterns set: the patterns must change nothing."""
    return settings_env(tmp_path_factory, "off", {**CONFIG, "ask_patterns": PATIENT})


@pytest.fixture(scope="module")
def on_env(tmp_path_factory):
    return settings_env(tmp_path_factory, "on", {**CONFIG, "allow_in_remote_roots": True})


@pytest.fixture(scope="module")
def phi_env(tmp_path_factory):
    settings = {**CONFIG, "allow_in_remote_roots": True, "ask_patterns": PATIENT}
    return settings_env(tmp_path_factory, "phi", settings)


@pytest.mark.parametrize("command", PASS + ASK + DENY)
def test_the_setting_off_changes_no_output(command, absent_env, env, off_env):
    """Off (false, with or without ask patterns) gives the same bytes and exit code as no
    setting at all, for every case."""
    done = [run(HOOK, command, e) for e in (absent_env, env, off_env)]
    output = [(d.returncode, d.stdout, d.stderr) for d in done]
    assert output[0] == output[1] == output[2]


@pytest.mark.parametrize("command", IN_WORK)
def test_without_the_setting_work_folder_commands_still_ask(command, absent_env, env, off_env):
    assert [decision(command, e) for e in (absent_env, env, off_env)] == ["ask"] * 3


@pytest.mark.parametrize("command", IN_WORK)
def test_the_setting_allows_what_stays_in_the_work_folders(command, on_env):
    assert decision(command, on_env) == "allow"


@pytest.mark.parametrize("command", STILL_ASK)
def test_the_setting_still_asks_outside_the_work_folders(command, env, on_env):
    assert decision(command, env) == "ask"
    assert decision(command, on_env) == "ask"


@pytest.mark.parametrize("command", DENY)
def test_the_setting_never_allows_a_denied_command(command, on_env):
    assert decision(command, on_env) == "deny"


@pytest.mark.parametrize("command", PASS)
def test_the_setting_leaves_status_and_unrelated_commands_alone(command, on_env):
    assert decision(command, on_env) == "pass"


@pytest.mark.parametrize("command", ASK)
def test_the_setting_turns_an_ask_only_into_allow(command, on_env):
    assert decision(command, on_env) in ("ask", "allow")


def test_an_allow_uses_the_hook_output_format(on_env):
    done = run(HOOK, f"{SSH} 'ls {WORK}'", on_env)
    assert done.returncode == 0 and not done.stderr
    assert json.loads(done.stdout) == {
        "hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": "allow",
            "permissionDecisionReason": "HPC rules: inside the work folders (allow_in_remote_roots)",
        }
    }


@pytest.mark.parametrize("value", ["true", 1, 0, None, [], {}])
def test_a_non_boolean_setting_denies_every_ssh(value, tmp_path):
    """allow_in_remote_roots is true or false; anything else is a settings problem."""
    config = tmp_path / "guard.json"
    config.write_text(json.dumps({**CONFIG, "allow_in_remote_roots": value}))
    env = {**os.environ, "HOME": str(tmp_path), "HPC_GUARD_CONFIG": str(config)}
    done = run(HOOK, f"{SSH} 'ls {WORK}'", env)
    reason = json.loads(done.stdout)["hookSpecificOutput"]["permissionDecisionReason"]
    assert decision(f"{SSH} 'ls {WORK}'", env) == "deny" and "allow_in_remote_roots" in reason
    assert decision("ls -la", env) == "pass"


# ask_patterns: with allow_in_remote_roots on, patient data in the work folders still asks.

JOB = f"--job-name={PREFIX}a --partition={SHORT} --time=60"

PATIENT_ASK = [
    # A patient folder or file named by path, by a cd, or relative to one.
    f"{SSH} 'ls {WORK}/inputs/patient'",
    f"{SSH} 'cd {WORK}/inputs/prod && ls'",
    f"{SSH} 'cat {WORK}/inputs/prod/.private/ids.tsv'",
    f"{SSH} 'ls -la {WORK}/kit/.private'",
    f"{SSH} 'head -n 5 {WORK}/runs/p1/sample-a.maf'",
    f"{SSH} 'cd {WORK}/runs/p1 && zcat calls.vcf.gz'",
    f"{SSH} 'stat {WORK}/runs/p1/A.BAM'",  # matched without regard to case
    f"{SSH} 'mkdir -p {WORK}/inputs/patient/new'",
    # Copies: a patient file as the source, a patient folder as the destination.
    f"{RSYNC} {HOST}:{WORK}/runs/p1/sample-a.maf {LOCAL}/out/",
    f"{RSYNC} {LOCAL}/kit/ {HOST}:{WORK}/inputs/patient/",
    # A job that reads patient data is a new kind of job: by path, by name, in an option.
    sbatch(JOB, f"{WORK}/kit/a.sbatch {WORK}/inputs/prod/a.bam"),
    sbatch(JOB, f"{WORK}/kit/a.sbatch a.bam"),
    sbatch(f"{JOB} --export=ALL,BAM={WORK}/inputs/patient/a.cram"),
    # Once patterns are set: a glob, or a read that goes below the folders it names.
    f"{SSH} 'ls {WORK}/inputs/*/'",
    f"{SSH} 'tail -n 5 {WORK}/runs/_slurm/*.out'",
    f"{SSH} \"find {WORK}/runs -name '*.tsv'\"",
    f"{SSH} 'grep -r ALT {WORK}/runs'",
    f"{SSH} 'du -h {WORK}'",
    f"{SSH} 'ls -R {WORK}/runs'",
]

PATIENT_ALLOW = [
    f"{SSH} 'tail -n 50 {WORK}/runs/_slurm/{PREFIX}a-123.out'",  # a job log
    f"{RSYNC} {HOST}:{WORK}/runs/p1/summary.tsv {LOCAL}/out/",  # an aggregate copied back
    f"{RSYNC} {LOCAL}/kit/ {HOST}:{WORK}/kit/",
    f"{SSH} 'ls {WORK}/inputs'",  # the folder's names, not what the patient folders hold
    f"{SSH} 'du -sh {WORK}/runs'",
    f"{SSH} 'mkdir -p {WORK}/runs/_slurm'",
    sbatch(f"{JOB} --output={WORK}/runs/_slurm/%x-%j.out", f"{WORK}/kit/a.sbatch {WORK}/runs/p1 x"),
]


@pytest.mark.parametrize("command", PATIENT_ASK)
def test_ask_patterns_keep_patient_data_asking(command, on_env, phi_env):
    """Allowed by allow_in_remote_roots alone; the patterns are what make it ask."""
    assert decision(command, on_env) == "allow"
    assert decision(command, phi_env) == "ask"


@pytest.mark.parametrize("command", PATIENT_ALLOW)
def test_ask_patterns_still_allow_the_rest_of_the_work_folders(command, phi_env):
    assert decision(command, phi_env) == "allow"


@pytest.mark.parametrize(
    "command,expected", [(c, "pass") for c in PASS] + [(c, "deny") for c in DENY]
)
def test_ask_patterns_never_allow_or_deny_by_themselves(command, expected, phi_env):
    assert decision(command, phi_env) == expected


@pytest.mark.parametrize("value", ["/inputs/patient", [1], [""], ["("], None, {}])
def test_bad_ask_patterns_deny_every_ssh(value, tmp_path):
    """ask_patterns is a list of regular expressions; anything else is a settings problem."""
    config = tmp_path / "guard.json"
    config.write_text(json.dumps({**CONFIG, "allow_in_remote_roots": True, "ask_patterns": value}))
    env = {**os.environ, "HOME": str(tmp_path), "HPC_GUARD_CONFIG": str(config)}
    done = run(HOOK, f"{SSH} 'ls {WORK}'", env)
    reason = json.loads(done.stdout)["hookSpecificOutput"]["permissionDecisionReason"]
    assert decision(f"{SSH} 'ls {WORK}'", env) == "deny"
    assert "ask_patterns" in reason or "ask pattern" in reason
    assert decision("ls -la", env) == "pass"
