"""The PHI read guard hook (.claude/hooks/phi-read-guard.py) keeps patient data out of tool output.

Each case is a tool call the model might make and the decision the hook must give: "pass" (no
objection) or "deny". Most cases call the hook's `check` function in this process; a few run
the hook as the tool does, as a subprocess on made-up JSON. The PHI list is a made-up file,
with HOME pointed at a folder that does not exist, so the operator's own list
(~/.claude/phi-paths.txt) never reaches the tests. Every path and name below is invented. Only
the link cases touch the disk, and they build their own tree under a temporary folder.
"""

import ast
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

HOOK = Path(__file__).resolve().parents[1] / ".claude" / "hooks" / "phi-read-guard.py"
EXAMPLE = HOOK.parent / "phi-paths.example.txt"

HOME = "/srv/fakecohort/fakehome"
BASE = "/srv/fakecohort/harness"
OUT = f"{BASE}/hpc/out"  # copy-backs: <OUT>/<run label>/<arm>
PLAN = f"{BASE}/plan"  # local inputs
KIT = f"{BASE}/hpc/kit"
SAN = f"{KIT}/phi_sanitize.py"  # the sanitizer, named in the list
DEFAULT_SAN = f"{HOME}/test/gbcms/harness/hpc/kit/phi_sanitize.py"  # the sanitizer, by default
ELSEWHERE = "/srv/fakeproject"  # the working folder of a call that is not in a PHI folder
PATIENT_ARMS = ["panel", "across", "wes", "germ", "gatex"]
PUBLIC_ARMS = ["c246", "tr", "sub", "nocl", "delins", "ch_a1", "ch_a2", "ch_both"]
LIST = f"""\
# A made-up list, as the operator would write one.
{OUT}/*/panel
{OUT}/*/across            # a patient arm under any run label
{OUT}/*/wes
{OUT}/*/germ
{OUT}/*/gatex
{PLAN}/across
{PLAN}/out_*              # a pattern inside a folder name
~/fakehome-phi
sanitizer {SAN}
"""
P = f"{OUT}/lab1/panel"  # a PHI folder
F = f"{P}/x.maf"  # a made-up file in it
RUN = f"{ELSEWHERE}/trace.py"  # a made-up script
PIPE = f"2>&1 | python3 -I {SAN}"


def load_hook(listing, home=HOME):
    """The hook as a module, loaded with this list file, never the operator's."""
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("HOME", home)
        mp.setenv("PHI_PATHS_FILE", str(listing))
        spec = importlib.util.spec_from_file_location("phi_read_guard", HOOK)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    return module


def listing_in(tmp_path, text):
    path = tmp_path / "phi-paths.txt"
    path.write_text(text)
    return path


@pytest.fixture(scope="module")
def guard(tmp_path_factory):
    return load_hook(listing_in(tmp_path_factory.mktemp("list"), LIST))


def verdict(module, tool, tool_input, cwd=ELSEWHERE):
    return "pass" if module.check(tool, tool_input, cwd) is None else "deny"


def decision(command, guard, cwd=ELSEWHERE):
    return verdict(guard, "Bash", {"command": command}, cwd)


# --------------------------------------------------------------------------- Bash: allowed

PASS = [
    # (1) File operations that print no content.
    f"rsync -a claude-hpc:/work/runs/lab1/panel/ {P}/",
    f"rsync -a -e 'ssh -a -x' claude-hpc:/work/runs/lab1/panel/ {P}/",
    f'rsync -az -e "ssh -a -x" --exclude="*.tmp" claude-hpc:/work/runs/lab1/panel/ {P}/',
    f"rsync -a --partial --stats claude-hpc:/work/runs/lab1/panel/ {P}/",
    f"rsync -a {P}/ claude-hpc:/work/copy/panel/",  # to another host
    f"rsync -a {OUT}/lab1/wes/ {OUT}/lab2/wes/",  # into another PHI folder
    f"rsync -an --stats claude-hpc:/work/runs/lab1/panel/ {P}/ 2>&1 | grep -E 'Number of files|Total transferred'",
    f"rsync -a --stats claude-hpc:/work/runs/lab1/panel/ {P}/ 2>&1 | grep -Ei 'total transferred'",
    f"rsync -a --stats claude-hpc:/work/runs/lab1/panel/ {P}/ 2>&1 | tail -15",
    f"rsync -a --stats claude-hpc:/work/runs/lab1/panel/ {P}/ 2>&1 | tail -n 15",
    f"rsync -a --stats claude-hpc:/work/runs/lab1/panel/ {P}/ | tail -n 15",
    f"rsync -a {P}/ {OUT}/lab2/panel/ > /dev/null 2>&1",
    f"cp -R {P} {OUT}/lab2/panel",
    f"cp {RUN} {P}/trace.py",  # a script going in
    f"cp -p {F} {OUT}/lab2/panel/y.maf",
    f"cp -R {OUT}/lab1 {OUT}/lab2",  # a run folder to a run folder: the arms keep their names
    f"mv {OUT}/lab1 {OUT}/lab1_old",
    f"rsync -a {OUT}/lab1/ claude-hpc:/work/copy/lab1/",
    f"cp {ELSEWHERE}/a.txt {OUT}/lab1/notes.txt",
    f"mkdir -p {OUT}/lab5",
    f"mkdir {OUT}/lab5/panel",
    f"trash {OUT}/lab5",
    f"cp -- {F} {OUT}/lab2/panel/y.maf",
    f"cp -t {OUT}/lab2/panel {F}",
    f"cp {F}",  # cp itself refuses: no destination
    f"mv {P}/old.maf {P}/new.maf",
    f"mv {F} {OUT}/lab2/wes/",
    f"ln -s {P} {OUT}/lab2/panel",
    f"ln {F} {P}/y.maf",
    f"mkdir {OUT}/lab3/panel",
    f"mkdir -p {OUT}/lab3/panel/sub/deeper",
    f"mkdir -m 755 {OUT}/lab3/panel",
    f"trash {P}/old",
    f"trash -- {P}/old {OUT}/lab1/wes/older",
    f"cd {P} && cp x.maf y.maf",
    f"cd {OUT}/lab1 && mv panel/x.maf panel/y.maf",
    f"cd {OUT}/lab1 && mkdir panel/sub",
    # (2) A script whose whole output goes through the sanitizer, as the last stage.
    f"python3 -I {RUN} {P} {PIPE}",
    f"python3 -I {RUN} {F} --out {P}/rows.tsv {PIPE}",
    f"python3 -I {RUN} {F} {PIPE} --names-from {BASE}/runs.tsv",
    f"python3 -I {RUN} {F} {PIPE} --table",
    f"python3 -I {RUN} {F} {PIPE} --names-from {BASE}/runs.tsv --table",
    f"python3.12 -I {RUN} {F} {PIPE}",
    f"/opt/fakepython/bin/python -I {RUN} {F} {PIPE}",
    f"python3 -I -u {RUN} {F} {PIPE}",
    f"python3 -I -B {RUN} {F} {PIPE}",
    f"python3 -I {P}/probe.py {P} {PIPE}",  # a script written inside the PHI folder
    f"PYTHONHASHSEED=0 python3 -I {RUN} {F} {PIPE}",
    f"python3 -I {RUN} --in={F} {PIPE}",
    f"python3 -I {RUN} --in {F} {PIPE}",
    f"python3 -I {RUN} {F} 2>&1 | python3 -I {DEFAULT_SAN}",  # the default sanitizer path
    f"python3 -I {RUN} {F} 2>&1 | python3 -I ~/test/gbcms/harness/hpc/kit/phi_sanitize.py",
    f"python3 -I {RUN} {F} 2>&1 | python3 -I $HOME/test/gbcms/harness/hpc/kit/phi_sanitize.py",
    f"python3 -I {RUN} {F} 2>&1 | python3 -I {KIT}/../kit/phi_sanitize.py",  # the same file, spelt differently
    f"python3 -I {RUN} {F} 2>&1 | python3 -I '{SAN}'",
    f"cd {P} && python3 -I {RUN} x.maf {PIPE}",
    f"cd {KIT} && python3 -I {RUN} {F} 2>&1 | python3 -I phi_sanitize.py",  # a relative sanitizer path
    f"python3 -I {RUN} {F} 2>/dev/null | python3 -I {SAN}",  # stderr dropped: nothing printed
    f"python3 -I {RUN} {F} |& python3 -I {SAN}",
    f"echo start; python3 -I {RUN} {F} {PIPE}",  # a command may come before the pipeline
    f"cd {ELSEWHERE} && python3 -I {RUN} {F} {PIPE}",
    # (3) One command with stdout and stderr both in a file under a PHI root.
    f"python3 {RUN} {F} > {P}/out.txt 2>&1",
    f"python3 {RUN} {F} >{P}/out.txt 2>&1",
    f"python3 {RUN} {F} &> {P}/out.txt",
    f"python3 {RUN} {F} >& {P}/out.txt",
    f"python3 {RUN} {F} >> {P}/out.txt 2>&1",
    f"python3 {RUN} {F} 2>{P}/err.txt >{P}/out.txt",
    f"python3 {RUN} {F} > {P}/out.txt 2> {P}/out.txt",
    f"python3 {RUN} {F} > {P}/out.txt 2>/dev/null",
    f"python3 {RUN} {F} > /dev/null 2>&1",
    f"python3 -c 'import sys' > {P}/out.txt 2>&1",
    f"grep -c ALT {F} > {P}/n.txt 2>&1",
    f"echo done > {P}/marker.txt 2>&1",
    f"cd {P} && python3 {RUN} x.maf > out.txt 2>&1",  # a relative output name
    f"cd {OUT}/lab1 && python3 {RUN} panel/x.maf > panel/out.txt 2>&1",
    f"PHI_RUN=1 python3 {RUN} {F} > {P}/out.txt 2>&1",
    f"python3 {RUN} {F} > '{P}/out.txt' 2>&1",
    f"python3 {RUN} {F} > {OUT}/lab1/panel/../panel/out.txt 2>&1",
    # A command that names no PHI path.
    "ls -la",
    "git status",
    "cat README.md",
    "grep -rn foo src/ tests/",
    f"python3 -I {RUN} 2>&1 | head -5",
    f"python3 -I {RUN} {ELSEWHERE}/a.maf 2>&1 | python3 -I {SAN} | head",
    f"sed -n 1,5p {ELSEWHERE}/notes.md",
    f"cat {KIT}/phi_sanitize.py",  # the sanitizer is a script, not patient data
    f"find {KIT} -name '*.py'",
    f"ls {OUT}",  # a folder above the arms: run labels only
    f"ls {OUT}/lab1",  # the arm names, no content
    f"ls {OUT}/lab1 {OUT}/lab2",
    f"cat {OUT}/notes.md",
    f"cat {PLAN}/across_notes.md",  # a name that only starts with a PHI folder's name
    f"cat {PLAN}/across.md",
    f"cat {PLAN}/across-old/x",
    f"cat {OUT}/lab1/panel2/x",  # an arm that is not on the list
    f"cat {OUT}/lab1/panelx",
    "git commit -m 'docs: the guard' && git push origin chore/x",
    "echo hello",
    "echo 'it is done'; printf '%s\\n' ok",
    "true && cd /srv/fakeproject && ls",
    f"cd {OUT}/lab1 && ls",
    f"cd {OUT}/lab1/c246 && cat x.maf",
    f"cd {OUT} && ls lab1",
    f"sleep 1; date; cd {PLAN} && pwd",
    "cd - && ls",
    f"python3 -I {RUN} 2>&1 | python3 -I {SAN}",
    "unbalanced 'quote",  # unreadable, but it names no PHI path
    "cat <<",
    "",
    "   ",
    # A search that stops at or above the PHI folder lists no patient content.
    f"find {OUT} -maxdepth 1",
    f"find {OUT} -maxdepth 2 -type d",
    f"find {OUT}/lab1 -maxdepth 1",
    f"find {BASE} -maxdepth 2",
    f"find {BASE} -maxdepth 1 -name 'src_*' -type d",
    f"find {OUT} -name '*' -maxdepth 1",
    f"cd {OUT} && find . -maxdepth 1",
    f"cd {OUT} && find -maxdepth 1",
    f"du -sh {OUT}",
    f"du -s {BASE}",
    f"du -sh {OUT}/lab1",
    f"du -d 1 {OUT}",
    f"du --max-depth=2 {OUT}",
    f"du --summarize {BASE}",
    f"cd {OUT}/lab1 && du -sh",
    f"cd {OUT} && du -sh *",
    f"tree -L 2 {OUT}",
    f"tree -L 1 {BASE}",
]
# The public (HG002) arms are not on the list, under any run label.
PASS += [f"cat {OUT}/lab1/{arm}/x.maf" for arm in PUBLIC_ARMS]
PASS += [f"head -5 {OUT}/lab9/{arm}/out/a.maf" for arm in PUBLIC_ARMS]
PASS += [f"ls -R {OUT}/lab1/{arm}" for arm in PUBLIC_ARMS]
PASS += [f"find {OUT}/lab1/{arm} -name '*.maf'" for arm in PUBLIC_ARMS]
PASS += [f"grep -r ALT {OUT}/lab1/{arm}" for arm in PUBLIC_ARMS]
PASS += [f"cd {OUT}/lab1/{arm} && cat x.maf" for arm in PUBLIC_ARMS]


@pytest.mark.parametrize("command", PASS)
def test_a_command_the_rules_allow_passes(command, guard):
    assert decision(command, guard) == "pass", command


# --------------------------------------------------------------------------- Bash: denied

DENY = [
    # Reading tools on a PHI file or folder.
    f"cat {F}",
    f"head -5 {F}",
    f"tail -n 5 {F}",
    f"tail -f {F}",
    f"less {F}",
    f"more {F}",
    f"sed -n 1,5p {F}",
    f"awk '{{print $1}}' {F}",
    f"grep ALT {F}",
    f"grep -c ALT {F}",
    f"egrep ALT {F}",
    f"rg ALT {F}",
    f"ls {P}",
    f"ls -la {P}",
    f"ls {P}/",
    f"ls -R {P}",
    f"find {P} -name '*.bam'",
    f"find {P} -type f",
    f"find {P} -exec cat {{}} +",
    f"du -sh {P}",
    f"wc -l {F}",
    f"stat {F}",
    f"file {F}",
    f"zcat {F}.gz",
    f"gzip -dc {F}.gz",
    f"xxd {F}",
    f"strings {F}",
    f"samtools view {P}/FAKE-0001.bam",
    f"samtools view -c {P}/FAKE-0001.bam",
    f"tree {P}",
    f"sort {F}",
    f"cut -f1 {F}",
    f"diff {F} {F}",
    f"jq . {P}/a.json",
    f"touch {P}/x",
    f"chmod +x {P}/run.py",
    f"rm {F}",
    f"open {F}",
    f"pbcopy < {F}",
    f"wc -l < {F}",
    f"cat < {F}",
    f"cat <<< {F}",
    f"echo {F}",  # even a path that is only typed: the rules name no echo form
    f"echo 'copied to {P}'",
    f"printf '%s\\n' {P}",
    f"echo '$(cat {F})'",
    # A script without the sanitizer, or with something wrong around it.
    f"python3 -c \"print(open('{F}').read())\"",
    f"python3 -c 'import os; print(os.listdir(\"{P}\"))'",
    f"python3 {RUN} {F}",
    f"python3 -I {RUN} {F}",
    f"python3 -I {RUN} {F} 2>&1",
    f"python3 -I {RUN} {F} 2>&1 | head",
    f"python3 -I {RUN} {F} 2>&1 | tee {P}/out.txt",
    f"python3 -I {RUN} {F} 2>&1 | cat",
    f"python3 -I {RUN} {F} | python3 -I {SAN}",  # stderr would reach the screen
    f"python3 -I {RUN} {F} 2>err.txt | python3 -I {SAN}",
    f"python3 -I {RUN} {F} 2>&1 >{ELSEWHERE}/o.txt | python3 -I {SAN}",  # stdout goes to a file instead
    f"python3 {RUN} {F} 2>&1 | python3 -I {SAN}",  # the script runs without -I
    f"python3 -I {RUN} {F} 2>&1 | python3 {SAN}",  # the sanitizer runs without -I
    f"python3 -I {RUN} {F} 2>&1 | python3 -I {KIT}/other_sanitize.py",  # not the sanitizer
    f"python3 -I {RUN} {F} 2>&1 | python3 -I {ELSEWHERE}/phi_sanitize.py",  # a copy elsewhere
    f"python3 -I {RUN} {F} 2>&1 | python3 -I phi_sanitize.py",  # relative, from a folder that has none
    f"python3 -I {RUN} {F} 2>&1 | python3 -I $SANITIZER",
    f"python3 -I {RUN} {F} 2>&1 | python3 -I {SAN} | head",  # a second pipe after the sanitizer
    f"python3 -I {RUN} {F} 2>&1 | python3 -I {SAN} | tee {ELSEWHERE}/o.txt",
    f"python3 -I {RUN} {F} 2>&1 | python3 -I {SAN} 2>&1 | cat",
    f"python3 -I {RUN} {F} 2>&1 | python3 -I {SAN} > {ELSEWHERE}/o.txt",  # a redirect after it
    f"python3 -I {RUN} {F} 2>&1 | python3 -I {SAN} 2> {ELSEWHERE}/e.txt",
    f"python3 -I {RUN} {F} 2>&1 | python3 -I {SAN}; cat {F}",  # a second command after it
    f"python3 -I {RUN} {F} 2>&1 | python3 -I {SAN}; echo done",
    f"python3 -I {RUN} {F} 2>&1 | python3 -I {SAN} && echo done",
    f"python3 -I {RUN} {F} 2>&1 | python3 -I {SAN} || echo failed",
    f"python3 -I {RUN} {F} 2>&1 | python3 -I {SAN} &",
    f"python3 -I {RUN} {F} 2>&1 | python3 -I {SAN} --names-from $(cat {F})",
    f"python3 -I $(echo {RUN}) {F} {PIPE}",
    f"python3 -I {RUN} $(cat {F}) {PIPE}",
    f"python3 -I -c 'print(1)' {F} {PIPE}",
    f"python3 -I -m trace {F} {PIPE}",
    f"python3 -I -W error {RUN} {F} {PIPE}",
    f"python3 -I {P}/run.sh {F} {PIPE}",  # a script that is not Python
    f"bash {RUN} {F} {PIPE}",
    f"sh -c 'cat {F}' {PIPE}",
    f"env python3 -I {RUN} {F} {PIPE}",
    f"time python3 -I {RUN} {F} {PIPE}",
    f"nohup python3 -I {RUN} {F} {PIPE}",
    f"cat {F} | python3 -I {SAN}",  # a reading tool is not a script
    f"head -5 {F} 2>&1 | python3 -I {SAN}",
    f"tail -n 5 {F} | python3 -I {SAN}",
    f"{ELSEWHERE}/trace.sh {F} {PIPE}",
    # A redirect that does not keep every stream inside the PHI folders.
    f"python3 {RUN} {F} > {P}/out.txt",  # stderr reaches the screen
    f"python3 {RUN} {F} 2> {P}/err.txt",  # stdout reaches the screen
    f"python3 {RUN} {F} > {ELSEWHERE}/out.txt 2>&1",  # a file outside the PHI folders
    f"python3 {RUN} {F} > {P}/out.txt 2> {ELSEWHERE}/err.txt",
    f"python3 {RUN} {F} 2>&1 > {P}/out.txt",  # the order sends stderr to the screen
    f"python3 {RUN} {F} > {P}/out.txt 2>&1 | head",
    f"python3 {RUN} {F} > {P}/out.txt 2>&1 | cat",
    f"python3 {RUN} {F} 2>&1 | tee {P}/out.txt",
    f"python3 {RUN} {F} > {P}/out.txt 2>&1 && cat {P}/out.txt",  # reading the result back
    f"python3 {RUN} {F} > {P}/out.txt 2>&1; head {P}/out.txt",
    f"python3 {RUN} {F} > {P}/out.txt 2>&1 && head -3 {P}/out.txt",
    f"python3 {RUN} {F} > $OUT 2>&1",
    f"python3 {RUN} {F} > $(mktemp) 2>&1",
    f"python3 {RUN} {F} > {P}/out.txt 2>&1 &",
    f"python3 {RUN} $(cat {F}) > {P}/out.txt 2>&1",
    f"python3 {RUN} {F} > /dev/stdout 2>&1",
    f"python3 {RUN} {F} > /dev/tty 2>&1",
    f"python3 {RUN} {F} > '{ELSEWHERE}/out.txt' 2>&1",
    f"python3 {RUN} {F} >&2",
    f"echo done > {P}/marker.txt",
    # Copies that print names, or that put patient data where nothing guards it.
    f"rsync -av claude-hpc:/work/runs/lab1/panel/ {P}/",
    f"rsync -a -v claude-hpc:/work/runs/lab1/panel/ {P}/",
    f"rsync -avP claude-hpc:/work/runs/lab1/panel/ {P}/",
    f"rsync -a --progress claude-hpc:/work/runs/lab1/panel/ {P}/",
    f"rsync -a --verbose claude-hpc:/work/runs/lab1/panel/ {P}/",
    f"rsync -ai claude-hpc:/work/runs/lab1/panel/ {P}/",
    f"rsync -a --itemize-changes claude-hpc:/work/runs/lab1/panel/ {P}/",
    f"rsync -a --out-format=%n claude-hpc:/work/runs/lab1/panel/ {P}/",
    f"rsync -a --info=progress2 claude-hpc:/work/runs/lab1/panel/ {P}/",
    f"rsync -a --info=name claude-hpc:/work/runs/lab1/panel/ {P}/",
    f"rsync -a --list-only {P}/",
    f"rsync --list-only claude-hpc:/work/runs/lab1/panel/ {P}/",
    f"rsync -a {P}/",  # one path lists it
    f"rsync -a {P}/ {ELSEWHERE}/copy/",  # patient data would leave the PHI folders
    f"rsync -a {OUT}/lab1/panel/ ~/Desktop/panel/",
    f"rsync -a --delete claude-hpc:/work/runs/lab1/panel/ {P}/",
    f"rsync -a --rsync-path='cat /etc/passwd' claude-hpc:/work/runs/lab1/panel/ {P}/",
    f"rsync -a -e 'bash -c cat' claude-hpc:/work/runs/lab1/panel/ {P}/",
    f"rsync -a --rsh='sh' claude-hpc:/work/runs/lab1/panel/ {P}/",
    f"rsync -a claude-hpc:/work/runs/lab1/panel/ {P}/ > {ELSEWHERE}/log.txt",
    f"rsync -a claude-hpc:/work/runs/lab1/panel/ {P}/ 2>{ELSEWHERE}/log.txt",
    f"rsync -a --stats claude-hpc:/work/runs/lab1/panel/ {P}/ 2>&1 | tail -200",  # more than a summary
    f"rsync -a --stats claude-hpc:/work/runs/lab1/panel/ {P}/ 2>&1 | tail -n 200",
    f"rsync -a --stats claude-hpc:/work/runs/lab1/panel/ {P}/ 2>&1 | tail {P}/x",
    f"rsync -a --stats claude-hpc:/work/runs/lab1/panel/ {P}/ 2>&1 | head",
    f"rsync -a --stats claude-hpc:/work/runs/lab1/panel/ {P}/ 2>&1 | cat",
    f"rsync -a --stats claude-hpc:/work/runs/lab1/panel/ {P}/ 2>&1 | grep ALT",
    f"rsync -a --stats claude-hpc:/work/runs/lab1/panel/ {P}/ 2>&1 | grep -E 'Number of files|bam'",
    f"rsync -a --stats claude-hpc:/work/runs/lab1/panel/ {P}/ 2>&1 | grep -E 'Number of files' {P}/x",
    f"rsync -a --stats claude-hpc:/work/runs/lab1/panel/ {P}/ 2>&1 | grep -E 'Number of files|Total transferred' | head",
    f"rsync -a --stats claude-hpc:/work/runs/lab1/panel/ {P}/ 2>&1 | tail -15 | head",
    f"rsync -a --stats claude-hpc:/work/runs/lab1/panel/ {P}/ 2>&1 | tail -15 > {ELSEWHERE}/o.txt",
    f"rsync -a --stats claude-hpc:/work/runs/lab1/panel/ {P}/ 2>&1 | sed p",
    f"rsync -a --stats claude-hpc:/work/runs/lab1/panel/ {P}/ 2>&1 | tail -15; cat {F}",
    f"rsync -a --stats claude-hpc:/work/runs/lab1/panel/ {P}/ 2>&1 | grep -E 'Number of files' && cat {F}",
    f"rsync -a claude-hpc:/work/runs/lab1/panel/ {P}/ $(cat {F})",
    f"rsync -a $SRC {P}/",
    f"cp -v {F} {P}/y.maf",
    f"cp -iv {F} {P}/y.maf",
    f"cp --verbose {F} {P}/y.maf",
    f"cp --interactive {F} {P}/y.maf",
    f"cp {F} {ELSEWHERE}/y.maf",  # out of the PHI folders
    f"cp -R {P} /tmp/panel_copy",
    f"cp -t {ELSEWHERE} {F}",
    f"cp --target-directory={ELSEWHERE} {F}",
    f"cp {F} ~/Desktop",
    f"cp {F} {BASE}/copy.maf",  # beside the PHI folders, not inside one
    f"cp {F} {OUT}/lab2/notes.maf",  # a run folder that holds PHI folders is no PHI folder
    f"cp -R {OUT}/lab1 {ELSEWHERE}/lab1_copy",  # a run folder holds the patient arms
    f"cp -R {OUT}/lab1 /tmp/lab1_copy",
    f"mv {OUT}/lab1 {ELSEWHERE}/lab1_moved",
    f"mv {OUT}/lab1 ~/lab1_moved",
    f"rsync -a {OUT}/lab1/ {ELSEWHERE}/lab1_copy/",
    f"rsync -a {OUT}/ /tmp/out_copy/",
    f"cp -R {BASE} {ELSEWHERE}/harness_copy",
    f"ln -s {OUT}/lab1 {ELSEWHERE}/lab1_link",
    f"cp {F} /dev/stdout",
    f"cp {F} /dev/tty",
    f"cp {F} c:y.maf",  # a name with a colon is a local file for cp, not another host
    f"mv -v {F} {P}/y.maf",
    f"mv {F} {ELSEWHERE}/y.maf",
    f"mv {F} ~/y.maf",
    f"mv -t {ELSEWHERE} {F}",
    f"ln -s {F} {ELSEWHERE}/link",
    f"ln -sv {F} {P}/y.maf",
    f"ln -s {F} /dev/stdout",
    f"mkdir -v {OUT}/lab3/panel",
    f"mkdir -pv {OUT}/lab3/panel",
    f"mkdir --verbose {OUT}/lab3/panel",
    f"trash -l {P}",
    f"trash -e {P}",
    f"trash -v {P}/old",
    f"trash {P}/old > {ELSEWHERE}/o.txt",
    f"cp {F} {P}/y.maf | cat",
    f"cp {F} {P}/y.maf > {ELSEWHERE}/o.txt",
    f"cp -R $(cat {ELSEWHERE}/list.txt) {P}",
    f"cp $SRC {P}/y.maf",
    f"cp {F} {P}/y.maf <<< x",
    # Commands run in a PHI folder, whatever they name.
    f"cd {P} && cat x.maf",
    f"cd {P} && ls",
    f"cd {P} && head x.maf",
    f"cd {P} && grep ALT x.maf",
    f"cd {P} && cat ./x.maf",
    f"cd {P} && cat ../panel/x.maf",
    f"cd {P} && ls .",
    f"cd {P} && ls ..",
    f"cd {P} && python3 -I {RUN} 2>&1 | head",
    f"cd {P} && python3 -I {RUN} x.maf",
    f"cd {P} && python3 {RUN} x.maf > {ELSEWHERE}/out.txt 2>&1",
    f"cd {P}; cat x.maf",
    f"cd {P}\ncat x.maf",
    f"cd {P} && cd .. && cat panel/x.maf",
    f"cd {OUT}/lab1 && cat panel/x.maf",
    f"cd {OUT}/lab1 && ls panel",
    f"cd {OUT}/lab1 && cat ./panel/x.maf",
    f"cd {OUT}/lab1 && grep -r ALT panel",
    f"cd {OUT}/lab1 && cat */x.maf",
    f"cd {OUT}/lab1 && cat p*/x.maf",
    f"cd {OUT} && cat lab1/panel/x.maf",
    f"cd {OUT} && cat */panel/x.maf",
    f"cd '{P}' && cat x.maf",
    f'cd "{P}" && head x.maf',
    f"cd {P}/.. && cat panel/x.maf",
    f"cd {OUT}/lab1/c246/../panel && cat x.maf",
    f"cd {PLAN}/across && ls",
    "cd ~/fakehome-phi && cat x",
    "cd $HOME/fakehome-phi && cat x",
    f"cd {BASE} && cat hpc/out/lab1/panel/x.maf",
    f"cd {ELSEWHERE} && cd {P} && cat x.maf",
    f"cd {ELSEWHERE}; cd {P}; cat x.maf",
    f"pushd {P} && cat x.maf",
    # Ways of writing the same path.
    f'cat "{F}"',
    f"cat '{F}'",
    f"cat {OUT}/lab1/pan''el/x.maf",
    f'cat {OUT}/lab1/pa"nel"/x.maf',
    f"cat {OUT}/lab1/pa\\nel/x.maf",
    f"cat {OUT}/lab1/./panel/x.maf",
    f"cat {OUT}/lab1//panel//x.maf",
    f"cat {OUT}/lab1/c246/../panel/x.maf",
    f"cat {OUT}/lab1/c246/../../lab1/panel/x.maf",
    f"cat {OUT}/lab1/PANEL/x.maf",  # the file system may ignore case
    f"cat {OUT}/Lab1/Panel/x.maf",
    f"cat {BASE}/HPC/out/lab1/panel/x.maf",
    f"cat {OUT}/lab1/panel/",
    f"cat {OUT}/lab1/panel",
    f"cat {OUT}/lab1/panel/.hidden",
    f"cat {OUT}/lab1/{{panel,c246}}/x.maf",  # a brace list that holds an arm
    f"cat {OUT}/lab{{1,2}}/panel/x.maf",
    f"cat {OUT}/lab1/p*/x.maf",  # patterns
    f"cat {OUT}/lab1/pan?l/x.maf",
    f"cat {OUT}/lab1/pan[e]l/x.maf",
    f"cat {OUT}/lab1/*/x.maf",
    f"cat {OUT}/*/panel/x.maf",
    f"cat {OUT}/*/*/x.maf",
    f"cat {OUT}/lab*/panel/*",
    f"cat {OUT}/lab1/*",
    f"cat {OUT}/*/*",
    f"cat {PLAN}/*/x",
    f"cat {PLAN}/ac*/x",
    f"cat {PLAN}/out_a/x",  # a pattern inside a folder name
    f"cat {PLAN}/out_label7/rows/x",
    f"cat {PLAN}/out_/x",
    f"cat {PLAN}/across/x",
    f"cat {PLAN}/across",
    f"cat {PLAN}/ACROSS/x",
    "cat ~/fakehome-phi",
    "cat ~/fakehome-phi/x",
    "cat $HOME/fakehome-phi/x",
    "cat ${HOME}/fakehome-phi/x",
    'cat "$HOME/fakehome-phi/x"',
    'cat "${HOME}/fakehome-phi/x"',
    f"cat '{HOME}/fakehome-phi/x'",
    f"cat {HOME}/fakehome-phi/x",
    f"cat {HOME}/./fakehome-phi/x",
    f"cat {HOME}/test/../fakehome-phi/x",
    "cat ~/test/../fakehome-phi/x",
    "cat ~/fakehome-phi/sub/deeper/x",
    "FOO=~/fakehome-phi/x cat",
    f"cat --number={F}",  # an option's value
    f"cat -n{F}",
    f"grep -f {F} {ELSEWHERE}/a.txt",
    f"grep --file={F} {ELSEWHERE}/a.txt",
    f"awk -f{F} {ELSEWHERE}/a.txt",
    f"diff {ELSEWHERE}/a.txt {F}",
    f"cat {ELSEWHERE}/a.txt {F}",
    f"cat {ELSEWHERE}/a.txt:{F}",
    f"cat {ELSEWHERE}/a.txt,{F}",
    f"{F}",  # the command itself is a PHI path
    f"{P}/run.sh",
    f"{P}/run.sh > {ELSEWHERE}/out.txt 2>&1",
    # A variable that holds a PHI path cannot be followed, so the assignment is refused.
    f"X={F} cat $X",
    f"X={F}",
    f"X={F}; cat $X",
    f"export X={F}",
    f"read -r X < {F}",
    # Chains, groups, substitutions and shells.
    f"echo a; cat {F}",
    f"echo a && cat {F}",
    f"false || cat {F}",
    f"echo a\ncat {F}",
    f"echo a &\ncat {F}",
    f"echo a | cat {F}",
    f"ls {ELSEWHERE} | xargs cat {F}",
    f"echo '{F}' | xargs cat",
    f"printf '%s\\n' {F} | xargs head",
    f"cat {F} &",
    f"cat {F} && echo ok",
    f"cat {F} ; echo ok",
    f"(cat {F})",
    f"( cd {P} ; ls )",
    f"{{ cat {F}; }}",
    f"if true; then cat {F}; fi",
    f"for f in {P}/*; do cat $f; done",
    f"for f in {P}/a {P}/b; do head -1 $f; done",
    f"while read l; do echo $l; done < {F}",
    f"echo $(cat {F})",
    f'echo "$(cat {F})"',
    f"echo `cat {F}`",
    f'echo "`head -1 {F}`"',
    f"echo $(echo {F})",
    f"echo $(ls {P})",
    f"cat $(echo {F})",
    f"echo $(echo $(cat {F}))",
    f"echo $(echo $(echo $(echo $(echo {F}))))",
    f"diff <(cat {F}) {ELSEWHERE}/a.txt",
    f"cat <(head -1 {F})",
    f"echo >(cat {F})",
    f"echo $(( $(wc -l < {F}) + 1 ))",
    f"bash -c 'cat {F}'",
    f'bash -c "cat {F}"',
    f"sh -c 'head {F}'",
    f"zsh -c 'cat {F}'",
    f"eval 'cat {F}'",
    f"eval cat {F}",
    f"source {F}",
    f". {F}",
    f"env cat {F}",
    f"env X=1 cat {F}",
    f"xargs cat < {F}",
    f"xargs -a {F} cat",
    f"time cat {F}",
    f"nohup cat {F}",
    f"sudo cat {F}",
    f"command cat {F}",
    f"builtin cd {P}",
    f"exec cat {F}",
    f"\\cat {F}",
    f"/bin/cat {F}",
    f"/usr/bin/head {F}",
    f"watch cat {F}",
    f"timeout 5 cat {F}",
    f"script -q /dev/null cat {F}",
    f"git commit -m 'docs: read {F}'",
    f"gh issue comment 1 --body 'see {F}'",
    f'sh -c "cd {P} && cat x.maf"',
    f"python3 -c \"import os; os.system('cat {F}')\"",
    # Text that the shell ignores can still hide a script.
    f"cat > {ELSEWHERE}/s.py <<'EOF'\nprint(open('{F}').read())\nEOF",
    f"cat > {ELSEWHERE}/s.py <<EOF\nopen('{F}')\nEOF",
    f"cat > {ELSEWHERE}/s.py <<-EOF\n\topen('{F}')\n\tEOF",
    f"python3 - <<'EOF'\nprint(open('{F}').read())\nEOF",
    f"python3 -I - <<'EOF' {PIPE}\nprint(open('{F}').read())\nEOF",
    f"echo hi # {F}",
    f"echo hi\n# {F}\necho ok",
    # A command that cannot be read and may name a PHI path fails closed.
    f"cat '{F}",
    f'cat "{F}',
    f"cat {F} `",
    f"cat $(cat {F}",
    f"cat {F} |",
    f"cat {F} >",
    f"echo {F} | | cat",
    f"cat {F} >&",
]
# Every patient arm, under any run label, in any command that reads.
DENY += [f"cat {OUT}/lab1/{arm}/x.maf" for arm in PATIENT_ARMS]
DENY += [f"cat {OUT}/lab9/{arm}/out/y.trace" for arm in PATIENT_ARMS]
DENY += [f"ls {OUT}/lab1/{arm}" for arm in PATIENT_ARMS]
DENY += [f"head -2 {OUT}/run_with_a_long-label.v2/{arm}/x" for arm in PATIENT_ARMS]
DENY += [f"cd {OUT}/lab1/{arm} && head x" for arm in PATIENT_ARMS]
DENY += [f"python3 -I {RUN} {OUT}/lab1/{arm} 2>&1 | head" for arm in PATIENT_ARMS]
# A folder that holds a PHI folder is guarded against the tools that search below it.
DENY += [
    f"find {OUT} -name '*.maf'",
    f"find {OUT}/lab1 -type f",
    f"find {BASE} -name x",
    f"find {OUT}/ -iname '*.bam'",
    f"find {OUT} {OUT}/lab1/c246 -name x",
    f"find {OUT}/lab1/c246 {OUT}/lab1 -name x",
    f"find -L {OUT} -name x",
    f"du -a {BASE}",
    f"tree {OUT}",
    f"tree -L 2 {OUT}/lab1",
    f"ls -R {OUT}",
    f"ls -lR {OUT}/lab1",
    f"ls --recursive {OUT}",
    f"grep -r ALT {OUT}",
    f"grep -rn ALT {OUT}/lab1",
    f"grep -R ALT {OUT}",
    f"grep --recursive ALT {OUT}",
    f"grep -rl ALT {OUT}/lab1 {OUT}/lab2",
    f"grep -r --include='*.maf' ALT {OUT}",
    f"rg ALT {OUT}",
    f"rg ALT {OUT}/lab1",
    f"rg -g '*.maf' ALT {BASE}",
    f"ag ALT {OUT}",
    f"fd maf {OUT}",
    f"diff -r {OUT}/lab1 {OUT}/lab2",
    f"diff -rq {OUT}/lab1 {ELSEWHERE}",
    f"cd {OUT}/lab1 && find . -name x",
    f"cd {OUT}/lab1 && rg ALT",
    f"cd {OUT}/lab1 && grep -r ALT .",
    f"cd {OUT}/lab1 && ls -R",
    f"cd {OUT} && find . -name x",
    f"cd {OUT} && grep -r ALT lab1",
    f"cd {OUT}/lab1 && du -a",
    f"cd {OUT}/lab1 && du -h",
    f"cd {OUT} && du -a *",
    # ... unless the search stops at or above the PHI folder: it then lists no content.
    f"find {OUT} -maxdepth 3",
    f"find {OUT} -maxdepth 9 -type f",
    f"find {OUT} -maxdepth x",
    f"find {OUT} -maxdepth",
    f"find {OUT} -mindepth 1",
    f"find {BASE} -maxdepth 3",  # the plan folder holds a PHI folder two levels down
    f"find {BASE} -maxdepth 5",
    f"find {OUT} -maxdepth 2 {OUT}/lab1 -maxdepth 4",
    f"find {OUT} -maxdepth 3 -exec cat {{}} +",
    f"du -a {OUT}",
    f"du -h {OUT}",
    f"du -d 3 {OUT}",
    f"du --max-depth=3 {OUT}",
    f"du -d x {OUT}",
    f"tree -L 3 {OUT}",
    f"tree -L x {OUT}",
    f"ls -R {OUT}",
]


@pytest.mark.parametrize("command", DENY)
def test_a_command_outside_the_rules_is_denied(command, guard):
    assert decision(command, guard) == "deny", command


def test_a_command_in_a_phi_folder_is_denied_whatever_it_names(guard):
    """A call that starts in a PHI folder (an earlier `cd`) reads there with any relative name."""
    cwd = f"{OUT}/lab1/panel"
    for command in (
        "ls",
        "cat x.maf",
        "head -2 ../panel/x.maf",
        "git status",
        "python3 -I s.py 2>&1 | head",
    ):
        assert decision(command, guard, cwd) == "deny", command
    # ... but the forms the rules allow stay allowed, and a move out is not blocked.
    assert decision(f"python3 -I {RUN} x.maf {PIPE}", guard, cwd) == "pass"
    assert decision("cp x.maf y.maf", guard, cwd) == "pass"
    assert decision("python3 run.py > out.txt 2>&1", guard, cwd) == "pass"
    assert decision("echo done", guard, cwd) == "pass"
    assert decision("cd /srv/fakeproject", guard, cwd) == "pass"
    assert decision("cd /srv/fakeproject && ls", guard, cwd) == "pass"
    assert decision("cd ~ && pwd", guard, cwd) == "pass"


def test_relative_names_are_resolved_with_the_call_folder(guard):
    cwd = f"{OUT}/lab1"
    assert decision("cat panel/x.maf", guard, cwd) == "deny"
    assert decision("cat ./panel/x.maf", guard, cwd) == "deny"
    assert decision("cat ../lab1/panel/x.maf", guard, cwd) == "deny"
    assert decision("cat c246/../panel/x.maf", guard, cwd) == "deny"
    assert decision("cat p*/x.maf", guard, cwd) == "deny"
    assert decision("cat c246/x.maf", guard, cwd) == "pass"
    assert decision("cat ../lab2/c246/x.maf", guard, cwd) == "pass"
    assert decision("ls", guard, cwd) == "pass"
    assert decision("ls panel", guard, cwd) == "deny"
    assert decision("find . -name x", guard, cwd) == "deny"
    assert decision("grep -r ALT .", guard, cwd) == "deny"
    assert decision("rg ALT", guard, cwd) == "deny"
    assert decision("du -a", guard, cwd) == "deny"
    assert decision("du -sh", guard, cwd) == "pass"
    assert decision("find . -maxdepth 1", guard, cwd) == "pass"
    assert decision("find . -maxdepth 2", guard, cwd) == "deny"
    assert decision("find c246 -name x", guard, cwd) == "pass"
    assert decision("cp panel/x.maf panel/y.maf", guard, cwd) == "pass"
    assert decision(f"python3 -I {RUN} panel {PIPE}", guard, cwd) == "pass"
    assert decision("python3 -I run.py panel 2>&1 | head", guard, cwd) == "deny"
    assert decision("python3 run.py panel > panel/out.txt 2>&1", guard, cwd) == "pass"
    assert (
        decision("python3 run.py panel > out.txt 2>&1", guard, cwd) == "deny"
    )  # lab1 is not a PHI folder
    # A relative name from an unrelated folder is not patient data.
    assert decision("cat panel/x.maf", guard, ELSEWHERE) == "pass"


def test_the_reason_says_what_to_do_instead(guard):
    reason = guard.check("Bash", {"command": f"cat {F}"}, ELSEWHERE)
    assert reason.startswith("`cat` prints patient file content")
    assert F not in reason  # the reason never repeats the path that was asked for
    how = guard.message(reason)
    assert "python3 -I ~/" not in how and SAN in how
    assert "phi_sanitize.py" in how and "> <PHI folder>/out.txt 2>&1" in how
    assert (
        "rsync, cp, mv, ln, mkdir, trash" in how and "security.md" in how and "phi-paths.txt" in how
    )
    assert "rsync option -av" in guard.check(
        "Bash", {"command": f"rsync -av claude-hpc:/w/ {P}/"}, ELSEWHERE
    )
    assert "outside the PHI folders" in guard.check(
        "Bash", {"command": f"cp {F} /tmp/x"}, ELSEWHERE
    )
    assert "substitution" in guard.check("Bash", {"command": f"echo $(cat {F})"}, ELSEWHERE)
    heredoc = f"cat > {ELSEWHERE}/s.py <<'EOF'\nopen('{F}')\nEOF"
    assert "Write tool" in guard.check("Bash", {"command": heredoc}, ELSEWHERE)


# --------------------------------------------------------------------------- the file tools


def test_read_and_notebook_read_may_not_touch_a_phi_path(guard):
    for tool, key in (("Read", "file_path"), ("NotebookRead", "notebook_path")):
        for path in (
            F,
            P,
            f"{OUT}/lab1/wes/a/b/c.txt",
            "~/fakehome-phi/x.ipynb",
            "$HOME/fakehome-phi/x",
        ):
            assert verdict(guard, tool, {key: path}) == "deny", (tool, path)
        assert verdict(guard, tool, {key: f"{OUT}/lab1/c246/x.maf"}) == "pass"
        assert verdict(guard, tool, {key: f"{KIT}/phi_sanitize.py"}) == "pass"
        assert verdict(guard, tool, {key: "relative/notes.md"}) == "pass"
        assert verdict(guard, tool, {key: f"{PLAN}/across_notes.md"}) == "pass"
        assert verdict(guard, tool, {key: "panel/x.maf"}, f"{OUT}/lab1") == "deny"
        assert verdict(guard, tool, {key: "x.maf"}, P) == "deny"
        assert verdict(guard, tool, {key: "../lab1/panel/x.maf"}, f"{OUT}/lab2") == "deny"
        assert verdict(guard, tool, {key: "x.maf"}, f"{OUT}/lab1/c246") == "pass"
        assert verdict(guard, tool, {key: f"{OUT}/lab1/c246/../panel/x.maf"}) == "deny"
        assert verdict(guard, tool, {key: f"{OUT}/lab1/PANEL/x.maf"}) == "deny"
    assert verdict(guard, "Read", {"path": F}) == "deny"
    assert verdict(guard, "Read", {"file_path": F, "offset": 1, "limit": 5}) == "deny"
    assert verdict(guard, "Read", {}) == "pass"


def test_grep_may_not_search_a_phi_path_or_a_folder_that_holds_one(guard):
    pattern = {"pattern": "ALT"}
    assert verdict(guard, "Grep", {**pattern, "path": P}) == "deny"
    assert verdict(guard, "Grep", {**pattern, "path": F}) == "deny"
    assert verdict(guard, "Grep", {**pattern, "path": "~/fakehome-phi"}) == "deny"
    assert verdict(guard, "Grep", {**pattern, "path": "panel"}, f"{OUT}/lab1") == "deny"
    assert verdict(guard, "Grep", {**pattern}, P) == "deny"  # no path: the call folder
    assert verdict(guard, "Grep", {**pattern}, f"{P}/sub") == "deny"
    # a folder that holds a PHI folder: the search reaches below it
    assert verdict(guard, "Grep", {**pattern, "path": f"{OUT}/lab1"}) == "deny"
    assert verdict(guard, "Grep", {**pattern, "path": OUT}) == "deny"
    assert verdict(guard, "Grep", {**pattern, "path": BASE}) == "deny"
    assert verdict(guard, "Grep", {**pattern}, f"{OUT}/lab1") == "deny"
    assert verdict(guard, "Grep", {**pattern}, OUT) == "deny"
    assert verdict(guard, "Grep", {**pattern, "path": ELSEWHERE, "glob": f"{P}/*.maf"}) == "deny"
    # folders that hold none
    assert verdict(guard, "Grep", {**pattern, "path": f"{OUT}/lab1/c246"}) == "pass"
    assert verdict(guard, "Grep", {**pattern, "path": "c246"}, f"{OUT}/lab1") == "pass"
    assert verdict(guard, "Grep", {**pattern, "path": ELSEWHERE}) == "pass"
    assert verdict(guard, "Grep", {**pattern, "path": KIT}) == "pass"  # beside it, not above it
    assert verdict(guard, "Grep", {**pattern}) == "pass"
    assert verdict(guard, "Grep", {**pattern, "path": ELSEWHERE, "glob": "*.maf"}) == "pass"
    assert (
        verdict(guard, "Grep", {"pattern": P, "path": ELSEWHERE}) == "pass"
    )  # the pattern is text


def test_glob_may_not_list_a_phi_path_and_a_pattern_says_where_it_looks(guard):
    star = {"pattern": "*.maf"}
    assert verdict(guard, "Glob", {**star, "path": P}) == "deny"
    assert verdict(guard, "Glob", {**star, "path": "~/fakehome-phi"}) == "deny"
    assert verdict(guard, "Glob", {**star, "path": "panel"}, f"{OUT}/lab1") == "deny"
    assert verdict(guard, "Glob", {**star}, P) == "deny"  # no path: the call folder
    assert verdict(guard, "Glob", {"pattern": f"{OUT}/*/panel/*.maf"}) == "deny"
    assert verdict(guard, "Glob", {"pattern": f"{OUT}/lab1/panel/**"}) == "deny"
    assert verdict(guard, "Glob", {"pattern": "panel/*.maf", "path": f"{OUT}/lab1"}) == "deny"
    assert verdict(guard, "Glob", {"pattern": "*/panel/*", "path": OUT}) == "deny"
    assert (
        verdict(guard, "Glob", {"pattern": "*", "path": f"{OUT}/lab1"}) == "deny"
    )  # lists the arms
    assert verdict(guard, "Glob", {"pattern": "**/*.bam", "path": OUT}) == "deny"
    assert verdict(guard, "Glob", {"pattern": "**/*.bam", "path": BASE}) == "deny"
    assert verdict(guard, "Glob", {"pattern": "**/*.bam"}, BASE) == "deny"
    assert verdict(guard, "Glob", {"pattern": "**/*.maf"}, f"{OUT}/lab1") == "deny"
    # a pattern that stays out of every PHI folder is fine whatever folder it starts in
    assert verdict(guard, "Glob", {**star, "path": f"{OUT}/lab1"}) == "pass"
    assert verdict(guard, "Glob", {**star, "path": OUT}) == "pass"
    assert verdict(guard, "Glob", {**star}, f"{OUT}/lab1") == "pass"
    assert verdict(guard, "Glob", {"pattern": f"{OUT}/lab1/c246/*.maf"}) == "pass"
    assert verdict(guard, "Glob", {"pattern": "c246/*.maf", "path": f"{OUT}/lab1"}) == "pass"
    assert verdict(guard, "Glob", {"pattern": "hpc/kit/*.py", "path": BASE}) == "pass"
    assert verdict(guard, "Glob", {"pattern": "**/*.py", "path": KIT}) == "pass"
    assert verdict(guard, "Glob", {"pattern": "src/**/*.py"}) == "pass"
    assert verdict(guard, "Glob", {**star, "path": ELSEWHERE}) == "pass"


def test_write_and_edit_of_a_script_in_a_phi_folder_are_allowed(guard):
    """Write and Edit show only what the model wrote, and refuse a file that was not read first."""
    for tool in ("Write", "Edit", "MultiEdit", "NotebookEdit"):
        assert verdict(guard, tool, {"file_path": f"{P}/probe.py", "content": "print(1)"}) == "pass"
        assert (
            verdict(guard, tool, {"file_path": f"{PLAN}/across/run.py", "new_string": "x"})
            == "pass"
        )


@pytest.mark.parametrize(
    "tool", ["WebFetch", "Task", "Agent", "TodoWrite", "mcp__x__y", "", None, "Skill"]
)
def test_other_tools_are_not_the_guards_business(tool, guard):
    assert verdict(guard, tool, {"command": f"cat {F}", "file_path": F, "path": F}) == "pass"


# ----------------------------------------------------------------- the list, and its failures


def test_a_fault_in_the_hook_fails_closed_only_for_a_phi_path(guard):
    # A command the tokenizer cannot read: closed when it names a PHI path, open when it does not.
    assert verdict(guard, "Bash", {"command": f"cat '{F}"}) == "deny"
    assert verdict(guard, "Bash", {"command": "cat 'unclosed"}) == "pass"
    assert verdict(guard, "Bash", {"command": "echo $(cat"}) == "pass"
    assert verdict(guard, "Bash", {"command": 123}) == "pass"
    assert verdict(guard, "Bash", {"command": None}) == "pass"
    assert verdict(guard, "Read", {"file_path": 123}) == "pass"
    assert verdict(guard, "Grep", {"pattern": "x", "path": None}) == "pass"


def test_a_long_command_is_still_checked(guard):
    long = "echo " + "x" * 200_000
    assert decision(long, guard) == "pass"
    assert decision(long + f"; cat {F}", guard) == "deny"
    assert decision(" ".join(["ls -la;"] * 2000), guard) == "pass"
    assert decision("echo $(" * 40 + f"cat {F}" + ")" * 40, guard) == "deny"


def test_an_empty_list_guards_nothing(tmp_path):
    module = load_hook(listing_in(tmp_path, "# nothing yet\n\n   \n"))
    assert verdict(module, "Bash", {"command": f"cat {F}"}) == "pass"
    assert verdict(module, "Read", {"file_path": F}) == "pass"


@pytest.mark.parametrize(
    "text",
    [
        "relative/path\n",
        f"{OUT}/*/panel\nrelative/path\n",
        "/\n",
        "//\n",
        "/*/panel\n",
        "/srv/**/panel\n",
        "sanitizer relative/phi_sanitize.py\n",
        "~relative/x\n",
    ],
)
def test_a_list_with_a_bad_line_denies_every_call_until_fixed(text, tmp_path):
    """A list that cannot guard anything must not look like one that does."""
    module = load_hook(listing_in(tmp_path, text))
    for tool, tool_input in (
        ("Bash", {"command": "ls"}),
        ("Bash", {"command": f"cat {F}"}),
        ("Read", {"file_path": "/srv/fakeproject/a.txt"}),
        ("Grep", {"pattern": "x"}),
        ("Glob", {"pattern": "*"}),
        ("NotebookRead", {"notebook_path": "a.ipynb"}),
    ):
        reason = module.check(tool, tool_input, ELSEWHERE)
        assert reason and "unusable" in reason and "operator" in reason, (tool, text)
    assert verdict(module, "Write", {"file_path": "/srv/fakeproject/a.txt"}) == "pass"


def test_a_named_list_that_is_missing_or_unreadable_denies_every_call(tmp_path):
    missing = load_hook(tmp_path / "none.txt")
    assert verdict(missing, "Bash", {"command": "ls"}) == "deny"
    assert verdict(missing, "Read", {"file_path": "/srv/fakeproject/a.txt"}) == "deny"
    (tmp_path / "folder").mkdir()
    unreadable = load_hook(tmp_path / "folder")  # a folder, not a file
    assert verdict(unreadable, "Bash", {"command": "ls"}) == "deny"


def test_list_comments_blank_lines_and_spacing_are_ignored(tmp_path):
    module = load_hook(
        listing_in(
            tmp_path, f"\n   # note\n\t{OUT}/*/panel   # trailing note\n\n  {PLAN}/across  \n"
        )
    )
    assert decision(f"cat {F}", module) == "deny"
    assert decision(f"cat {PLAN}/across/x", module) == "deny"
    assert decision(f"cat {OUT}/lab1/wes/x", module) == "pass"
    assert decision("cat '# note'", module) == "pass"


def test_an_entry_may_be_a_root_a_pattern_or_a_pattern_in_a_name(tmp_path):
    module = load_hook(listing_in(tmp_path, f"{PLAN}/across\n{OUT}/*/wes\n{BASE}/inputs/gate_*\n"))
    for path in (
        f"{PLAN}/across",
        f"{PLAN}/across/a/b/c",
        f"{OUT}/anylabel/wes/x",
        f"{OUT}/any.label-2/wes/x",
        f"{BASE}/inputs/gate_pi/x",
        f"{BASE}/inputs/gate_local/x",
    ):
        assert decision(f"cat {path}", module) == "deny", path
    for path in (
        f"{PLAN}/other/x",
        f"{PLAN}",
        f"{OUT}/anylabel/panel/x",
        f"{OUT}/wes/x",  # an arm folder with no run label above it
        f"{OUT}/a/b/wes/x",  # two levels: the pattern is one level
        f"{BASE}/inputs/gate/x",
        f"{BASE}/inputs/other_gate_pi/x",
    ):
        assert decision(f"cat {path}", module) == "pass", path


def test_a_second_sanitizer_may_be_named_in_the_list(tmp_path):
    other, other2 = f"{ELSEWHERE}/tools/mask.py", f"{ELSEWHERE}/tools/mask2.py"
    module = load_hook(
        listing_in(tmp_path, f"{OUT}/*/panel\nsanitizer {other}\nsanitizer = {other2}\n")
    )
    for san in (other, other2, DEFAULT_SAN):
        assert decision(f"python3 -I {RUN} {F} 2>&1 | python3 -I {san}", module) == "pass", san
    assert decision(f"python3 -I {RUN} {F} 2>&1 | python3 -I {SAN}", module) == "deny"
    assert (
        decision(f"python3 -I {RUN} {F} 2>&1 | python3 -I {ELSEWHERE}/tools/mask3.py", module)
        == "deny"
    )


def test_the_deny_reason_names_the_listed_sanitizer(tmp_path):
    module = load_hook(
        listing_in(tmp_path, f"{OUT}/*/panel\nsanitizer {ELSEWHERE}/tools/mask.py\n")
    )
    assert f"python3 -I {ELSEWHERE}/tools/mask.py" in module.HOW_TO.replace(
        "{san}", module.SANITIZERS[-1]
    )


# ------------------------------------------------------------------ symlinks and the real disk


def tree(tmp_path):
    """A made-up tree: real/<label>/{panel,c246}, links to an arm, a label and a public arm."""
    root = (tmp_path / "real").resolve()
    for arm in ("panel", "c246"):
        (root / "lab1" / arm).mkdir(parents=True)
        (root / "lab1" / arm / "x.maf").write_text("FAKE\n")
    links = tmp_path / "links"
    links.mkdir()
    (links / "to_arm").symlink_to(root / "lab1" / "panel")
    (links / "to_label").symlink_to(root / "lab1")
    (links / "to_public").symlink_to(root / "lab1" / "c246")
    (tmp_path / "linkroot").symlink_to(root)
    return root, links


def test_a_link_to_a_phi_folder_is_a_phi_path(tmp_path):
    root, links = tree(tmp_path)
    module = load_hook(listing_in(tmp_path, f"{root}/*/panel\n"))
    assert decision(f"cat {links}/to_arm/x.maf", module) == "deny"
    assert decision(f"ls {links}/to_arm", module) == "deny"
    assert decision(f"cat {links}/to_label/panel/x.maf", module) == "deny"
    assert decision(f"cd {links}/to_arm && cat x.maf", module) == "deny"
    assert decision(f"find {links}/to_label -name x.maf", module) == "deny"
    assert decision(f"cat {links}/to_public/x.maf", module) == "pass"
    assert decision("cat x.maf", module, str(links / "to_arm")) == "deny"
    assert decision("cat x.maf", module, str(links / "to_public")) == "pass"
    assert verdict(module, "Read", {"file_path": f"{links}/to_arm/x.maf"}) == "deny"
    assert verdict(module, "Read", {"file_path": f"{links}/to_public/x.maf"}) == "pass"
    assert verdict(module, "Grep", {"pattern": "F", "path": f"{links}/to_arm"}) == "deny"
    assert verdict(module, "Glob", {"pattern": "*", "path": f"{links}/to_arm"}) == "deny"
    # a link that makes a PHI folder reachable elsewhere is a copy that escapes the list
    assert decision(f"ln -s {root}/lab1/panel {tmp_path}/new_link", module) == "deny"


def test_a_list_root_reached_through_a_link_guards_the_real_folder_too(tmp_path):
    root, _ = tree(tmp_path)
    module = load_hook(
        listing_in(tmp_path, f"{tmp_path}/linkroot/*/panel\n")
    )  # the list spells the link
    assert decision(f"cat {tmp_path}/linkroot/lab1/panel/x.maf", module) == "deny"
    assert decision(f"cat {root}/lab1/panel/x.maf", module) == "deny"
    assert decision(f"cat {root}/lab1/c246/x.maf", module) == "pass"
    assert verdict(module, "Read", {"file_path": f"{root}/lab1/panel/x.maf"}) == "deny"


# ------------------------------------------------------------------------ the list's format


def test_parse_list_reads_entries_and_sanitizers(guard):
    entries, sanitizers = guard.parse_list(f"{OUT}/*/panel\n# c\n~/x\nsanitizer {SAN}\n")
    assert entries == [
        ["srv", "fakecohort", "harness", "hpc", "out", "*", "panel"],
        ["srv", "fakecohort", "fakehome", "x"],
    ]
    assert sanitizers[-1] == SAN and sanitizers[0] == DEFAULT_SAN


@pytest.mark.parametrize("text", ["rel/x\n", "/\n", "/a/**/b\n", "/*/x\n", "sanitizer rel.py\n"])
def test_parse_list_rejects_a_line_it_cannot_use(guard, text):
    with pytest.raises(ValueError):
        guard.parse_list(text)


def test_the_example_list_parses_and_guards(tmp_path):
    text = EXAMPLE.read_text()
    module = load_hook(listing_in(tmp_path, text))
    assert module.PROBLEM is None
    for arm in PATIENT_ARMS:
        assert any(e[-1] == arm for e in module.ENTRIES), f"the example lists the {arm} arm"
    for arm in PUBLIC_ARMS:
        assert not any(
            arm in e for e in module.ENTRIES
        ), f"the example must not list the public {arm} arm"
    lines = [ln.split("#")[0].strip() for ln in text.splitlines()]
    patient = next(ln for ln in lines if ln.endswith("/*/panel")).replace("*", "run1")
    public = patient.replace("/panel", "/c246")
    plain = next(ln for ln in lines if ln.endswith("/inputs/across"))
    assert decision(f"cat {patient}/x.maf", module) == "deny"
    assert decision(f"cat {plain}/x.maf", module) == "deny"
    assert decision(f"cat {public}/x.maf", module) == "pass"
    assert verdict(module, "Read", {"file_path": f"{patient}/x.maf"}) == "deny"
    assert decision("ls", module) == "pass"
    # the pattern inside a folder name, and the home folder
    assert decision(f"cat {plain.rsplit('/', 1)[0]}/out_label1/x", module) == "deny"
    assert decision("cat ~/example-project/patient-slices/x", module) == "deny"


def test_the_example_list_uses_made_up_roots_only():
    text = EXAMPLE.read_text().lower()
    for real in ("shahr2", "/users/", "msk"):
        assert real not in text, real


# ------------------------------------------------------------------ the hook as the tool runs it


def run(payload, tmp_path, listing=None, home=HOME):
    """The hook as a subprocess: JSON on stdin, the decision on stdout."""
    env = {k: v for k, v in os.environ.items() if k not in ("PHI_PATHS_FILE", "PYTHONPATH")}
    env["HOME"] = home
    if listing is not None:
        env["PHI_PATHS_FILE"] = str(listing)
    done = subprocess.run(
        [sys.executable, "-I", str(HOOK)],
        input=payload if isinstance(payload, str) else json.dumps(payload),
        capture_output=True,
        text=True,
        env=env,
    )
    assert done.returncode == 0, done.stderr
    assert not done.stderr.strip(), done.stderr
    return json.loads(done.stdout)


def test_the_hook_denies_with_the_json_the_tool_reads(tmp_path):
    listing = listing_in(tmp_path, LIST)
    call = {"tool_name": "Bash", "tool_input": {"command": f"cat {F}"}, "cwd": ELSEWHERE}
    out = run(call, tmp_path, listing)
    verdict_ = out["hookSpecificOutput"]
    assert verdict_["hookEventName"] == "PreToolUse" and verdict_["permissionDecision"] == "deny"
    reason = verdict_["permissionDecisionReason"]
    assert reason.startswith("PHI guard: `cat` prints patient file content.")
    assert SAN in reason and "phi-paths.txt" in reason and F not in reason
    for tool, tool_input in (
        ("Read", {"file_path": F}),
        ("Grep", {"pattern": "x", "path": P}),
        ("Glob", {"pattern": f"{P}/*"}),
    ):
        out = run(
            {"tool_name": tool, "tool_input": tool_input, "cwd": ELSEWHERE}, tmp_path, listing
        )
        assert out["hookSpecificOutput"]["permissionDecision"] == "deny", tool


def test_the_hook_prints_an_empty_object_to_allow(tmp_path):
    listing = listing_in(tmp_path, LIST)
    assert (
        run(
            {"tool_name": "Bash", "tool_input": {"command": "ls"}, "cwd": ELSEWHERE},
            tmp_path,
            listing,
        )
        == {}
    )
    assert (
        run(
            {"tool_name": "Write", "tool_input": {"file_path": F}, "cwd": ELSEWHERE},
            tmp_path,
            listing,
        )
        == {}
    )


@pytest.mark.parametrize(
    "text",
    [
        "",
        "not json",
        "[]",
        "null",
        "7",
        '"cat"',
        "{}",
        '{"tool_name": "Bash"}',
        '{"tool_name": "Bash", "tool_input": null}',
        '{"tool_name": "Bash", "tool_input": "cat"}',
        '{"tool_name": "Read", "tool_input": []}',
        '{"tool_name": 5}',
    ],
)
def test_the_hook_passes_when_the_input_is_not_a_tool_call(text, tmp_path):
    assert run(text, tmp_path, listing_in(tmp_path, LIST)) == {}


def test_a_call_with_no_folder_uses_the_hook_folder(tmp_path):
    call = {"tool_name": "Bash", "tool_input": {"command": f"cat {F}"}}
    out = run(call, tmp_path, listing_in(tmp_path, LIST))
    assert out["hookSpecificOutput"]["permissionDecision"] == "deny"


def test_without_a_list_the_hook_guards_nothing(tmp_path):
    """No list on the machine means no PHI to guard: nothing changes."""
    call = {"tool_name": "Bash", "tool_input": {"command": f"cat {F}"}, "cwd": ELSEWHERE}
    assert run(call, tmp_path, None, home=str(tmp_path)) == {}  # no .claude/phi-paths.txt here


def test_the_list_in_the_home_folder_is_read_when_no_other_is_named(tmp_path):
    (tmp_path / ".claude").mkdir()
    (tmp_path / ".claude" / "phi-paths.txt").write_text(f"{OUT}/*/panel\n")
    home = str(tmp_path)
    call = {"tool_name": "Bash", "tool_input": {"command": f"cat {F}"}, "cwd": ELSEWHERE}
    assert (
        run(call, tmp_path, None, home=home)["hookSpecificOutput"]["permissionDecision"] == "deny"
    )
    call["tool_input"]["command"] = f"cat {OUT}/lab1/c246/x"
    assert run(call, tmp_path, None, home=home) == {}
    # the default sanitizer sits below the home folder
    san = f"{tmp_path}/test/gbcms/harness/hpc/kit/phi_sanitize.py"
    call["tool_input"]["command"] = f"python3 -I {RUN} {F} 2>&1 | python3 -I {san}"
    assert run(call, tmp_path, None, home=home) == {}


def test_the_environment_variable_picks_the_list(tmp_path):
    other = listing_in(tmp_path, f"{PLAN}/across\n")
    call = {"tool_name": "Bash", "tool_input": {"command": f"cat {F}"}, "cwd": ELSEWHERE}
    assert run(call, tmp_path, other) == {}
    call["tool_input"]["command"] = f"cat {PLAN}/across/x"
    assert run(call, tmp_path, other)["hookSpecificOutput"]["permissionDecision"] == "deny"


def test_a_broken_list_denies_through_the_process_too(tmp_path):
    call = {"tool_name": "Bash", "tool_input": {"command": "ls"}, "cwd": ELSEWHERE}
    out = run(call, tmp_path, listing_in(tmp_path, "not/absolute\n"))
    assert "unusable" in out["hookSpecificOutput"]["permissionDecisionReason"]


def test_the_hook_source_is_plain_stdlib_python_for_old_interpreters():
    """The tool runs `python3` from the PATH; the macOS system Python is 3.9."""
    tree_ = ast.parse(HOOK.read_text(), feature_version=(3, 9))
    imported = {
        a.name.split(".")[0] for n in ast.walk(tree_) if isinstance(n, ast.Import) for a in n.names
    }
    imported |= {
        n.module.split(".")[0]
        for n in ast.walk(tree_)
        if isinstance(n, ast.ImportFrom) and n.module
    }
    assert imported <= set(sys.stdlib_module_names)
