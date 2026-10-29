"""The agent harness stays true to the code.

Agents read the skills (`.agents/skills/*/SKILL.md`) and the rules (`.agents/rules/*.md`)
on trust. Two checks keep that trust earned:

* every skill's frontmatter parses as YAML and carries its whole description (an unquoted
  ``#`` or ``: `` truncates it, and the skill's "use when" trigger never reaches the model);
  the same checker runs in ``.githooks/pre-commit``;
* every backticked name in a skill or rule resolves: a path exists, a flag is in the
  source, an identifier is in ``src/`` or ``rust/src/``, and the module map in
  ``architecture.md`` lists the files that exist (``scripts/check_agent_names.py``).
"""

import importlib.util
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
SKILLS = sorted((ROOT / ".agents" / "skills").glob("*/SKILL.md"))
CHECKER = ROOT / "scripts" / "check_agent_names.py"
FRONTMATTER = ROOT / ".githooks" / "check_skill_frontmatter.py"


def _checker():
    spec = importlib.util.spec_from_file_location("check_agent_names", CHECKER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _head(skill: Path) -> dict:
    text = skill.read_text()
    assert text.startswith("---\n"), f"{skill}: no frontmatter"
    return yaml.safe_load(text[4:].split("\n---\n", 1)[0])


def test_the_skills_are_found():
    """A glob that matches nothing would pass every check below."""
    assert len(SKILLS) >= 8, f"expected the project's skills under {ROOT / '.agents'}, got {SKILLS}"


@pytest.mark.parametrize("skill", SKILLS, ids=lambda p: p.parent.name)
def test_skill_frontmatter_parses_with_the_hook_checker(skill):
    """The pre-commit checker accepts every skill as committed."""
    done = subprocess.run(
        [sys.executable, "-I", str(FRONTMATTER), str(skill.relative_to(ROOT))],
        input=skill.read_text(),
        capture_output=True,
        text=True,
    )
    assert done.returncode == 0, done.stderr


@pytest.mark.parametrize("skill", SKILLS, ids=lambda p: p.parent.name)
def test_skill_description_is_whole(skill):
    """The description reaches its "Use when" trigger and is not cut at a YAML comment."""
    meta = _head(skill)
    assert meta["name"] == skill.parent.name
    description = meta["description"]
    assert isinstance(description, str)
    assert re.search(
        r"\buse (when|before)\b", description.lower()
    ), f"{skill.parent.name}: no 'Use when' trigger"
    assert description.rstrip().endswith((".", ")")), f"{skill.parent.name}: cut short?"


def test_checker_rejects_an_unquoted_description(tmp_path):
    """The frontmatter checker refuses what truncates a description."""
    bad = "---\nname: x\ndescription: Reference for #[pyo3(get)] fields. Use when editing.\n---\n"
    done = subprocess.run(
        [sys.executable, "-I", str(FRONTMATTER), ".agents/skills/x/SKILL.md"],
        input=bad,
        capture_output=True,
        text=True,
    )
    assert done.returncode != 0 and "quote the description" in done.stderr


def test_every_name_in_the_skills_and_rules_resolves():
    done = subprocess.run([sys.executable, str(CHECKER)], capture_output=True, text=True, cwd=ROOT)
    assert done.returncode == 0, "names that do not resolve:\n" + done.stdout


def test_name_checker_flags_a_missing_name_flag_path_and_constant(tmp_path):
    """The checker would fail on drift: each kind of dead name is reported, live ones pass."""
    skill = tmp_path / "x" / "SKILL.md"
    skill.parent.mkdir()
    skill.write_text(
        '---\nname: x\ndescription: "d. Use when x."\n---\n'
        "Live: `check_allele_with_qual`, `--umi-tag`, `rust/src/shared/fragment.rs`, "
        "`BIN_WINDOW = 10_000`, `window::allele_kind(ref, alt)`.\n"
        "Dead: `no_such_function_anywhere`, `--no-such-flag`, `rust/src/no_such_file.rs`, "
        "`BIN_WINDOW = 7`, `NoSuchTypeAnywhere`.\n"
    )
    done = subprocess.run(
        [sys.executable, str(CHECKER), str(skill)], capture_output=True, text=True, cwd=ROOT
    )
    assert done.returncode == 1
    reported = done.stdout
    for dead in (
        "`no_such_function_anywhere`",
        "`--no-such-flag`",
        "`rust/src/no_such_file.rs`",
        "constant BIN_WINDOW is not 7",
        "`NoSuchTypeAnywhere`",
    ):
        assert dead in reported, (dead, reported)
    for live in ("check_allele_with_qual", "--umi-tag", "fragment.rs", "allele_kind"):
        assert f"`{live}" not in reported, (live, reported)
    assert "BIN_WINDOW = 10_000" not in reported


def test_the_module_map_lists_the_modules_that_exist():
    """architecture.md names every Rust and Python module file, and only those."""
    checker = _checker()
    text = (ROOT / ".agents" / "rules" / "architecture.md").read_text()
    problems = checker.module_map_problems(text, checker.Source().names)
    assert not problems, problems


# ── a shortened flag, a removed duplicate-name module, staged against unstaged ──────────


def _check(skill: Path, *flags: str, cwd: Path = ROOT, script: Path = CHECKER):
    return subprocess.run(
        [sys.executable, str(script), *flags, str(skill)],
        capture_output=True,
        text=True,
        cwd=cwd,
    )


def test_a_flag_is_matched_whole_not_as_the_start_of_a_longer_flag(tmp_path):
    """`--mfsd-report-min` must not pass because `--mfsd-report-min-alt` exists."""
    skill = tmp_path / "x" / "SKILL.md"
    skill.parent.mkdir()
    skill.write_text(
        '---\nname: x\ndescription: "d. Use when x."\n---\n'
        "Live: `--mfsd-report-min-alt`, `--mfsd-report`, `--no-baq` (spelled `--apply-baq/--no-baq`).\n"
        "Dead: `--mfsd-report-min`, `--mfsd-rep`.\n"
    )
    done = _check(skill)
    assert done.returncode == 1
    assert "flag --mfsd-report-min is not in the code" in done.stdout
    assert "flag --mfsd-rep is not in the code" in done.stdout
    for live in ("`--mfsd-report-min-alt`", "`--mfsd-report`", "`--no-baq`"):
        assert live not in done.stdout, done.stdout


MODULE_MAP = """## Project Layout

```
gbcms/
├── rust/
│   └── src/
│       ├── counting/
│       │   └── engine.rs        # counting loop
│       └── normalize/
│           ├── mod.rs, engine.rs  # entry, preparation engine
│           └── types.rs
└── src/gbcms/
    └── cli.py
```
"""


def test_the_module_map_tells_two_files_of_one_name_apart():
    """`counting/engine.rs` and `normalize/engine.rs` are two files; neither stands in for the other."""
    checker = _checker()
    everything = [
        "rust/src/counting/engine.rs",
        "rust/src/normalize/engine.rs",
        "rust/src/normalize/mod.rs",
        "rust/src/normalize/types.rs",
        "src/gbcms/cli.py",
    ]
    assert checker.module_map_problems(MODULE_MAP, everything) == []
    assert checker.map_paths(MODULE_MAP.split("```")[1]) == {
        "rust/src/counting/engine.rs",
        "rust/src/normalize/engine.rs",
        "rust/src/normalize/mod.rs",
        "rust/src/normalize/types.rs",
        "src/gbcms/cli.py",
    }
    # The normalize module is gone but counting/engine.rs remains: matching by name alone
    # would hide it.
    without = [n for n in everything if n != "rust/src/normalize/engine.rs"]
    assert checker.module_map_problems(MODULE_MAP, without) == [
        "module map lists rust/src/normalize/engine.rs, which is not a file"
    ]
    # And the other way: a new normalize/engine.rs, with only counting/engine.rs in the map.
    short_map = MODULE_MAP.replace(
        "           ├── mod.rs, engine.rs  # entry, preparation engine\n", "           ├── mod.rs\n"
    )
    assert checker.module_map_problems(short_map, everything) == [
        "module map omits rust/src/normalize/engine.rs"
    ]


GOOD_SKILL = '---\nname: foo\ndescription: "Reference. Use when testing."\n---\nCalls `{name}`.\n'


def _git(repo: Path, *args: str, check: bool = True):
    return subprocess.run(["git", *args], cwd=repo, capture_output=True, text=True, check=check)


@pytest.fixture
def repo(tmp_path):
    """A throwaway repo with the hook, the checker and one function in src/ (no real data)."""
    pytest.importorskip("yaml")
    if shutil.which("git") is None or shutil.which("bash") is None:
        pytest.skip("needs git and bash")
    _git(tmp_path, "init", "-q")
    for key, value in (
        ("user.email", "t@example.org"),
        ("user.name", "t"),
        ("commit.gpgsign", "false"),
        ("core.hooksPath", ".githooks"),
    ):
        _git(tmp_path, "config", key, value)
    shutil.copytree(ROOT / ".githooks", tmp_path / ".githooks")
    (tmp_path / "scripts").mkdir()
    shutil.copy(CHECKER, tmp_path / "scripts" / "check_agent_names.py")
    (tmp_path / "src" / "gbcms").mkdir(parents=True)
    (tmp_path / "src" / "gbcms" / "x.py").write_text("def real_function():\n    return 1\n")
    _git(tmp_path, "add", "scripts", "src")
    first = _git(tmp_path, "commit", "-q", "-m", "base", check=False)
    assert first.returncode == 0, first.stderr
    return tmp_path


def _skill(repo: Path, name: str) -> Path:
    path = repo / ".agents" / "skills" / "foo" / "SKILL.md"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(GOOD_SKILL.format(name=name))
    return path


def test_the_hook_judges_the_staged_code_not_the_working_tree(repo):
    """A function added only in the working tree does not satisfy a committed skill."""
    skill = _skill(repo, "unstaged_only_function")
    _git(repo, "add", str(skill.relative_to(repo)))
    code = repo / "src" / "gbcms" / "x.py"
    code.write_text(code.read_text() + "\n\ndef unstaged_only_function():\n    return 2\n")
    # The working tree has the function, so the working-tree check passes ...
    assert _check(skill, cwd=repo, script=repo / "scripts" / "check_agent_names.py").returncode == 0
    # ... the staged check does not, and nor does the hook.
    staged = _check(skill, "--staged", cwd=repo, script=repo / "scripts" / "check_agent_names.py")
    assert staged.returncode == 1 and "unstaged_only_function" in staged.stdout
    blocked = _git(repo, "commit", "-q", "-m", "skill", check=False)
    assert blocked.returncode != 0 and "does not exist" in blocked.stderr, blocked.stderr


def test_the_hook_judges_the_staged_skill_not_the_working_tree_copy(repo):
    """Fixing the skill only in the working tree does not let the staged one through."""
    skill = _skill(repo, "missing_function")
    _git(repo, "add", str(skill.relative_to(repo)))
    skill.write_text(GOOD_SKILL.format(name="real_function"))  # an unstaged fix
    assert _check(skill, cwd=repo, script=repo / "scripts" / "check_agent_names.py").returncode == 0
    staged = _check(skill, "--staged", cwd=repo, script=repo / "scripts" / "check_agent_names.py")
    assert staged.returncode == 1 and "missing_function" in staged.stdout
    blocked = _git(repo, "commit", "-q", "-m", "skill", check=False)
    assert blocked.returncode != 0, "the hook committed a skill naming a missing function"


def test_the_hook_lets_a_good_staged_skill_through_whatever_the_working_tree_holds(repo):
    """The staged skill and the staged code agree, so the commit goes in."""
    skill = _skill(repo, "real_function")
    _git(repo, "add", str(skill.relative_to(repo)))
    code = repo / "src" / "gbcms" / "x.py"
    code.write_text("def renamed_in_the_working_tree():\n    return 1\n")  # unstaged
    done = _git(repo, "commit", "-q", "-m", "skill", check=False)
    assert done.returncode == 0, done.stderr
