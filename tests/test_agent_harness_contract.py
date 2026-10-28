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
    problems = checker.module_map_problems(
        ROOT / ".agents" / "rules" / "architecture.md", checker.tracked()
    )
    assert not problems, problems
