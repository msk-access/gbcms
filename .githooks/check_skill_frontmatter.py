"""Check that a skill's SKILL.md frontmatter parses as YAML. Used by .githooks/pre-commit.

Usage: git show :<path> | python3 -I check_skill_frontmatter.py <path>

An unquoted `#` or `: ` in a description (`#[pyo3(get)]`, `NOTE: `) silently truncates
it (the rest becomes a YAML comment) or breaks the file, so the skill's "use when"
trigger never reaches the model. Two checks follow. With PyYAML (a test dependency) the
frontmatter must parse with `yaml.safe_load`; without it a minimal check applies. Either
way the description must be written as a quoted string (or a block scalar): a plain
scalar can be truncated without any parse error, so it is refused outright. Exits
non-zero, with the reason on stderr, on any problem.
"""

import re
import sys


def fail(path: str, reason: str) -> None:
    """Report the problem and stop."""
    sys.exit(f"{path}: {reason}")


def main(path: str) -> None:
    text = sys.stdin.read()
    if not text.startswith("---\n") or "\n---\n" not in text[4:]:
        fail(path, "no YAML frontmatter between two --- lines")
    head = text[4:].split("\n---\n", 1)[0]
    try:
        import yaml
    except ImportError:
        yaml = None
    if yaml is not None:
        try:
            meta = yaml.safe_load(head)
        except yaml.YAMLError as err:
            fail(path, "frontmatter is not valid YAML: " + str(err).replace("\n", " "))
        if not isinstance(meta, dict):
            fail(path, "frontmatter is not a mapping")
        name, desc = meta.get("name"), meta.get("description")
    else:
        fields = dict(re.findall(r"^(name|description):[ ]*(.*)$", head, re.M))
        name, desc = fields.get("name"), fields.get("description", "")
    if not isinstance(name, str) or not name.strip():
        fail(path, "frontmatter has no name")
    if not isinstance(desc, str) or not desc.strip():
        fail(path, "frontmatter has no description (an unquoted # or colon may have cut it off)")
    raw = re.search(r"^description:[ ]*(.*)$", head, re.M)
    if raw is None or not re.match(r"[\"'|>]", raw.group(1)):
        fail(path, "quote the description (double quotes): a plain scalar is cut at ' #'")
    if yaml is None and not re.fullmatch(
        r"\"(?:[^\"\\]|\\.)*\"|'[^']*'|[|>].*", raw.group(1).strip()
    ):
        fail(
            path, "PyYAML is not installed, and the quoted description is not one well-formed line"
        )
    if name != path.split("/")[-2]:
        fail(path, f"name {name!r} does not match its directory")


if __name__ == "__main__":
    main(sys.argv[1])
