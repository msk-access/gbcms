"""The user docs stay true to the code.

Three checks, each parsing the docs the way a reader copies from them:

* every fenced ```bash example that runs ``gbcms dna|rna|merge|normalize|convert``
  parses against the Typer app: every flag exists, every required option is there;
* the MAF column and VCF INFO/FORMAT names the Output Formats tables list are the
  names the writers emit for the options each table is about, and the other way
  round (a field the writer emits and the table lacks is a finding too);
* the Nextflow parameters table matches ``nextflow.config``: names and defaults.

The flag vocabulary beside them (QC flags) is pinned in
``test_merge_outputs_group5_contract.py``.
"""

import pathlib
import re
import shlex

import pytest
import typer.main

from gbcms.cli import app
from gbcms.io.output import MafWriter, VcfWriter

ROOT = pathlib.Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
OUTPUT_FORMATS = DOCS / "reference" / "output-formats.md"
NEXTFLOW_CONFIG = ROOT / "nextflow" / "nextflow.config"
NEXTFLOW_PARAMS = DOCS / "nextflow" / "parameters.md"


# ── (a) fenced gbcms examples parse against the Typer app ─────────────────────

_SUBCOMMANDS = ("dna", "rna", "merge", "normalize", "convert")
_ANNOTATION = re.compile(r"\s*#\s*\(\d+\)!\s*$")  # mkdocs code annotation: `# (1)!`
_PLACEHOLDER = re.compile(r"^<[^>]+>(\.\.\.)?$")  # a synopsis operand: `<FILE>`, `<NAME:PATH>...`
_SHELL_BREAK = {"|", "||", "&&", ";", "&", ">", ">>", "<", "2>&1", "2>", "&>"}
# Docs that quote commands which are not examples to run: the changelog is history.
_NOT_EXAMPLES = {"CHANGELOG.md"}


def _doc_pages() -> list[pathlib.Path]:
    pages = [p for p in sorted(DOCS.rglob("*.md")) if p.name not in _NOT_EXAMPLES]
    return [*pages, ROOT / "README.md"]


def _bash_commands(page: pathlib.Path) -> list[tuple[int, str]]:
    """``(line, command)`` for each logical command line of the page's ```bash blocks:
    backslash continuations joined, comments and mkdocs annotations dropped."""
    commands: list[tuple[int, str]] = []
    in_bash = False
    in_other = False
    pending: list[str] = []
    start = 0
    for number, raw in enumerate(page.read_text().splitlines(), 1):
        line = raw.strip()
        if line.startswith("```"):
            if in_bash or in_other:
                if pending:  # a command left open at the fence: keep what it had
                    commands.append((start, " ".join(pending)))
                    pending = []
                in_bash = in_other = False
            else:
                in_bash = line[3:].strip().startswith("bash")
                in_other = not in_bash
            continue
        if not in_bash:
            continue
        line = _ANNOTATION.sub("", line)
        if not pending and (not line or line.startswith("#")):
            continue
        if not pending:
            start = number
        continued = line.endswith("\\")
        pending.append(line[:-1].rstrip() if continued else line)
        if not continued:
            commands.append((start, " ".join(pending)))
            pending = []
    return commands


def _gbcms_invocation(command: str) -> tuple[str, list[str], bool] | None:
    """``(subcommand, args, elided)`` for a command that runs one of the subcommands, from
    its ``gbcms <sub>`` token on (so ``docker run ... gbcms dna ...`` counts), cut at the
    first shell operator. ``elided`` when the example abbreviates its options with
    ``...``. Synopsis placeholders (``[OPTIONS]``, ``<FILE>``) become dummy values."""
    try:
        tokens = shlex.split(command, comments=True)
    except ValueError:
        return None
    for i, token in enumerate(tokens[:-1]):
        if token != "gbcms" or tokens[i + 1] not in _SUBCOMMANDS:
            continue
        args: list[str] = []
        elided = False
        for arg in tokens[i + 2 :]:
            if arg in _SHELL_BREAK or arg.startswith(("2>", ">")):
                break
            if arg == "...":
                elided = True
            elif arg == "[OPTIONS]":
                continue
            elif _PLACEHOLDER.match(arg):
                args.append("PLACEHOLDER")
            else:
                args.append(arg)
        return tokens[i + 1], args, elided
    return None


@pytest.fixture
def commands(monkeypatch):
    """The app's commands with their path existence checks off: an example names files
    that are not here. (Typer bundles its own click, so its errors are matched by
    name, not by class.)"""
    group = typer.main.get_command(app)
    for command in group.commands.values():  # type: ignore[attr-defined]
        for param in command.params:
            if getattr(param.type, "exists", False):
                monkeypatch.setattr(param.type, "exists", False)
    return group.commands  # type: ignore[attr-defined]


def _check_example(commands, sub: str, args: list[str], elided: bool) -> str | None:
    """Why the app would refuse the example, or None."""
    try:
        ctx = commands[sub].make_context(sub, list(args))
    except Exception as e:  # a usage error of whichever click typer runs on
        if type(e).__name__ == "Exit":  # `--help`: the flag exists
            return None
        if not hasattr(e, "format_message"):
            raise
        if type(e).__name__ == "MissingParameter":
            return None if elided else f"missing {getattr(e.param, 'name', e)}"
        return e.format_message()
    if sub in ("dna", "rna") and not elided:
        if not (ctx.params.get("bam_files") or ctx.params.get("bam_list")):
            return "needs --bam or --bam-list"
    return None


def test_the_doc_examples_find_the_commands():
    """The scan is not vacuous: it sees the examples of every command and the
    synopses."""
    seen: dict[str, int] = {}
    for page in _doc_pages():
        for _, command in _bash_commands(page):
            found = _gbcms_invocation(command)
            if found:
                seen[found[0]] = seen.get(found[0], 0) + 1
    assert all(seen.get(sub, 0) >= 1 for sub in _SUBCOMMANDS), seen
    assert seen["dna"] > 10 and seen["rna"] > 5, seen


def test_every_gbcms_example_in_the_docs_parses(commands):
    """Every flag exists and every required option is there: a copied example runs
    as far as argument parsing (an example abbreviating with ``...`` may omit
    required options, never use an unknown flag)."""
    failures = []
    for page in _doc_pages():
        for line, command in _bash_commands(page):
            found = _gbcms_invocation(command)
            if found is None:
                continue
            why = _check_example(commands, *found)
            if why:
                failures.append(f"{page.relative_to(ROOT)}:{line}: {why}: {command[:90]}")
    assert failures == [], "\n".join(failures)


# ── (b) documented output names exist in the writers' sets ─────────────────────


def _doc_tables(page: pathlib.Path) -> dict[tuple[str, str], list[str]]:
    """The first-cell names of the page's tables, by ``(heading, tab)``: the nearest
    heading above (``##`` or ``###``; a ``####`` sub-heading stays in its section), and
    the nearest ``=== "tab"`` or ``??? note "title"`` below it."""
    tables: dict[tuple[str, str], list[str]] = {}
    heading = tab = ""
    for line in page.read_text().splitlines():
        h = re.match(r"^(#{2,6})\s+(.*)", line)
        if h:
            tab = ""
            if len(h.group(1)) <= 3:
                heading = h.group(2).strip()
            continue
        t = re.match(r'^\s*(?:===|\?\?\?)\s+(?:note\s+)?"(.*)"', line)
        if t:
            tab = t.group(1)
            continue
        # A name, not a worked example (`461 TAA > T`): `{prefix}` is the column prefix.
        row = re.match(r"^\s*\|\s*`([\w{}]+)`\s*\|", line)
        if row:
            tables.setdefault((heading, tab), []).append(row.group(1))
    assert tables, f"no tables parsed from {page}"
    return tables


def _table(tables, heading: str, tab: str = "") -> set[str]:
    """The names under a heading (prefix match) and tab (prefix match), unioned."""
    found: set[str] = set()
    for (h, t), names in tables.items():
        if h.startswith(heading) and t.startswith(tab):
            found.update(names)
    return found  # empty when the page has no such table: every emitted name is then reported


def _vcf_header_ids(tmp_path, maf_input: bool = False, **options) -> tuple[set[str], set[str]]:
    """The ``##INFO`` and ``##FORMAT`` IDs the VCF writer declares for these options."""
    path = tmp_path / "header.vcf"
    writer = VcfWriter(path, **options)
    writer._write_header(maf_input=maf_input)
    writer.close()
    text = path.read_text()
    return (
        set(re.findall(r"^##INFO=<ID=([^,>]+)", text, re.M)),
        set(re.findall(r"^##FORMAT=<ID=([^,>]+)", text, re.M)),
    )


def test_the_documented_vcf_info_fields_are_the_headers_the_writer_declares(tmp_path):
    tables = _doc_tables(OUTPUT_FORMATS)

    def info(**kw):
        return _vcf_header_ids(tmp_path, **kw)[0]

    base = info()
    rna = info(mode="rna")
    cases = {
        "Always present": base,
        "--rescue-mnp only": info(rescue_mnp=True) - base,
        "RNA mode only": rna - base,
        "--mfsd only": info(mfsd=True) - base,
        "--gtf only": info(mode="rna", has_gtf=True) - rna,
        "--show-normalization only": info(show_normalization=True) - base,
        "MAF input": info(maf_input=True) - base,
    }
    problems = []
    for tab, expected in cases.items():
        documented = _table(tables, "INFO Fields", tab)
        problems += [f"{tab}: documented, not declared: {n}" for n in sorted(documented - expected)]
        problems += [f"{tab}: declared, not documented: {n}" for n in sorted(expected - documented)]
    assert problems == [], "\n".join(problems)


def test_the_documented_vcf_format_fields_are_the_headers_the_writer_declares(tmp_path):
    tables = _doc_tables(OUTPUT_FORMATS)
    dna = _vcf_header_ids(tmp_path)[1]
    rna = _vcf_header_ids(tmp_path, mode="rna")[1]
    problems = []
    for tab, expected in {"DNA mode": dna, "RNA mode": rna - dna}.items():
        documented = _table(tables, "FORMAT Fields", tab)
        problems += [f"{tab}: documented, not declared: {n}" for n in sorted(documented - expected)]
        problems += [f"{tab}: declared, not documented: {n}" for n in sorted(expected - documented)]
    assert problems == [], "\n".join(problems)


def test_the_documented_maf_columns_are_the_columns_the_writer_emits():
    tables = _doc_tables(OUTPUT_FORMATS)
    names = MafWriter.column_names
    dna = set(names())
    rna = set(names(mode="rna"))
    cases = {
        ("Two Output Paths", "VCF → MAF"): set(MafWriter.vcf_input_headers()),
        ("gbcms Count Columns", "Default (no prefix)"): set(names(rescue_mnp=True)),
        ("RNA-Specific MAF Columns", ""): rna - dna,
        ("GTF-Aware MAF Columns", ""): set(names(mode="rna", has_gtf=True)) - rna,
        ("Library Type", "mFSD MAF Columns"): set(names(mfsd=True)) - dna,
        ("Library Type", "Normalization MAF Columns"): set(names(show_normalization=True)) - dna,
    }
    problems = []
    for (heading, tab), expected in cases.items():
        documented = {n.replace("{prefix}", "") for n in _table(tables, heading, tab)}
        label = f"{heading} / {tab}" if tab else heading
        problems += [
            f"{label}: documented, not emitted: {n}" for n in sorted(documented - expected)
        ]
        problems += [
            f"{label}: emitted, not documented: {n}" for n in sorted(expected - documented)
        ]
    assert problems == [], "\n".join(problems)


# Backticked names of the column families that are not MAF columns of a counting run:
# the `gbcms normalize` TSV's columns, a normalization parameter, a GTF attribute.
_NOT_MAF_COLUMNS = {"norm_pos", "norm_ref", "norm_alt", "norm_window", "transcript_id"}


def test_every_column_name_the_docs_mention_is_one_gbcms_writes():
    """Prose too, not only the tables: a renamed or misremembered column (`norm_chrom`,
    `rna_splice_spanning_count`) in an example's note or a feature list."""
    emitted = set(
        MafWriter.column_names(
            mfsd=True, mode="rna", has_gtf=True, rescue_mnp=True, show_normalization=True
        )
    )
    family = re.compile(r"`((?:rna|norm|asjd|mfsd|exon|transcript)_\w+)`")
    unknown = []
    for page in _doc_pages():
        for n, line in enumerate(page.read_text().splitlines(), 1):
            unknown += [
                f"{page.relative_to(ROOT)}:{n}: {name}"
                for name in family.findall(line)
                if name not in emitted and name not in _NOT_MAF_COLUMNS
            ]
    assert unknown == [], unknown


def test_the_docs_count_the_mfsd_fields_the_writers_emit(tmp_path):
    """The numbers the prose quotes ("41 mFSD columns", "13 `MFSD_*` INFO fields")."""
    maf = len(set(MafWriter.column_names(mfsd=True)) - set(MafWriter.column_names()))
    info = len(_vcf_header_ids(tmp_path, mfsd=True)[0] - _vcf_header_ids(tmp_path)[0])
    quoted = [
        (maf, re.compile(r"\b(\d+) mFSD (?:MAF )?columns")),
        (info, re.compile(r"\b(\d+) `MFSD_\*` INFO fields")),
    ]
    wrong = []
    for page in _doc_pages():
        for n, line in enumerate(page.read_text().splitlines(), 1):
            for emitted, pattern in quoted:
                wrong += [
                    f"{page.relative_to(ROOT)}:{n}: {count}, emitted {emitted}"
                    for count in pattern.findall(line)
                    if int(count) != emitted
                ]
    assert (maf, info) == (41, 13)
    assert wrong == [], wrong


# ── (c) the Nextflow parameters table matches nextflow.config ─────────────────

# Declared in the config, documented nowhere on purpose: Nextflow's own `--help`
# convention; nothing in the pipeline reads it.
_UNDOCUMENTED_PARAMS = {"help"}


def _config_params() -> dict[str, str]:
    """``params { ... }`` of nextflow.config: name -> default, quotes and comments off."""
    params: dict[str, str] = {}
    inside = False
    for line in NEXTFLOW_CONFIG.read_text().splitlines():
        if re.match(r"^params\s*\{", line):
            inside = True
            continue
        if inside and line.startswith("}"):
            break
        m = re.match(r"^\s+(\w+)\s*=\s*(.*?)\s*(?://.*)?$", line) if inside else None
        if m:
            params[m.group(1)] = m.group(2).strip().strip("'\"")
    return params


def _documented_params() -> dict[str, str | None]:
    """The parameters table rows: ``--name`` -> default (None where the table has no
    Default column, the required parameters)."""
    documented: dict[str, str | None] = {}
    has_default = False
    for line in NEXTFLOW_PARAMS.read_text().splitlines():
        if not line.startswith("|"):
            has_default = False
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if cells[0] == "Parameter":
            has_default = len(cells) > 2 and cells[1] == "Default"
            continue
        m = re.match(r"^`--(\w+)`$", cells[0])
        if m:
            default = cells[1].strip("`").strip("'\"") if has_default else None
            documented[m.group(1)] = default
    return documented


def test_the_nextflow_parameters_table_matches_the_config():
    config = _config_params()
    documented = _documented_params()
    assert len(config) > 40 and len(documented) > 40
    problems = [
        f"documented, not in nextflow.config: --{n}" for n in sorted(set(documented) - set(config))
    ]
    problems += [
        f"in nextflow.config, not documented: --{n}"
        for n in sorted(set(config) - set(documented) - _UNDOCUMENTED_PARAMS)
    ]
    for name, default in documented.items():
        if default is not None and name in config and default != config[name]:
            problems.append(f"--{name}: documented default {default!r}, config {config[name]!r}")
    assert problems == [], "\n".join(problems)


@pytest.mark.parametrize("name", ["input", "variants", "fasta"])
def test_the_required_nextflow_parameters_have_no_default(name):
    assert _config_params()[name] == "null"
    assert _documented_params()[name] is None
