"""The Nextflow modules pass the CLI's own defaults unless a user sets a param.

Every value a module hands the CLI from a param (``--flag ${params.x}``), every
boolean it switches on a param (``asBool(params.x) ? "--a" : "--b"``) and every
option it emits only off a default (``params.x != 'd' ? "--flag ..."``) must, at the
``nextflow.config`` default, leave the command with the CLI default of that option
for its mode (``gbcms dna`` or ``gbcms rna``). A param shared by both modules can
only hold one default, so where the modes' CLI defaults differ the RNA module reads
its own ``rna_*`` param (``rna_min_mapq``, the RNA gap penalties, ``rna_apply_baq``).
"""

import enum
import re
from pathlib import Path

import pytest
from typer.main import get_command

from gbcms.cli import app

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = ROOT / "nextflow"


def _value(text: str):
    text = text.strip()
    if text in ("true", "false"):
        return text == "true"
    if text[:1] in "'\"":
        return text[1:-1]
    try:
        return int(text)
    except ValueError:
        pass
    try:
        return float(text)
    except ValueError:
        return text


def _params() -> dict:
    """The ``params { ... }`` defaults of nextflow.config (literal values only)."""
    lines = (NEXTFLOW / "nextflow.config").read_text().splitlines()
    start = lines.index("params {")
    end = next(i for i in range(start + 1, len(lines)) if lines[i] == "}")
    out = {}
    for line in lines[start + 1 : end]:
        m = re.match(r"\s*(\w+)\s*=\s*('[^']*'|\"[^\"]*\"|[^\s/]+)", line)
        if m:
            out[m.group(1)] = _value(m.group(2))
    return out


def _options(mode: str) -> dict:
    """``--flag`` → (click option, the value the flag sets: None for a value option)."""
    out = {}
    for p in get_command(app).commands[mode].params:
        for o in p.opts:
            out[o] = (p, True if getattr(p, "is_flag", False) else None)
        for o in getattr(p, "secondary_opts", []):
            out[o] = (p, False)
    return out


def _default(option):
    d = option.default
    return d.value if isinstance(d, enum.Enum) else d


def _same(nf, cli) -> bool:
    if isinstance(cli, bool) or isinstance(nf, bool):
        return nf is cli
    if isinstance(cli, int | float) and isinstance(nf, int | float):
        return float(nf) == float(cli)
    return str(nf) == str(cli)


def _divergences(mode: str) -> list[str]:
    script = (NEXTFLOW / "modules" / "local" / "gbcms" / mode / "main.nf").read_text()
    params, options = _params(), _options(mode)
    out = []
    # --flag ${params.x}: the param's default is the value passed (an empty default
    # behind `params.x ? "--flag ${params.x}" : ""` passes nothing: the CLI's applies)
    for flag, name in re.findall(r"(--[\w-]+) \$\{params\.(\w+)\}", script):
        option, _ = options[flag]
        if params[name] == "" and re.search(rf'params\.{name} \? "{flag} ', script):
            continue
        if not _same(params[name], _default(option)):
            out.append(f"{flag} ${{params.{name}}}: {params[name]!r} vs CLI {_default(option)!r}")
    # asBool(params.x) ? "<on>" : "<off>": the branch the default takes
    for neg, name, on, off in re.findall(
        r'(!?)asBool\(params\.(\w+)\)\s*\?\s*"([^"]*)"\s*:\s*"([^"]*)"', script
    ):
        taken = on if bool(params[name]) != bool(neg) else off
        for flag in re.findall(r"--[\w-]+", taken):
            option, sets = options[flag]
            if sets is not None and sets is not _default(option):
                out.append(f"params.{name}={params[name]!r} passes {flag}, off the CLI default")
    # params.x != 'd' ? "--flag ...": emitted only off 'd', which must be the CLI default
    for name, literal, flag in re.findall(r"params\.(\w+) != '([^']*)' \? \"(--[\w-]+)", script):
        option, _ = options[flag]
        if not (params[name] == literal and _same(literal, _default(option))):
            out.append(
                f"{flag}: params.{name}={params[name]!r}, omitted at {literal!r}, CLI {_default(option)!r}"
            )
    return out


@pytest.mark.parametrize("mode", ["dna", "rna"])
def test_module_defaults_are_the_cli_defaults(mode):
    assert _divergences(mode) == []


def test_rna_reads_its_own_params_where_the_modes_differ():
    script = (NEXTFLOW / "modules" / "local" / "gbcms" / "rna" / "main.nf").read_text()
    for shared in ("min_mapq", "apply_baq", "gap_open_prob", "gap_extend_prob"):
        assert re.findall(rf"params\.{shared}\b", script) == [], f"rna passes params.{shared}"
    params = _params()
    rna = _options("rna")
    assert params["rna_apply_baq"] is _default(rna["--apply-baq"][0]) is True
    assert _same(params["rna_gap_open_prob"], _default(rna["--gap-open-prob"][0]))
    assert _same(params["rna_gap_extend_prob"], _default(rna["--gap-extend-prob"][0]))


def test_the_cli_defaults_are_the_models():
    """The RNA CLI defaults the Nextflow params follow are the config model's too."""
    from gbcms.models.core import AlignmentConfig, GbcmsRnaConfig

    rna = _options("rna")
    fields = GbcmsRnaConfig.model_fields
    alignment = fields["alignment"].default_factory()
    assert _default(rna["--gap-open-prob"][0]) == alignment.hmm_gap_open
    assert _default(rna["--gap-extend-prob"][0]) == alignment.hmm_gap_extend
    assert _default(rna["--repeat-gap-open-prob"][0]) == alignment.hmm_gap_open_repeat
    assert _default(rna["--repeat-gap-extend-prob"][0]) == alignment.hmm_gap_extend_repeat
    assert _default(rna["--llr-threshold"][0]) == alignment.hmm_llr_threshold
    assert AlignmentConfig.validate_backend(_default(rna["--alignment-backend"][0])) == (
        alignment.backend
    )
    assert _default(rna["--apply-baq"][0]) is fields["apply_baq"].default
    quality = fields["quality"].default_factory()
    assert _default(rna["--min-mapq"][0]) == quality.min_mapping_quality
