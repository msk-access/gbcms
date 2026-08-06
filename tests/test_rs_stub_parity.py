"""`_rs.pyi` must match the runtime pyo3 signatures (AGENTS.md invariant 5).

A stub that grants a default the extension does not actually have is invisible to every tool
that would normally catch it: mypy accepts the call, ruff accepts it, the import succeeds, and
it fails only at runtime with ``missing N required positional arguments``. That drift shipped
twice before this test existed — ``count_bam`` and ``count_bam_binned`` each declared **nine**
phantom defaults (``min_mapq`` … ``threads``) that the pyo3 signature never granted.

Keeping the stub honest was previously a manual review item. pyo3 publishes the real signature
as ``__text_signature__``, so it can simply be checked.

What is compared: parameter **names**, their **order**, and **which** parameters carry a
default. Not the default's rendered value — pyo3 renders non-literal defaults as ``...``
(``sibling_variants=Vec::new()``), so comparing values would fail on correct stubs. Names,
order, and default-presence are exactly what a caller's code depends on.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from gbcms import _rs

STUB_PATH = Path(__file__).resolve().parents[1] / "src" / "gbcms" / "_rs.pyi"

# `count_bam` is behind the default-on `legacy-parity` cargo feature and is absent from the
# shipped wheel. Any *other* missing name is a genuine stub/runtime mismatch, not a build
# variant, so only this one may be skipped.
FEATURE_GATED = {"count_bam"}

# Guards against the test quietly degenerating to comparing nothing (e.g. a stub-parse change
# that yields an empty list would otherwise "pass").
MIN_FUNCTIONS_COMPARED = 6


def _split_params(sig: str) -> list[str]:
    """Split a ``__text_signature__`` body on top-level commas.

    Depth-aware so a default containing a comma (a tuple or list) does not split the
    parameter in half. pyo3 emits no such default today; this keeps the parser from being
    quietly wrong if one appears.
    """
    out: list[str] = []
    depth = 0
    quote: str | None = None
    current = ""
    for ch in sig:
        if quote:
            current += ch
            if ch == quote:
                quote = None
            continue
        if ch in "\"'":
            quote = ch
            current += ch
        elif ch in "([{":
            depth += 1
            current += ch
        elif ch in ")]}":
            depth -= 1
            current += ch
        elif ch == "," and depth == 0:
            out.append(current.strip())
            current = ""
        else:
            current += ch
    if current.strip():
        out.append(current.strip())
    return out


def _runtime_signature(name: str) -> tuple[list[str], set[str]]:
    """(ordered parameter names, names carrying a default) from pyo3's `__text_signature__`."""
    text = getattr(_rs, name).__text_signature__
    assert text, f"{name} has no __text_signature__ to check against"
    params = _split_params(text.strip().removeprefix("(").removesuffix(")"))
    names, defaulted = [], set()
    for p in params:
        if p in ("/", "*") or p.startswith(("*", "$")):
            continue  # positional/keyword markers and `$self` are not parameters
        pname, _, has_default = p.partition("=")
        pname = pname.strip()
        names.append(pname)
        if has_default or p.endswith("="):
            defaulted.add(pname)
    return names, defaulted


def _stub_signatures() -> dict[str, tuple[list[str], set[str]]]:
    """(ordered parameter names, names carrying a default) for each `def` in the stub."""
    tree = ast.parse(STUB_PATH.read_text())
    out: dict[str, tuple[list[str], set[str]]] = {}
    for node in tree.body:
        if not isinstance(node, ast.FunctionDef):
            continue
        args = node.args
        assert not args.kwonlyargs, (
            f"{node.name}: this check does not model keyword-only params; extend it rather "
            f"than letting them go unchecked"
        )
        names = [a.arg for a in args.args]
        # ast aligns defaults to the TAIL of the positional list
        defaulted = {a.arg for a in args.args[len(args.args) - len(args.defaults) :]}
        out[node.name] = (names, defaulted)
    return out


STUBS = _stub_signatures()


def test_every_stub_function_exists_at_runtime():
    """A stub for a function the extension does not export is dead documentation."""
    missing = [n for n in STUBS if not hasattr(_rs, n)]
    assert (
        missing == [] or set(missing) <= FEATURE_GATED
    ), f"declared in _rs.pyi but not exported by the extension: {missing}"


@pytest.mark.parametrize("name", sorted(STUBS))
def test_stub_matches_runtime_signature(name):
    """Parameter names, order, and which ones have defaults must agree exactly."""
    if not hasattr(_rs, name):
        pytest.skip(f"{name} is behind the legacy-parity feature and absent from this build")

    stub_names, stub_defaulted = STUBS[name]
    rt_names, rt_defaulted = _runtime_signature(name)

    assert (
        stub_names == rt_names
    ), f"{name}: parameter names/order differ.\n  stub:    {stub_names}\n  runtime: {rt_names}"

    phantom = sorted(stub_defaulted - rt_defaulted)
    assert not phantom, (
        f"{name}: the stub grants defaults the extension does not have: {phantom}. "
        f"mypy would accept a call omitting them; it fails at runtime with "
        f"'missing required positional arguments'."
    )
    undocumented = sorted(rt_defaulted - stub_defaulted)
    assert not undocumented, (
        f"{name}: the extension defaults these but the stub marks them required: "
        f"{undocumented}. Callers are forced to pass values they need not supply."
    )


def test_the_check_actually_covered_the_module():
    """Without this, a parsing regression that finds no functions would 'pass' silently."""
    compared = [n for n in STUBS if hasattr(_rs, n)]
    assert (
        len(compared) >= MIN_FUNCTIONS_COMPARED
    ), f"only {len(compared)} functions compared ({compared}); the stub parser is likely broken"
