"""Every gbcms module declares ``__all__``, and its names exist.

The code-quality rule is ``__all__`` in every module, matching its public API. Two
modules (``report.mfsd_report`` and ``rescue_audit``) had none, so ``from module import *``
and doc tools took every name, helpers and imports included. This pins the rule for the
whole package and the exact public names of those two.
"""

import importlib
import pkgutil

import pytest

import gbcms

# walk_packages yields the submodules only, so the package root is named here. The compiled
# extension (typed by _rs.pyi) is not Python source and has no __all__.
MODULES = [
    "gbcms",
    *sorted(
        m.name for m in pkgutil.walk_packages(gbcms.__path__, "gbcms.") if m.name != "gbcms._rs"
    ),
]


@pytest.mark.parametrize("name", MODULES)
def test_a_module_declares_all_and_every_name_in_it_exists(name):
    module = importlib.import_module(name)
    assert hasattr(module, "__all__"), f"{name} has no __all__"
    missing = [n for n in module.__all__ if not hasattr(module, n)]
    assert not missing, f"{name}.__all__ lists names it does not define: {missing}"
    assert len(set(module.__all__)) == len(module.__all__), f"{name}.__all__ repeats a name"


def test_the_rescue_audit_contract_lists_its_functions():
    """Like io.output, it lists the functions; its constants are imported by name. (The
    qc-flags docs test reads every upper-case string literal in the source as a flag the
    code emits, so the outcome constants' names cannot be listed without documenting them
    as flags.)"""
    from gbcms import rescue_audit

    assert sorted(rescue_audit.__all__) == [
        "format_rescue_audit",
        "rescue_outcome",
        "rescued_component",
    ]


def test_the_mfsd_report_lists_its_public_names():
    from gbcms.report import mfsd_report

    assert sorted(mfsd_report.__all__) == ["TOOLTIPS", "generate_mfsd_report"]
    # Helpers and imports stay out of the public list.
    assert not any(
        n.startswith("_") or n in ("json", "logging", "math") for n in mfsd_report.__all__
    )


def test_the_qc_summary_lists_its_classes_and_functions():
    """Classes and functions only, like the rescue audit: ``QC_FLAGS`` and the other
    upper-case constants would be read by the qc-flags docs test as flags the code emits."""
    from gbcms import qc_summary

    assert sorted(qc_summary.__all__) == [
        "Kind",
        "QcFlag",
        "QcSummary",
        "Tier",
        "log_qc_summary",
        "render_summary",
        "summarize_output",
    ]
