"""Bit-rot guard for the reproduction + supporting scripts.

These tests do NOT execute the (often heavy) scripts — they assert that every
script *compiles* and that every per-paper reproduction is *discoverable* by the
`nwt-repro` dispatcher. That catches the failure mode this package exists to
prevent: a cited script that has silently broken or gone missing. Numerical /
golden-value tests per reproduction are added per-paper over time.

Layout: `paperNN_*/` holds each paper's reproduction drivers (discoverable via
`nwt-repro`); `supporting/` holds shared / exploratory analysis scripts that are
not a single paper's reproduction (compiled here, but not dispatcher entries).
"""
from __future__ import annotations

import py_compile
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parent.parent / "src" / "nwt_analysis"
ALL_SCRIPTS = sorted(p for p in SRC.rglob("*/*.py") if p.name != "__init__.py")
PAPER_DRIVERS = sorted(p for p in SRC.rglob("paper*/*.py") if p.name != "__init__.py")


def test_scripts_present():
    assert len(ALL_SCRIPTS) >= 140, f"only {len(ALL_SCRIPTS)} scripts under {SRC}"


@pytest.mark.parametrize("path", ALL_SCRIPTS, ids=lambda p: f"{p.parent.name}/{p.name}")
def test_script_compiles(path: Path):
    py_compile.compile(str(path), doraise=True)


def test_dispatcher_discovers_all_scripts():
    # The dispatcher resolves both paper reproductions and supporting scripts so
    # `nwt-repro run <name>` works for any script cited in a paper.
    from nwt_analysis.cli import _discover
    found = _discover()
    assert len(found) == len(ALL_SCRIPTS), (
        f"dispatcher found {len(found)} but {len(ALL_SCRIPTS)} scripts on disk")
