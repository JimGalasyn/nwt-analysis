"""Bit-rot guard for the reproduction drivers.

These tests do NOT execute the (often heavy) drivers — they assert that every
reproduction module *compiles* and is *discoverable* by the `nwt-repro`
dispatcher. That catches the failure mode this package exists to prevent: a
cited script that has silently broken or gone missing. Numerical / golden-value
tests per reproduction are added per-paper over time.
"""
from __future__ import annotations

import py_compile
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parent.parent / "src" / "nwt_analysis"
DRIVERS = sorted(p for p in SRC.rglob("paper*/*.py") if p.name != "__init__.py")


def test_some_drivers_present():
    assert len(DRIVERS) >= 50, f"only {len(DRIVERS)} drivers found under {SRC}"


@pytest.mark.parametrize("path", DRIVERS, ids=lambda p: f"{p.parent.name}/{p.name}")
def test_driver_compiles(path: Path):
    py_compile.compile(str(path), doraise=True)


def test_dispatcher_discovers_all():
    from nwt_analysis.cli import _discover
    found = _discover()
    assert len(found) == len(DRIVERS), (
        f"dispatcher found {len(found)} but {len(DRIVERS)} on disk")
