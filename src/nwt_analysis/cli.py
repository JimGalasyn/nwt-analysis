"""`nwt-repro` — discover and run the per-paper reproduction drivers.

Each reproduction is a module under nwt_analysis.paperNN_*/. This dispatcher
avoids declaring 56 separate console scripts: it walks the installed package,
lists the reproductions, and runs one by key with runpy (so a module's
`if __name__ == "__main__":` block executes exactly as it would standalone).

    nwt-repro list                      # all reproductions, grouped by paper
    nwt-repro list --paper 20           # just Paper 20
    nwt-repro run k7_pmns_mass_hierarchy_test -- --help   # run one (args after --)
    nwt-repro run paper20/k7_triality_delta_CP            # disambiguate by paper
"""
from __future__ import annotations

import argparse
import importlib
import pkgutil
import runpy
import sys


def _discover() -> dict[str, str]:
    """stem -> fully-qualified module name, for every reproduction module."""
    import nwt_analysis
    found: dict[str, str] = {}
    for pkg in pkgutil.iter_modules(nwt_analysis.__path__):
        if not pkg.ispkg or not (pkg.name.startswith("paper") or pkg.name == "supporting"):
            continue
        sub = importlib.import_module(f"nwt_analysis.{pkg.name}")
        for mod in pkgutil.iter_modules(sub.__path__):
            if mod.name == "__init__":
                continue
            found[f"{pkg.name}/{mod.name}"] = f"nwt_analysis.{pkg.name}.{mod.name}"
    return found


def _resolve(key: str, found: dict[str, str]) -> str:
    if key in found:
        return found[key]
    # allow bare stem if unambiguous
    matches = [full for k, full in found.items() if k.split("/", 1)[1] == key]
    if len(matches) == 1:
        return matches[0]
    if not matches:
        sys.exit(f"no reproduction matches '{key}'. Try: nwt-repro list")
    sys.exit(f"'{key}' is ambiguous: {matches}. Qualify it as paperNN/<name>.")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="nwt-repro", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    pl = sub.add_parser("list", help="list reproductions")
    pl.add_argument("--paper", type=int, default=None)
    pr = sub.add_parser("run", help="run a reproduction by key")
    pr.add_argument("key")
    pr.add_argument("rest", nargs=argparse.REMAINDER,
                    help="args passed to the driver (after --)")
    a = ap.parse_args(argv)

    found = _discover()
    if a.cmd == "list":
        keys = sorted(found)
        if a.paper is not None:
            keys = [k for k in keys if k.startswith(f"paper{a.paper:02d}")]
        cur = None
        for k in keys:
            pkg = k.split("/", 1)[0]
            if pkg != cur:
                print(f"\n{pkg}"); cur = pkg
            print(f"  {k.split('/', 1)[1]}")
        print(f"\n{len(keys)} script(s) (paperNN_* = reproductions; supporting/ = shared analysis).")
        return 0

    mod = _resolve(a.key, found)
    rest = a.rest[1:] if a.rest and a.rest[0] == "--" else a.rest
    sys.argv = [mod] + rest
    runpy.run_module(mod, run_name="__main__")
    return 0


if __name__ == "__main__":
    sys.exit(main())
