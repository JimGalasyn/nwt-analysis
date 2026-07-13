# nwt-analysis

> ## ⛔ RETIRED — July 2026
>
> **The Null Worldtube Theory program is retired; its physical claims did not
> survive precision audit.** This repository is the papers' reproduction code —
> it still runs, but the numbers it reproduces are postdictions, not confirmed
> predictions. The 22 papers it reproduces are marked *obsolete* on Zenodo, and
> the physics library it builds on ([nwt-substrate](https://github.com/JimGalasyn/nwt-substrate))
> is deprecated. The one surviving maintained artifact of the program is
> [jax-solitons](https://github.com/JimGalasyn/jax-solitons).
>
> Full accounting: **[Retirement Retrospective](https://doi.org/10.5281/zenodo.21339662)**
> (Zenodo, DOI [10.5281/zenodo.21339662](https://doi.org/10.5281/zenodo.21339662)).
> This repo is archived read-only as part of the record.

Reproduction code for the **Null Worldtube Theory** papers. Each module under
`src/nwt_analysis/paperNN_*/` is a driver that computes a published paper's
figures and numbers, built on top of the NWT field-theory upstreams
[`nwt-substrate`](https://github.com/JimGalasyn/nwt-substrate) (physics library)
and [`jax-solitons`](https://github.com/JimGalasyn/jax-solitons) (JAX/GPU soliton
engine).

This package is the **single, citable home** for the papers' reproduction code.
Papers cite it by `pip install` + a Zenodo software DOI rather than linking loose
scripts — so a reader can reproduce a result without chasing files across repos.
See [`REPRODUCTIONS.md`](REPRODUCTIONS.md) for the full script → paper map.

## Scope and rigor (read before assuming uniform validation)

This repo holds **two kinds of code**, and they are **not** held to the same bar:

- **(a) Reproduction drivers over the engines** — configurations and analyses
  built on `nwt-substrate` / `jax-solitons`. These inherit the upstreams'
  validation (jax-solitons' exactness gates, nwt-substrate's cross-engine
  oracle); a green driver means the published number reproduces against a
  validated engine.
- **(b) Standalone numerical methods outside the engines' scope-by-design** —
  e.g. the BPS-Helmholtz κ eigenvalue solver lives here because
  `jax-solitons` is, by design (its DESIGN.md P3), an energy-functional-on-a-
  periodic-lattice engine and does **not** do eigenvalue problems on closed
  curved manifolds. These carry their own, independent validation.

Do not assume a result in (b) inherits an engine's rigor. Each module's
docstring/README states its validation status; treat per-module status as the
source of truth, not the repo as a whole.

## Install

```bash
pip install nwt-analysis            # core (numpy/scipy/matplotlib) drivers
pip install "nwt-analysis[engine]"  # + nwt-substrate & jax-solitons (field-theory drivers)
```

Until the upstreams publish to PyPI, install them from git (DOI-pinned tags):

```bash
pip install -e . -r requirements-engine.txt
```

Drivers that need the engine import it lazily, so the core install works on its own.

## Use

```bash
nwt-repro list                  # every reproduction, grouped by paper
nwt-repro list --paper 20       # just Paper 20 (neutrino sector)
nwt-repro run k7_pmns_mass_hierarchy_test -- --help
```

`nwt-repro run <name>` executes the driver exactly as it would standalone. Bare
names work when unambiguous; otherwise qualify as `paperNN/<name>`.

## Layout

```
src/nwt_analysis/
  paper06_form_factor/      paper17_qec/           paper20_neutrinos/
  paper07_knot_invariants/  paper18_gravity/       paper21_boson_rep/
  paper08_ewk_couplings/    paper19_substrate_figures/  paper22_baryogenesis/
  paper09_lepton_observables/  paper13_spectrum/
  supporting/               # shared / exploratory analysis, not one paper's reproduction
  cli.py                    # the nwt-repro dispatcher
```

`paperNN_*/` modules hold each paper's reproduction drivers (discoverable via
`nwt-repro`). `supporting/` holds shared and exploratory analysis scripts that
aren't tied to a single paper's reproduction (consolidated here from the old
`null-worldtube/analysis` directory). See [`REPRODUCTIONS.md`](REPRODUCTIONS.md).

## Relationship to the other repos

| repo | role |
|---|---|
| `nwt-substrate` | canonical physics library (primitives) |
| `jax-solitons` | JAX/GPU soliton engine |
| **`nwt-analysis`** | **per-paper reproduction drivers (this repo)** |
| `null-worldtube` | public paper site (links DOIs; hosts no code) |

The cross-repo `reconcile.py` gate asserts every script a *published* paper cites
is present here, so the data-availability links can never silently rot.

## License

MIT — see [LICENSE](LICENSE).
