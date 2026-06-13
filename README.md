# nwt-analysis

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
  cli.py                    # the nwt-repro dispatcher
```

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
