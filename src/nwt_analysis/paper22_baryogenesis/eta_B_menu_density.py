#!/usr/bin/env python3
"""eta_B candidate-menu density -- the sweep exhibit calculation.

Reconstructs the 17-formula sweep behind eta_B = (3/14) alpha^4
(baryon_asymmetry_memo.md, section "Candidate sweep verification") and measures
the look-elsewhere volume of its candidate menu. Three parts:

  PART 1  reproduce the memo's named top matches (validates alpha + arithmetic).
  PART 2  rank (3/14) alpha^4 against the FULL menu at the true Planck target,
          and count how many menu candidates "agree" within the 2-sigma bar --
          the direct look-elsewhere at eta_B itself.
  PART 3  the Auditor's mode-5: what fraction of RANDOM targets the same menu
          fits inside a Planck-width bar (full vs minimal integer menu). The
          eta_B analogue of paper6_menu_density.py; the Auditor reported 54.9%
          full / 23.8% minimal (verdict 4032d04).

The sweep artifact it accompanies is a prose summary (only 3 of 17 candidates
are written out, no enumeration was pinned); this script is the missing
reproducible enumeration. Menu = (p/q)*alpha^n over the memo's stated primitive
integers, plus the memo's QED-loop prefactors. The exact 17 the author hand-
picked are not recoverable from the memo; this enumerates the FULL menu those
sets imply, which is the honest look-elsewhere denominator.

Result (deterministic, seed 0, substrate alpha):
  PART 1: the three named candidates reproduce the memo exactly -- (3/14)a^4
          -0.38%, (7/32)a^4 +1.69%, a^4/5 -7.02%.
  PART 2: (3/14)a^4 IS the single closest member of the 519-candidate menu
          (rank 1) -- the memo's "top match" claim holds on its own terms.
  PART 3: yet the SAME menu (n=4) lands inside a 2-sigma Planck bar for ~55-63%
          of random targets (full menu) / ~19% (minimal 6-integer menu),
          reproducing the Auditor's mode-5 (54.9% / 23.8%, verdict 4032d04).
          Being the best of a menu that fits half of everything is not evidence.

Deterministic (seed fixed). stdlib + numpy only.

  python3 eta_B_menu_density.py
"""
import argparse
from math import pi, sqrt, gcd
import numpy as np

# --- constants ---------------------------------------------------------------
ALPHA_SUB = 1.0 / (25 * pi * sqrt(3) + 1)     # substrate alpha = 1/(25 pi sqrt3 + 1)
ALPHA_CODATA = 7.2973525693e-3                 # CODATA 2018
ETA_B = 6.10e-10                               # Planck 2018 central (memo's target)
ETA_B_SIG = 0.10e-10                           # Planck 2018 1-sigma (memo)

# memo's substrate-primitive integer menu and QED-loop prefactors
PRIMS = [1, 2, 3, 4, 5, 7, 8, 14, 16, 21, 28, 32, 35]
PRIMS_MIN = [1, 2, 3, 4, 5, 7]                 # a 6-integer "minimal" set (memo's exact 6 unspecified)
QED = [("1/pi", 1 / pi), ("1/4pi", 1 / (4 * pi)),
       ("1/16pi", 1 / (16 * pi)), ("1/32pi", 1 / (32 * pi))]


def menu_values(alpha, prims, powers, include_qed=True):
    """All candidate values (p/q)*alpha^n (+ QED-loop prefactors) over the sets."""
    vals = []
    for n in powers:
        an = alpha ** n
        for p in prims:
            for q in prims:
                vals.append((p / q) * an)
        if include_qed:
            for _, pref in QED:
                vals.append(pref * an)
    return np.array(vals)


def rel(a, b):
    return (a - b) / b


# --- PART 1: reproduce the memo's named top matches --------------------------
def part1(alpha):
    print("PART 1 -- named candidates vs Planck eta_B = 6.10e-10 (memo's top 3)")
    named = [("(3/14) a^4", 3 / 14, 4), ("(7/32) a^4", 7 / 32, 4), ("a^4 / 5", 1 / 5, 4)]
    print(f"  {'candidate':<14}{'value':>13}{'deviation':>12}   memo")
    memo_dev = {"(3/14) a^4": "0.38%", "(7/32) a^4": "1.69%", "a^4 / 5": "7.02%"}
    for lab, c, n in named:
        v = c * alpha ** n
        print(f"  {lab:<14}{v:>13.4e}{rel(v, ETA_B)*100:>+11.2f}%   (memo {memo_dev[lab]})")
    print()


# --- PART 2: rank at the true target + direct look-elsewhere ------------------
def part2(alpha, nsigma):
    powers = [3, 4, 5]
    vals = menu_values(alpha, PRIMS, powers, include_qed=True)
    hw = nsigma * ETA_B_SIG / ETA_B                        # relative half-width of the bar
    devs = np.abs(rel(vals, ETA_B))
    order = np.sort(devs)
    target = 3 / 14 * alpha ** 4
    rank = int(np.sum(devs < abs(rel(target, ETA_B)))) + 1
    inside = int(np.sum(devs <= hw))
    print(f"PART 2 -- full menu (p/q)*a^n, n in {powers}, {len(vals)} candidates")
    print(f"  (3/14)a^4 deviation {rel(target, ETA_B)*100:+.2f}%  ->  rank {rank} of {len(vals)} by |deviation|")
    print(f"  closest menu candidate: {order[0]*100:.3f}%")
    print(f"  So (3/14)a^4 IS the best menu member (the memo's claim holds) -- but that")
    print(f"  is not evidence until PART 3 shows how easily the menu fits ANY target.")
    print()


# --- PART 3: the Auditor's mode-5 (random-target coverage) --------------------
def coverage(alpha, prims, powers, hw, lo, hi, n_targets, seed):
    rng = np.random.default_rng(seed)
    vals = menu_values(alpha, prims, powers, include_qed=True)
    logv = np.sort(np.log(vals))
    targets = np.exp(rng.uniform(np.log(lo), np.log(hi), n_targets))
    hit = 0
    for t in targets:
        i = np.searchsorted(logv, np.log(t))
        near = logv[max(0, i - 1):i + 2]                  # nearest menu neighbours in log space
        if np.any(np.abs(np.exp(near) - t) / t <= hw):
            hit += 1
    return hit / n_targets


def part3(alpha, nsigma, lo, hi, n_targets, seed):
    hw = nsigma * ETA_B_SIG / ETA_B
    powers = [4]                                          # fix the scale; vary the INTEGER menu (Auditor's full vs minimal)
    print(f"PART 3 -- mode-5: fraction of random targets with a menu hit inside the "
          f"{nsigma}-sigma bar (+/-{hw*100:.2f}%)")
    print(f"  targets: {n_targets} log-uniform draws over [{lo:.0e}, {hi:.0e}], seed {seed}, n={powers}")
    full = coverage(alpha, PRIMS, powers, hw, lo, hi, n_targets, seed)
    mini = coverage(alpha, PRIMS_MIN, powers, hw, lo, hi, n_targets, seed)
    print(f"  FULL menu    ({len(PRIMS)} integers): {full*100:5.1f}%   (Auditor mode-5: 54.9%)")
    print(f"  MINIMAL menu ({len(PRIMS_MIN)} integers): {mini*100:5.1f}%   (Auditor mode-5: 23.8%)")
    print("  robustness across target ranges (full menu):")
    for a, b in [(3e-11, 3e-9), (1e-10, 1e-9), (3e-10, 1.2e-9)]:
        c = coverage(alpha, PRIMS, powers, hw, a, b, n_targets, seed)
        print(f"    [{a:.0e}, {b:.0e}]: {c*100:5.1f}%")
    print()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--alpha", choices=["substrate", "codata"], default="substrate")
    ap.add_argument("--nsigma", type=float, default=2.0)
    ap.add_argument("--targets", type=int, default=200000)
    ap.add_argument("--lo", type=float, default=1e-10)
    ap.add_argument("--hi", type=float, default=1e-9)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    alpha = ALPHA_SUB if args.alpha == "substrate" else ALPHA_CODATA
    print(f"alpha = {alpha:.10e}  ({args.alpha})   eta_B target = {ETA_B:.2e} +/- {ETA_B_SIG:.0e} (Planck 2018)")
    print(f"      1-sigma relative = {ETA_B_SIG/ETA_B*100:.2f}%   ({args.nsigma}-sigma bar = +/-{args.nsigma*ETA_B_SIG/ETA_B*100:.2f}%)\n")
    part1(alpha)
    part2(alpha, args.nsigma)
    part3(alpha, args.nsigma, args.lo, args.hi, args.targets, args.seed)
    print("Reading: a ppm/sub-percent match to eta_B is the expected output of this "
          "menu, not evidence about nature -- the look-elsewhere volume is large.")


if __name__ == "__main__":
    main()
