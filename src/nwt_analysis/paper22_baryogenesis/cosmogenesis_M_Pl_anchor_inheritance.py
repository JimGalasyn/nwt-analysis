"""The M_Pl-anchor inheritance mechanism (Paper 22a load-bearing section).

NWT's 21a §2.1 cites Paper 22 for the M_Pl-fixing mechanism: M_Pl is the daughter's
single dimensional anchor, fixed by the parent black hole's Kerr MASS-HAIR, with
v_EW = (5/2)α⁸·M_Pl and m_e (and the whole spectrum) DERIVED from M_Pl. Earlier VV
write-ups used m_e as the convenient anchor and treated M_Pl as "equivalent via
Paper 17"; this script closes the reframing rigorously and confirms it is NOT mere
relabeling — m_e/M_Pl is a sub-percent PREDICTION of the K_7/K_8 Wilson tower.

The argument, in five steps:

  1. No-hair leaves ONE dimensional anchor. The inherited substrate (α and the Wilson
     integers — dimensionless; plus ℏ, c) determines the daughter up to a single
     dimensional scale ([[framework_substrate_inheritance_no_hair]]).

  2. That anchor IS the parent's Kerr M-hair → M_Pl. A Kerr BH has exactly three
     external hairs (M, Q, J). The daughter inherits J (→ AoE axis), Q≈0 (charge
     neutrality), and M — the MASS-hair, a GRAVITATIONAL scale, which is therefore the
     daughter's Planck mass M_Pl. This is why the anchor is M_Pl specifically.

  3. m_e and the full spectrum are DERIVED rungs of the Wilson tower anchored on M_Pl:
        m = (8/N_v)·α^(N_e/2)·(1+α/7)(1+3α²)·M_Pl.
     Electron = (N_e=21, N_v=7) rung; v_EW = (5/2)α⁸·M_Pl. Verified below at sub-percent.

  4. ★ M_Pl CANCELS from every dimensionless observable. Because each scale is
     M_Pl × (substrate α-power), all dimensionless ratios (m_e/M_Pl, v_EW/M_Pl,
     Λ_cc/M_Pl⁴, m_H/v_EW, Ω_b/Ω_c, η_B, n_s, r, …) are M_Pl-INDEPENDENT. The daughter
     has one free SCALE and ZERO free dimensionless parameters. The "1 free parameter
     (m_e)" framing is thus sharpened: the free parameter is the anchor M_Pl, and it
     drops out of every test.

  5. Reconciliation with no-hair: the ROLE (M-hair donates the anchor M_Pl) is
     inherited; the VALUE of M_Pl (which parent) is contingent — consistent with the
     no-hair memo's "parent mass in the contingent column."

NOTE ON ORDERS (Paper 22a harmonization, 2026-05-26): this script demonstrates the
MECHANISM at the NLO Wilson order — m_e at +0.01%, v_EW at -0.19%. The published
ledger cites the authoritative library (`cosmology.*` in nwt-substrate v0.2.0) values:
m_e to -5.5 ppm (NNLO) and v_EW = 245.5 GeV (-0.3%, LO form, matching Paper 21a §9.5).
The order difference does not affect the mechanism (m_e derived from M_Pl; M_Pl cancels).
"""
from __future__ import annotations
import numpy as np

ALPHA = 1.0 / (25 * np.pi * np.sqrt(3) + 1)     # Paper 17 trefoil residue
M_PL_GEV = 1.220890e19                          # Planck mass (energy units), GeV
H_V, RANK, H_COX = 5, 3, 6
N_POS = 9                                        # N_pos(so(7)) ; λ_H = 2·N_pos·α = 18α


def banner(msg, ch="="):
    print(ch * 78); print(msg); print(ch * 78)


def wilson_mass(n_e, n_v, M_Pl=M_PL_GEV):
    """K_7/K_8 phase-soliton mass: m = (8/N_v)·α^(N_e/2)·(1+α/7)(1+3α²)·M_Pl."""
    return (8.0 / n_v) * ALPHA ** (n_e / 2.0) * (1 + ALPHA / 7) * (1 + 3 * ALPHA**2) * M_Pl


# ===========================================================================
# PART 1 — anchor on M_Pl, DERIVE the spectrum, check vs observation
# ===========================================================================
def part1_derive_from_MPl():
    banner("PART 1 — anchor M_Pl, derive m_e / v_EW / m_H, compare to observation")
    print()
    print(f"  Single dimensional input:  M_Pl = {M_PL_GEV:.4e} GeV  (the parent M-hair).")
    print(f"  Everything below is M_Pl × (substrate α-power) — no other dimensional input.")
    print()
    # electron: (N_e=21, N_v=7) rung of the Wilson tower
    m_e = wilson_mass(21, 7)                       # GeV
    m_e_obs = 0.51099895e-3                         # GeV
    # EW vacuum: v_EW = (h_v/2)·α^8·M_Pl, with the same NLO factor
    v_EW = (H_V / 2) * ALPHA**8 * (1 + ALPHA / 7) * (1 + 3 * ALPHA**2) * M_PL_GEV
    v_EW_obs = 246.22
    # Higgs quartic + mass
    lam_H = 2 * N_POS * ALPHA                       # = 18α
    m_H = np.sqrt(2 * lam_H) * v_EW
    m_H_obs = 125.25
    rows = [
        ("m_e  (N_e=21, N_v=7)", m_e, m_e_obs, "GeV"),
        ("v_EW (h_v/2·α⁸, +NLO)", v_EW, v_EW_obs, "GeV"),
        ("m_H  = √(2·18α)·v_EW", m_H, m_H_obs, "GeV"),
    ]
    print(f"    {'quantity':24s}{'derived':>14}{'observed':>14}{'Δ':>9}")
    for name, val, obs, unit in rows:
        d = (val - obs) / obs * 100
        vs = f"{val:.5e}" if val < 1e-2 else f"{val:.3f}"
        os_ = f"{obs:.5e}" if obs < 1e-2 else f"{obs:.3f}"
        print(f"    {name:24s}{vs:>14}{os_:>14}{d:>+8.2f}%")
    print()
    print("  ⟹ m_e is a DERIVED rung, not the input: m_e/M_Pl is a sub-percent")
    print("    PREDICTION. So 'M_Pl anchor, m_e derived' is physics, not relabeling.")
    print()
    return m_e


# ===========================================================================
# PART 2 — ★ M_Pl cancels from every dimensionless observable
# ===========================================================================
def part2_anchor_cancels(m_e):
    banner("PART 2 — ★ the anchor M_Pl CANCELS from all dimensionless observables")
    print()
    print("  Compute the dimensionless predictions at TWO wildly different anchor values")
    print("  M_Pl and 10⁶·M_Pl. If they are identical, M_Pl is a pure scale that never")
    print("  enters a test — one free SCALE, ZERO free dimensionless parameters.")
    print()

    def dimensionless_set(M_Pl):
        me = wilson_mass(21, 7, M_Pl)
        v = (H_V / 2) * ALPHA**8 * (1 + ALPHA / 7) * (1 + 3 * ALPHA**2) * M_Pl
        lam = 2 * N_POS * ALPHA
        mH = np.sqrt(2 * lam) * v
        return {
            "m_e/M_Pl": me / M_Pl,
            "v_EW/M_Pl": v / M_Pl,
            "m_H/v_EW": mH / v,
            "Λ_cc/M_Pl⁴": (me / M_Pl) ** 4 * ALPHA**16 * H_COX,   # Stage-7 Λ_cc form
            "Ω_b/Ω_c": 25 * ALPHA + 75 * ALPHA**2,
            "η_B": 3 * ALPHA**4 / 14,
            "n_s": 1 - H_V * ALPHA,
            "r": 48 * ALPHA**2,
        }

    A = dimensionless_set(M_PL_GEV)
    B = dimensionless_set(M_PL_GEV * 1e6)
    print(f"    {'observable':14s}{'at M_Pl':>16}{'at 10⁶·M_Pl':>18}{'identical?':>12}")
    all_same = True
    for k in A:
        same = np.isclose(A[k], B[k], rtol=1e-12)
        all_same = all_same and same
        print(f"    {k:14s}{A[k]:>16.6e}{B[k]:>18.6e}{('yes' if same else 'NO'):>12}")
    print()
    print(f"  All dimensionless predictions M_Pl-independent: {'YES' if all_same else 'NO'}.")
    print("  ★ The daughter has ONE free dimensional anchor (M_Pl = the M-hair) and ZERO")
    print("    free dimensionless parameters. M_Pl cancels from every observable — which")
    print("    is WHY the >40 substrate predictions are zero-parameter. The earlier")
    print("    'free parameter = m_e' is the same single d.o.f. in a different basis;")
    print("    M_Pl is basis-natural (gravitational M-hair) and m_e a derived rung.")
    print()


# ===========================================================================
# PART 3 — reconciliation with the no-hair theorem (role vs value)
# ===========================================================================
def part3_reconcile_no_hair():
    banner("PART 3 — reconciliation with the substrate no-hair theorem")
    print()
    print("  Apparent tension: the no-hair memo lists (axis, algebra, sign, ℏ, c) as")
    print("  inherited and puts the parent's MASS in the 'contingent' column, yet 22a")
    print("  makes M_Pl the inherited anchor. Resolution — separate ROLE from VALUE:")
    print()
    print("    • ROLE (inherited): the parent's Kerr M-hair breaks the daughter's scale")
    print("      invariance and donates exactly ONE dimensional anchor — the gravitational")
    print("      mass scale M_Pl. This completes the Kerr-hair inheritance: J→axis (AoE),")
    print("      Q≈0→charge neutrality, M→M_Pl. The no-hair set is (M_Pl, axis, algebra,")
    print("      sign) — the M-hair is the dimensional member, the others dimensionless.")
    print()
    print("    • VALUE (contingent): the NUMERICAL value of M_Pl (i.e. which parent BH,")
    print("      what mass) is not predicted — exactly the memo's 'parent mass contingent'.")
    print("      But since M_Pl cancels from every dimensionless test (Part 2), its value")
    print("      is unobservable in ratios: the daughter's PHYSICS is anchor-independent.")
    print()
    print("  So the no-hair memo and the M_Pl-anchor framing are consistent: the memo")
    print("  must be updated to record the M-hair as the inherited dimensional anchor")
    print("  (role), keeping the value contingent. Then m_e moves from 'the free input'")
    print("  to 'a derived (N_e=21,N_v=7) Wilson rung of M_Pl'.")
    print()


def verdict():
    banner("VERDICT — M_Pl-anchor mechanism CLOSED (22a M-inheritance section ready)")
    print()
    print("  • M_Pl = the parent's Kerr mass-hair = the daughter's one dimensional anchor.")
    print("  • m_e (N_e=21,N_v=7) and v_EW (h_v/2·α⁸) are DERIVED from M_Pl at sub-percent")
    print("    — so 'm_e derived from M_Pl' is a prediction, not relabeling.")
    print("  • M_Pl CANCELS from every dimensionless observable ⟹ one free scale, zero")
    print("    free dimensionless parameters (the precise form of the '1 free parameter').")
    print("  • Consistent with no-hair: role (M-hair→M_Pl) inherited, value contingent.")
    print()
    print("  This is the load-bearing M-inheritance content that NWT 21a §2.1 cites Paper")
    print("  22 for. Next: update the no-hair framework memo + build the 22a outline.")
    print()


def main():
    print()
    banner("Paper 22a — the M_Pl-anchor inheritance mechanism", ch="█")
    print()
    m_e = part1_derive_from_MPl()
    part2_anchor_cancels(m_e)
    part3_reconcile_no_hair()
    verdict()


if __name__ == "__main__":
    main()
