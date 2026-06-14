"""α_KL disk-radius derivation — PHASE B: test the Sugawara/coadjoint route.

Per research_plan_alpha_KL_radius, Phase B's primary route is Sugawara/coadjoint:
the inflaton lives in K_8 = the so(7) spinor module, so the disk-radius numerator
is read as a spinor quantity (the "most-reduced reading" stated in
cosmesis ... alpha_KL_radius.py was dim(spinor)/2), normalized by the Sugawara
denominator h_v.  Phase A fixed the acceptance gate: a real derivation must
reproduce the WHOLE cross-substrate family AND respect the bare-NHEK/Starobinsky
boundary (α_KL=1 for algebras with h_v=rank+1).

This script RUNS that route through the gate.  Result: it FAILS — a clean,
computed negative result — and the failure resolves an ambiguity that Phase A
glossed over.  No derivation is faked.

  FINDING 1 — the spinor route fails the boundary.  Numerator = dim(spinor)/2
    gives √α_KL = (dim spinor/2)/h_v, hence B_2 → 2/3 (NOT the Starobinsky 1)
    and B_4 → 8/7 > 1 (α_KL>1: wrong direction — less curved than Starobinsky).
    The bare Kerr-NHEK ratio is 1 for ANY algebra, so an algebra whose dressing
    should be trivial cannot land at 2/3.  ⟹ Route 1 (spinor coadjoint) RULED OUT.

  FINDING 2 — the numerator is ambiguous off B_3, and the boundary picks it.
    rank+1, dim(spinor)/2, ⟨ρ,θ⟩ all equal 4 ONLY at B_3 (the triple coincidence);
    they diverge elsewhere.  Only (rank+1)/h_v → 1 exactly when h_v=rank+1 (the
    A_n/C_n/B_2 class), i.e. only the AFFINE-NODE-COUNT numerator is consistent
    with the bare-NHEK=Starobinsky boundary.  ⟹ the "K_8 = spinor" reading of the
    numerator is a B_3 coincidence; the boundary-forced numerator is rank+1 =
    (# affine Dynkin nodes).  This CORRECTS the radius-script's reading.

  CONSEQUENCE — the radius numerator is (affine rank)/(dual Coxeter), which is
    NOT a spinor-orbit quantity (Route 1 dead) and not a naive KK-reduction
    Casimir/dimension (Route 2 would give those).  The mechanism must be
    affine-structural.  Redirect: Route 3 (JT/Schwarzian) or an explicit
    affine-WZW argument for why (affine rank)/h_v fixes the disk radius.
"""
from __future__ import annotations
from fractions import Fraction as F


def banner(msg, ch="="):
    print(ch * 78)
    print(msg)
    print(ch * 78)


def so_data(n):
    """B_n = so(2n+1): rank, dual Coxeter, 2ρ, and weight Casimirs (long root²=2)."""
    rank, h_v = n, 2 * n - 1
    rho = [F(2 * n - 2 * i + 1, 2) for i in range(1, n + 1)]
    ip = lambda a, b: sum(x * y for x, y in zip(a, b))
    C2 = lambda lam: ip(lam, [lam[i] + 2 * rho[i] for i in range(n)])
    spin = [F(1, 2)] * n
    theta = [F(1), F(1)] + [F(0)] * (n - 2)
    return {
        "rank": rank, "h_v": h_v, "rank+1": rank + 1,
        "dim_spin_half": F(2 ** n, 2),
        "rho_theta": ip(rho, theta),
        "C2_spin": C2(spin),
    }


# ==========================================================================
# SECTION 1 — the spinor-route candidate radius
# ==========================================================================

def section1_spinor_route():
    banner("SECTION 1 — Route 1 candidate: inflaton in K_8 = so(7) spinor module")
    print()
    print("    K_8 IS the 8-dim spinor module of so(7); the inflaton (⟨Y_0Y_7⟩ Higgs)")
    print("    lives in it.  The natural spinor numerator that matches B_3 (=4) is")
    print("    dim(spinor)/2 = #positive spinor weights.  Sugawara denominator = h_v.")
    print()
    print("      candidate:  √α_KL = (dim spinor/2) / h_v    (spinor coadjoint/Sugawara)")
    print()
    d = so_data(3)
    print(f"    B_3 check:  dim(spinor)/2 = {d['dim_spin_half']}, h_v = {d['h_v']}  ⟹  √α_KL = {d['dim_spin_half']/d['h_v']}  ✓ (= 4/5)")
    print("    (B_3 works — but B_3 is exactly where the numerator is ambiguous.)")
    print()


# ==========================================================================
# SECTION 2 — run it through the family gate: FAILS
# ==========================================================================

def section2_gate_test():
    banner("SECTION 2 — gate test: spinor route vs the Starobinsky boundary → FAILS")
    print()
    print(f"    {'B_n':>4} {'√α_KL spinor':>14} {'√α_KL (rk+1)/hv':>16} {'boundary expects':>17}")
    fail = False
    for n in [2, 3, 4, 5]:
        d = so_data(n)
        spin_r = d["dim_spin_half"] / d["h_v"]
        node_r = F(d["rank+1"], d["h_v"])
        expect = "1 (Starobinsky)" if node_r == 1 else f"{node_r} (dressed)"
        flag = ""
        if node_r == 1 and spin_r != 1:
            flag = "  ✗ spinor≠1"; fail = True
        if spin_r > 1:
            flag += "  ✗ α_KL>1 wrong-dir"; fail = True
        print(f"    B_{n}  {str(spin_r):>13} {str(node_r):>16} {expect:>17}{flag}")
    print()
    print("  ★ FINDING 1: the spinor route FAILS the gate.")
    print("    • B_2: spinor gives 2/3, but bare-NHEK (algebra-independent) = 1, and the")
    print("      boundary class h_v=rank+1 must be undressed (Starobinsky).  2/3 ≠ 1.")
    print("    • B_4: spinor gives 8/7 ⟹ α_KL=(8/7)²>1, i.e. LESS curved than")
    print("      Starobinsky (r larger) — wrong direction; the substrate deviates the")
    print("      other way (B_3: 48 < 75).")
    print("    ⟹ Route 1 (Sugawara/coadjoint via the K_8 spinor) is RULED OUT.")
    print()
    return fail


# ==========================================================================
# SECTION 3 — resolve the numerator: the boundary forces affine-node-count
# ==========================================================================

def section3_numerator():
    banner("SECTION 3 — numerator resolution (affine nodes, not spinor)")
    print()
    print("    Three candidate numerators, all = 4 at B_3 (the triple coincidence),")
    print("    DIVERGE elsewhere:")
    print()
    print(f"    {'B_n':>4} {'rank+1':>7} {'dim(spin)/2':>12} {'⟨ρ,θ⟩':>7} {'(rk+1)/hv→1?':>13}")
    for n in [2, 3, 4, 5]:
        d = so_data(n)
        node1 = (F(d["rank+1"], d["h_v"]) == 1)
        print(f"    B_{n}  {d['rank+1']:>7} {str(d['dim_spin_half']):>12} {str(d['rho_theta']):>7}  {'yes' if node1 else 'no':>13}")
    print()
    print("    Only the AFFINE-NODE-COUNT numerator (rank+1) gives √α_KL→1 exactly when")
    print("    h_v=rank+1 (A_n/C_n/B_2) — the boundary-consistent class.  dim(spinor)/2")
    print("    and ⟨ρ,θ⟩ do not.  ⟹ the boundary FORCES numerator = rank+1 = #affine")
    print("    Dynkin nodes; the 'K_8 spinor' reading was a B_3 coincidence.")
    print()
    print("  ★ FINDING 2: CORRECTS the radius-script 'most-reduced reading'.  The")
    print("    numerator is the AFFINE NODE COUNT (rank+1), NOT dim(spinor)/2.  They")
    print("    agree only because B_3 has rank+1 = dim(spinor)/2 = ⟨ρ,θ⟩ = 4.")
    print()


# ==========================================================================
# SECTION 4 — redirect
# ==========================================================================

def section4_redirect():
    banner("SECTION 4 — consequence + redirect")
    print()
    print("    The boundary-forced radius is √α_KL = (rank+1)/h_v = (#affine Dynkin")
    print("    nodes)/(dual Coxeter).  This is NOT:")
    print("      • a spinor-orbit quantity  → Route 1 (Sugawara/coadjoint) DEAD;")
    print("      • a naive KK-reduction Casimir/dimension/index (those give the")
    print("        spinor-type ratios just ruled out) → naive Route 2 disfavored.")
    print()
    print("    The numerator (affine rank) is a DYNKIN-DIAGRAM combinatorial quantity.")
    print("    A genuine derivation must explain why (affine rank)/(dual Coxeter) sets")
    print("    the disk radius.  Candidate mechanisms to try next:")
    print("      • Route 3 — near-extremal JT/Schwarzian: the dilaton/Schwarzian")
    print("        coupling vs the so(7) current algebra normalization.")
    print("      • affine-WZW structural: (rank+1) = #level-0 simple roots + the")
    print("        affine root; ratio to h_v at level k=1 (k+h_v=h_Cox).")
    print()
    print("    The Phase C gate is now well-posed: the family is ((rank+1)/h_v)² with")
    print("    rank+1 = affine node count (boundary-fixed), and any derivation must")
    print("    reproduce THAT (A_n/C_n/B_2→1, B_3→16/25, B_4→25/49, ...).")
    print()


def verdict(failed):
    banner("PHASE B — verdict")
    print()
    print(f"  • Route 1 (Sugawara/coadjoint via K_8 spinor): RULED OUT.  Gate failed: {failed}")
    print("    spinor numerator dim(spinor)/2 gives B_2→2/3 (≠Starobinsky 1) and")
    print("    B_4→α_KL>1 (wrong direction).  Computed, not asserted.")
    print("  • Numerator RESOLVED: boundary forces rank+1 = #affine Dynkin nodes; the")
    print("    'K_8 = spinor' reading (dim spinor/2) is a B_3 coincidence — radius")
    print("    script reading corrected.")
    print("  • Open question REDIRECTED + SHARPENED: derive why (affine rank)/(dual")
    print("    Coxeter) sets the disk radius — affine-structural, not spinor/Casimir.")
    print("    Next: Route 3 (JT/Schwarzian) or explicit affine-WZW. Phase C gate now")
    print("    well-posed (affine-node family).")
    print()
    print("  HONEST: Phase B is a rigorous NEGATIVE result + a resolution, not a")
    print("  derivation of 4/5.  The primary route is eliminated for a computed reason.")
    print()


def main():
    print()
    banner("α_KL disk-radius — PHASE B (Sugawara/coadjoint route → gate test)", ch="█")
    print()
    section1_spinor_route()
    failed = section2_gate_test()
    section3_numerator()
    section4_redirect()
    verdict(failed)


if __name__ == "__main__":
    main()
