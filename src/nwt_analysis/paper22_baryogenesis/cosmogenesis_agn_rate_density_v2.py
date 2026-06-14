"""Cosmogenic AGN rate density — v2 (rigorous refinement of the open audit item).

v1 (`cosmogenesis_agn_rate_density.py`, committed b0da18d) gave an order-of-magnitude
bracket R ~ 1e-7..1e-5 Mpc^-3 Gyr^-1 with hand-picked astrophysical factors. This v2
sharpens the three soft levers and, crucially, recasts the deliverable around the ONE
quantity that is actually observable — the candidate-parent AGN population density —
rather than the (intrinsically unobservable) nucleation rate.

WHAT THE SUBSTRATE FIXES (zero free parameters — the SHARP part):
  • mass band  M* · [1/√2, √2] = [6.86e7, 1.37e8] M_⊙   (Landauer ledger, Phase E.2)
  • spin       a*_eq = 0.998052 = Thorne accretion-equilibrium = K_7-viability spin
  • width      δa*  = 1.02e-4   (Aretakis resonance, gauntlet Stage 2/6)
These pin WHERE cosmogenesis can happen. Astrophysics supplies HOW OFTEN (the rate
normalization) and carries the residual uncertainty.

v2 refinements:
  1. n_SMBH(band): integrate a real local SMBH mass function (Schechter-in-log,
     literature-typical params) over the band — replaces the flat Φ≈3e-3 guess.
  2. f_spin: treat the Thorne spin as a coherent-accretion ATTRACTOR (our band is in
     the coherent/high-spin regime), giving a physically-justified central + bracket.
  3. Deliverable: the OBSERVABLE candidate-parent density n_cand and survey counts —
     the falsifiable proxy. The rate density and cousins/parent follow from it.

Astrophysical inputs are labelled literature-typical / order-of-magnitude; none are
fabricated precision. Runs clean with numpy only.
"""
from __future__ import annotations
import numpy as np

# ---------------------------------------------------------------------------
# Substrate-fixed selection (the sharp, zero-free-parameter part)
# ---------------------------------------------------------------------------
alpha       = 1 / (25 * np.pi * np.sqrt(3) + 1)
M_star      = 9.66e7                       # M_⊙  (Landauer ledger, Phase E.1/E.2)
band_lo     = M_star / np.sqrt(2)
band_hi     = M_star * np.sqrt(2)
band_dex    = np.log10(band_hi / band_lo)  # = log10(2) ≈ 0.301
a_star_eq   = 0.998052
delta_a     = 1.02e-4


def banner(msg, ch="="):
    print(ch * 78); print(msg); print(ch * 78)


# ===========================================================================
# SECTION 1 — local SMBH mass function integrated over the band
# ===========================================================================
def schechter_log(logM, logMknee, phi_star, slope):
    """Schechter function per dex in log10(M):
       Φ(logM) = ln10 · φ* · x^(1+slope) · exp(-x),  x = M/Mknee   [Mpc^-3 dex^-1]."""
    x = 10.0 ** (logM - logMknee)
    return np.log(10) * phi_star * x ** (1 + slope) * np.exp(-x)


def section1_bhmf():
    banner("SECTION 1 — n_SMBH(band) from a real local BH mass function")
    print()
    # Literature-typical local SMBH MF (Shankar+2009 / Vika+2009 / Marconi+2004 range):
    # Schechter-in-log, knee near 10^8.4, φ* ~ few×10^-3 Mpc^-3 dex^-1, faint slope ~ -0.6.
    logMknee, phi_star, slope = 8.4, 3.0e-3, -0.6
    # integrate over the band in log-space
    logs = np.linspace(np.log10(band_lo), np.log10(band_hi), 400)
    phis = schechter_log(logs, logMknee, phi_star, slope)
    _trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))
    n_band = _trapz(phis, logs)            # Mpc^-3
    phi_at_band = schechter_log(np.log10(M_star), logMknee, phi_star, slope)
    print(f"  Schechter-in-log MF (literature-typical local SMBH MF):")
    print(f"     log(Mknee/M_⊙)={logMknee}, φ*={phi_star:.1e} Mpc⁻³dex⁻¹, slope={slope}")
    print(f"  Φ at M* (10⁸ M_⊙)        ≈ {phi_at_band:.2e} Mpc⁻³ dex⁻¹")
    print(f"  band = [{band_lo:.2e}, {band_hi:.2e}] M_⊙  (width {band_dex:.3f} dex)")
    print(f"  ∫ Φ d(logM) over band     = n_SMBH(band) ≈ {n_band:.2e} Mpc⁻³")
    print(f"  (v1 flat estimate was {3e-3*band_dex:.2e} Mpc⁻³ — consistent; systematic ~×3)")
    print()
    return n_band


# ===========================================================================
# SECTION 2 — duty cycle f_active and spin fraction f_spin
# ===========================================================================
def section2_duty_and_spin():
    banner("SECTION 2 — duty cycle f_active and spin-resonance fraction f_spin")
    print()
    # --- f_active: AGN duty cycle at ~10^8 M_⊙ (fraction radiating at appreciable λ_Edd)
    f_active = 0.05                        # literature 1-10%; mass/z dependent
    print(f"  f_active ≈ {f_active}  (AGN duty cycle at 10⁸ M_⊙; literature 1–10%).")
    print()
    # --- f_spin: the Thorne spin is a coherent-accretion ATTRACTOR ---------
    print("  f_spin — the spin is NOT a random draw.  Coherently-accreting BHs are")
    print("  DRIVEN to the Thorne equilibrium a*_eq≈0.998 (a fixed point) and STAY")
    print("  there while accretion continues.  Our band (10⁸ M_⊙) is in the coherent-")
    print("  accretion / high-spin regime (X-ray reflection: most ~10⁷–10⁸ AGN have")
    print("  a*>0.9; chaotic-accretion low-spin sets in only at higher mass/mergers).")
    print()
    print("  So the relevant fraction is: of band AGN at high spin, what fraction sit")
    print("  within ±δa* of a*_eq?  = 2·δa*/σ_eq, where σ_eq is the scatter of the")
    print("  per-object equilibrium spin around 0.998 (disk-model/radiative-efficiency).")
    print()
    print(f"  {'σ_eq (equil. scatter)':24s}{'f_spin':>10s}   regime")
    for sig, reg in [(1e-4, "MAD/sharp equil → saturates"),
                     (5e-4, "tight (substrate-favoured)"),
                     (3e-3, "moderate disk-model spread"),
                     (1e-2, "loose (conservative)")]:
        fs = min(1.0, 2 * delta_a / sig)
        print(f"  {sig:>24.0e}{fs:>10.3f}   {reg}")
    print()
    # central: coherent-accretion attractor → tighter than v1's 0.02
    sigma_eq_central = 3e-3
    f_spin_central = min(1.0, 2 * delta_a / sigma_eq_central)
    f_spin_lo = min(1.0, 2 * delta_a / 1e-2)
    f_spin_hi = min(1.0, 2 * delta_a / 5e-4)
    print(f"  Central (σ_eq~3e-3): f_spin ≈ {f_spin_central:.3f}")
    print(f"  Bracket: f_spin ∈ [{f_spin_lo:.3f}, {f_spin_hi:.3f}]  ← DOMINANT uncertainty.")
    print(f"  (v1 central was 0.02; the attractor argument supports a ~3× larger central.)")
    print()
    return f_active, f_spin_central, (f_spin_lo, f_spin_hi)


# ===========================================================================
# SECTION 3 — the OBSERVABLE: candidate-parent population density
# ===========================================================================
def section3_candidate_population(n_band, f_active, f_spin, f_spin_bracket):
    banner("SECTION 3 — candidate-parent AGN density  (THE observable proxy)")
    print()
    print("  The nucleation rate is not observable.  The candidate-PARENT population")
    print("  is: active AGN in the mass band at near-extremal spin.  This is countable")
    print("  by mass (M-σ / reverberation) × activity (X-ray/optical) × spin (Fe Kα).")
    print()
    n_cand = n_band * f_active * f_spin
    n_cand_lo = n_band * f_active * f_spin_bracket[0]
    n_cand_hi = n_band * f_active * f_spin_bracket[1]
    print(f"  n_cand = n_SMBH(band) × f_active × f_spin")
    print(f"         = {n_band:.1e} × {f_active} × {f_spin:.3f}")
    print(f"         ≈ {n_cand:.2e} Mpc⁻³     (bracket [{n_cand_lo:.1e}, {n_cand_hi:.1e}])")
    print()
    # counts in a local survey volume (z<0.1: comoving R≈420 Mpc, V≈3.1e8 Mpc^3)
    V_z01 = (4 / 3) * np.pi * 420.0 ** 3
    N_z01 = n_cand * V_z01
    print(f"  Local volume z<0.1: V ≈ {V_z01:.1e} Mpc³")
    print(f"  ⟹ candidate parents within z<0.1:  N ≈ {N_z01:.0f}"
          f"  (bracket {n_cand_lo*V_z01:.0f}–{n_cand_hi*V_z01:.0f})")
    print()
    print("  ★ FALSIFIABLE PROXY: the substrate predicts an O(10²–10⁴) population of")
    print("    active ~10⁸ M_⊙ AGN at a*≈0.998 within z<0.1.  A near-total ABSENCE of")
    print("    high-spin AGN in this narrow mass band (if spin surveys grow) would")
    print("    disfavour the coherent-accretion-attractor picture the rate relies on.")
    print()
    return n_cand, (n_cand_lo, n_cand_hi)


# ===========================================================================
# SECTION 4 — rate density and cousins per parent
# ===========================================================================
def section4_rate_and_cousins(n_cand, n_cand_bracket):
    banner("SECTION 4 — nucleation rate density and cousins per parent")
    print()
    Gamma_trigger = 1.0          # Gyr^-1 per viable parent (Aretakis ε_disk>0.1 spike)
    R = n_cand * Gamma_trigger
    R_lo = n_cand_bracket[0] * Gamma_trigger
    R_hi = n_cand_bracket[1] * Gamma_trigger
    print(f"  Γ_trigger ≈ {Gamma_trigger:.1f} Gyr⁻¹ per viable parent (ε_disk>0.1 spike).")
    print(f"  R_cosmo = n_cand × Γ_trigger ≈ {R:.2e} Mpc⁻³ Gyr⁻¹"
          f"  (bracket [{R_lo:.1e}, {R_hi:.1e}])")
    print()
    # total events in observable universe
    V_obs, t_cosmic = 4e15, 13.8           # Mpc^3, Gyr
    N_events = R * V_obs * t_cosmic
    print(f"  Observable universe V≈{V_obs:.0e} Mpc³ over {t_cosmic} Gyr:")
    print(f"  ⟹ total nucleation events ≈ {N_events:.1e}  (≳10¹⁰ daughters; COMMON in")
    print(f"     absolute terms, HIGHLY SELECTIVE per-BH — narrow band × δa*=1e-4).")
    print()
    # cousins per parent: conditioned on being a viable parent (already at resonance)
    for life in (1.0, 3.0):
        cousins = Gamma_trigger * life * 0.5    # 0.5 = non-spin gauntlet residual
        print(f"  cousins/parent ({life:.0f} Gyr high-spin life) ≈ {cousins:.1f}")
    print()
    print("  ★ O(0.5–2) daughters per viable parent → genealogy is SHALLOW.  Cosmogenic-")
    print("    AGN population rank stays HIGH (~22–38, not collapsed to ~9): deep-")
    print("    genealogy does NOT rescue the 𝕆⊗𝕆-vs-no-hair degeneracy (consistent")
    print("    with the cousin-clustering 'observationally subtle' conclusion).")
    print()
    return R, (R_lo, R_hi), N_events


# ===========================================================================
# SECTION 5 — uncertainty budget + closure
# ===========================================================================
def section5_closure(R, R_bracket, N_events):
    banner("SECTION 5 — uncertainty budget + closure")
    print()
    print("  Uncertainty budget (multiplicative, on R and n_cand):")
    print(f"   • f_spin       : ×{1e-2/5e-4:.0f}   DOMINANT (σ_eq scatter around a*_eq)")
    print(f"   • n_SMBH(band) : ×3    (local BHMF normalisation/slope at 10⁸ M_⊙)")
    print(f"   • Γ_trigger    : ×3    (ε_disk>0.1 transient rate; ZTF/LSST can constrain)")
    print(f"   • f_active     : ×2    (duty cycle 1–10%)")
    print()
    print("  What is SUBSTRATE-FIXED (no uncertainty): the SELECTION — mass band")
    print("  M*·[1/√2,√2] and spin a*_eq=0.998052±1.02e-4. What carries the spread:")
    print("  the rate NORMALISATION (pure astrophysics). This cleanly separates the")
    print("  zero-parameter prediction (WHERE) from the astrophysical rate (HOW OFTEN).")
    print()
    print("  Testable NOW (coarse): band AGN should be high-spin (a*>0.99) — consistent")
    print("  with current X-ray reflection data. NOT testable yet (sharp): a* within")
    print("  1e-4 of 0.998 (spin precision is ~±0.05). DESI/LSST (mass+activity) +")
    print("  Athena/HEX-P (spin) sharpen n_cand; ZTF/LSST AGN variability constrains Γ.")
    print()
    print(f"  ★★ AUDIT ITEM RESOLVED (refined):")
    print(f"     R_cosmo ≈ {R:.1e} Mpc⁻³ Gyr⁻¹  (bracket [{R_bracket[0]:.0e}, {R_bracket[1]:.0e}])")
    print(f"     ≈ {N_events:.0e} daughters in the observable universe over cosmic time.")
    print(f"     Observable proxy = candidate-parent density (Section 3); cousins/parent")
    print(f"     O(1) → shallow genealogy. Selection substrate-fixed; rate astrophysical.")
    print()


def main():
    print()
    banner("Cosmogenic AGN rate density — v2 (refined)", ch="█")
    print()
    n_band = section1_bhmf()
    f_active, f_spin, f_spin_bracket = section2_duty_and_spin()
    n_cand, n_cand_bracket = section3_candidate_population(
        n_band, f_active, f_spin, f_spin_bracket)
    R, R_bracket, N_events = section4_rate_and_cousins(n_cand, n_cand_bracket)
    section5_closure(R, R_bracket, N_events)


if __name__ == "__main__":
    main()
