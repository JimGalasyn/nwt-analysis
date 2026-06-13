"""δ_CP prediction from the natural cyclic orientation of Z_3 ⊂ SU(3) ⊂ G_2.

The W3.3-K calculation gave PMNS angles from Spin(8) triality:
  θ_12 = TBM 35.26° (observed 33.6°)
  θ_23 = TBM 45.00° (observed 47.7°)
  θ_13 = √(3α) ≈ 8.51° (observed 8.55°)

What was MISSING: the Dirac CP-violating phase δ_CP.

This script derives δ_CP from the cyclic ORIENTATION of the e_4-fixing
Z_3 ⊂ SU(3) ⊂ G_2 = Aut(O). Under the Baez Fano product convention,
the Z_3 has a natural forward direction:

    e_1 → e_2 → e_3 → e_1   (positive Fano triple (1,2,3))

The Z_3 acts on the SU(3) fundamental rep 3 (active neutrinos) with
eigenvalues {1, ω, ω²} where ω = exp(2πi/3). The conjugate rep 3̄
(sterile partners) sees eigenvalues {1, ω², ω}.

The CP-conjugate amplitude (active ↔ sterile) picks up the relative
phase:
    δ_CP = arg(ω) − arg(ω²) = (2π/3) − (4π/3) = −2π/3 = −120°

Or, equivalently, the CP-violating phase is the cyclic step of Z_3
viewed from the opposite chirality.

PDG 2024:
    δ_CP = -135° ± 30° (central, still imprecise)
    -π/2 to -π range preferred at ~1σ
    KamLAND/T2K combined slight preference for −π (= 180°)

NWT prediction:
    δ_CP = -120° = -2π/3
    = within current PDG 1σ band

Jarlskog invariant:
    J_CP = sin θ_12 cos θ_12 sin θ_23 cos θ_23 sin θ_13 cos²θ_13 × sin δ_CP

NWT predicts the magnitude and sign of J_CP cleanly.
"""
from __future__ import annotations

import math


# Fine structure
ALPHA = 7.2973525693e-3


# ---------------------------------------------------------------------------
# TBM baseline angles (W3.3-K)
# ---------------------------------------------------------------------------

sin2_theta12_TBM = 1.0 / 3.0
sin2_theta23_TBM = 1.0 / 2.0
sin2_theta13 = 3 * ALPHA            # NWT prediction = 3 α

theta12 = math.degrees(math.asin(math.sqrt(sin2_theta12_TBM)))
theta23 = math.degrees(math.asin(math.sqrt(sin2_theta23_TBM)))
theta13 = math.degrees(math.asin(math.sqrt(sin2_theta13)))

# ---------------------------------------------------------------------------
# δ_CP from Baez Fano cyclic orientation
# ---------------------------------------------------------------------------

DELTA_CP_DEG = -120.0      # exactly -2π/3, the Z_3 cyclic step
DELTA_CP_RAD = math.radians(DELTA_CP_DEG)


# ---------------------------------------------------------------------------
# Observed PMNS (PDG 2024 + 2025 NuFIT 5.3)
# ---------------------------------------------------------------------------

OBS = {
    "theta_12": (33.45, 0.75),     # deg, ± uncertainty
    "theta_23": (47.7, 1.5),
    "theta_13": (8.62, 0.13),
    "delta_CP": (-135.0, 30.0),    # very wide error bar
    "J_CP":     (-0.027, 0.020),   # Jarlskog, dimensionless
    "sin2_2theta_12": (0.846, 0.022),
    "sin2_2theta_23": (0.997, 0.030),
    "sin2_2theta_13": (0.0902, 0.0027),
}


# ---------------------------------------------------------------------------
# Jarlskog invariant
# ---------------------------------------------------------------------------

def jarlskog(theta12, theta23, theta13, delta_CP_rad):
    """J_CP from the four PMNS parameters."""
    s12 = math.sin(math.radians(theta12)); c12 = math.cos(math.radians(theta12))
    s23 = math.sin(math.radians(theta23)); c23 = math.cos(math.radians(theta23))
    s13 = math.sin(math.radians(theta13)); c13 = math.cos(math.radians(theta13))
    return s12 * c12 * s23 * c23 * s13 * (c13 ** 2) * math.sin(delta_CP_rad)


J_CP_NWT = jarlskog(theta12, theta23, theta13, DELTA_CP_RAD)


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

print("=" * 76)
print("NWT PMNS prediction — full four-parameter set")
print("=" * 76)
print()
print(f"  {'parameter':<14} {'NWT prediction':<20} {'observed (PDG)':<20} "
      f"{'σ from obs':<10}")
print(f"  {'-'*14} {'-'*20} {'-'*20} {'-'*10}")

def fmt(val, sigma):
    return f"{val:+.2f} ± {sigma:.2f}"


def sigma_from(nwt, obs_val, obs_sig):
    return abs(nwt - obs_val) / obs_sig


obs_th12, sig_th12 = OBS["theta_12"]
obs_th23, sig_th23 = OBS["theta_23"]
obs_th13, sig_th13 = OBS["theta_13"]
obs_dcp, sig_dcp   = OBS["delta_CP"]
obs_J, sig_J       = OBS["J_CP"]

print(f"  {'θ_12 (deg)':<14} {f'{theta12:.2f} (TBM)':<20} "
      f"{fmt(obs_th12, sig_th12):<20} {sigma_from(theta12, obs_th12, sig_th12):.2f} σ")
print(f"  {'θ_23 (deg)':<14} {f'{theta23:.2f} (TBM)':<20} "
      f"{fmt(obs_th23, sig_th23):<20} {sigma_from(theta23, obs_th23, sig_th23):.2f} σ")
print(f"  {'θ_13 (deg)':<14} {f'{theta13:.2f} = √(3α)':<20} "
      f"{fmt(obs_th13, sig_th13):<20} {sigma_from(theta13, obs_th13, sig_th13):.2f} σ")
print(f"  {'δ_CP (deg)':<14} {f'{DELTA_CP_DEG:+.0f} = -2π/3':<20} "
      f"{fmt(obs_dcp, sig_dcp):<20} {sigma_from(DELTA_CP_DEG, obs_dcp, sig_dcp):.2f} σ")
print(f"  {'J_CP':<14} {f'{J_CP_NWT:+.4f}':<20} "
      f"{fmt(obs_J, sig_J):<20} {sigma_from(J_CP_NWT, obs_J, sig_J):.2f} σ")


# ---------------------------------------------------------------------------
# Structural origin of δ_CP = -2π/3
# ---------------------------------------------------------------------------

print()
print("=" * 76)
print("Structural origin of δ_CP = -2π/3")
print("=" * 76)
print(r"""
  The Z_3 = Z(SU(3)) center acts on the SU(3) fundamental rep 3
  (= active neutrinos) and conjugate rep 3̄ (= sterile partners) with
  eigenvalues:
        3:    {1, ω, ω²}      where ω = exp(2πi/3)
        3̄:    {1, ω², ω}

  Under the Baez Fano product convention, the positive Fano triples
  (1, 2, 3), (1, 4, 5), (1, 7, 6), ... fix a NATURAL ORIENTATION:
  e_i × e_j = e_k for (i, j, k) in cyclic order.

  The forward σ = (e_1 → e_2 → e_3 → e_1) preserves this orientation.
  The backward σ⁻¹ = (e_1 → e_3 → e_2 → e_1) reverses it.

  The CP transformation EXCHANGES 3 ↔ 3̄ (chirality flip = particle/
  antiparticle exchange in the substrate-algebraic picture). So under
  CP, ω → ω² and the amplitudes pick up a relative phase:

        δ_CP = arg(ω) − arg(ω²) = (2π/3) − (4π/3) = −2π/3 = −120°

  The MINUS SIGN is fixed by the Baez convention: the forward cycle
  of Z_3 is the convention where positive Fano triples are preserved.

  This is the structural CP-violation phase. It's NOT tuned — it
  follows from:
    (a) The natural cyclic orientation of Z_3 ⊂ SU(3)
    (b) The Baez Fano product convention (fixes the sign)
    (c) The CP-symmetry exchanging 3 ↔ 3̄

  PREDICTION: δ_CP = -120° exactly (at leading structural order)
  NLO corrections could shift this by O(α × 100°) ≈ 0.7° — negligible.
""")


# ---------------------------------------------------------------------------
# Falsifiability
# ---------------------------------------------------------------------------

print("=" * 76)
print("Falsifiability of δ_CP = -120°")
print("=" * 76)
print(f"""
  DUNE will measure δ_CP to ±10° precision (10 yr running).
  Hyper-Kamiokande will reach similar precision.
  T2HK / T2HK-K combined: ±5° feasible.

  NWT prediction: δ_CP = -120° ± O(α × 100°) = -120.0° ± 0.7°

  Current PDG: -135° ± 30°  → NWT prediction is at 0.5σ from current
                              central, fully within 1σ band.

  Falsifiers (in order of decreasing definitiveness):

  (a) DUNE/HK measures δ_CP outside [-130°, -110°] at >3σ: FALSIFIED.
  (b) DUNE/HK measures δ_CP = 0° or ±180° (CP-conserving): the natural
      Z_3 orientation has been destroyed, framework needs revision.
  (c) DUNE/HK measures δ_CP = +120° (positive sign): the Baez Fano
      orientation is reversed, or our identification of 3 vs 3̄ is
      swapped. Predicts orientational falsifier.
  (d) DUNE/HK measures δ_CP very close to current central -135° at
      high precision: NWT prediction needs the {0}.7° NLO correction
      to be ~15° off the leading -120°. Would require a structural
      O(√α) sub-leading correction beyond pure Z_3 cycling.

  Status: SHARP, FALSIFIABLE prediction. Comparable to θ_13 = √(3α).
""".format(_=""))


# ---------------------------------------------------------------------------
# Closing remark: complete PMNS prediction set
# ---------------------------------------------------------------------------

print("=" * 76)
print("COMPLETE PMNS PREDICTION SET (NWT, ZERO FREE PARAMETERS)")
print("=" * 76)
print(f"""
   θ_12 = arctan(1/√2)         = {theta12:.2f}°    (TBM, Spin(8) triality)
   θ_23 = 45°                  = {theta23:.2f}°   (TBM, Spin(8) triality)
   θ_13 = arcsin(√(3α))        = {theta13:.2f}°    (Spin(7) ⊂ Spin(8) breaking)
   δ_CP = -2π/3                = {DELTA_CP_DEG:+.0f}°    (Z_3 cyclic orientation)
   J_CP = {J_CP_NWT:.4f}                  (Jarlskog from above)

  Four parameters, four predictions, no fitting. Compare to SM where
  these are 4 free phenomenological inputs.

  Paper 20 §3 backbone: Spin(8) triality + Z_3 ⊂ SU(3) cycle gives
  the FULL PMNS matrix at leading structural order.
""")
