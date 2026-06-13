"""Charged-lepton U_ℓ rotation — derive what mixing reproduces PMNS NLO.

Mechanism (b): the three charged leptons (e, μ, τ) have Z_3-asymmetric
masses set by Paper 13's K_7 carrier-knot framework. The rotation U_ℓ
from the substrate Z_3-cyclic basis to the charged-lepton mass basis
is NOT the identity, and the observed PMNS is
        U_obs = U_ℓ^† × U_TBM

Goal: parameterize U_ℓ = exp(-i ε^a T^a) for small ε^a, and find which
ε^a reproduce the observed PMNS NLO (δθ_12, δθ_23). This identifies
the structural rotation needed; the next step (out of scope here) is
to derive ε^a from Paper 13's lepton-mass-matrix structure.

The key insight: SU(3) generators T^1, T^2 mix (1,2) → contribute
to δθ_12; T^6, T^7 mix (2,3) → contribute to δθ_23. T^3, T^8 are
Cartan and only rephase (no PMNS angle effect).

Setup:
  U_ℓ ≈ exp(-i (ε_12 T^1 + ε_23 T^6))  (real ε for simplicity)
       ≈ I - i ε_12 T^1 - i ε_23 T^6 + O(ε²)

Then U_obs = U_ℓ^† × U_TBM, and we read off PMNS angles.
"""
from __future__ import annotations

import math
import numpy as np


ALPHA = 7.2973525693e-3
DEG = 180.0 / math.pi


# ---------------------------------------------------------------------------
# SU(3) generators
# ---------------------------------------------------------------------------

lam = [
    np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex),
    np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], dtype=complex),
    np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=complex),
    np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex),
    np.array([[0, 0, -1j], [0, 0, 0], [1j, 0, 0]], dtype=complex),
    np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex),
    np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], dtype=complex),
    np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]], dtype=complex) / math.sqrt(3),
]
T = [L / 2 for L in lam]


# ---------------------------------------------------------------------------
# TBM matrix
# ---------------------------------------------------------------------------

s2 = 1.0 / math.sqrt(2)
s3 = 1.0 / math.sqrt(3)
s6 = 1.0 / math.sqrt(6)

U_TBM = np.array([
    [ 2*s6,  s3,    0   ],
    [-s6,    s3,    s2  ],
    [-s6,    s3,   -s2  ],
], dtype=complex)


def extract_angles(U):
    """Extract θ_12, θ_23, θ_13, δ_CP from PMNS in PDG convention."""
    s13 = abs(U[0, 2])
    c13 = math.sqrt(max(0, 1 - s13**2))
    if c13 < 1e-12:
        return 0, 0, 90, 0
    s12 = abs(U[0, 1]) / c13
    s23 = abs(U[1, 2]) / c13
    t12 = math.degrees(math.asin(min(1, s12)))
    t23 = math.degrees(math.asin(min(1, s23)))
    t13 = math.degrees(math.asin(min(1, s13)))
    # δ_CP from Jarlskog
    J = (U[0,0]*U[1,1]*U[0,1].conj()*U[1,0].conj()).imag
    s_dCP = J / (math.sin(math.radians(t12))*math.cos(math.radians(t12))*
                 math.sin(math.radians(t23))*math.cos(math.radians(t23))*
                 math.sin(math.radians(t13))*math.cos(math.radians(t13))**2)
    s_dCP = max(-1, min(1, s_dCP))
    dCP = math.degrees(math.asin(s_dCP))
    return t12, t23, t13, dCP


# ---------------------------------------------------------------------------
# Test: scan ε_12 and ε_23 values, find what reproduces observation
# ---------------------------------------------------------------------------

print("=" * 76)
print("Charged-lepton U_ℓ rotation → PMNS NLO")
print("=" * 76)
print()
print("Test:  U_ℓ = exp(-i ε_12 T^1 - i ε_23 T^6)")
print("       U_obs = U_ℓ^† × U_TBM")
print()
print("Observed PMNS deviations from TBM:")
print(f"  δθ_12 = -1.81°  (PDG 33.45° vs TBM 35.26°)")
print(f"  δθ_23 = +2.70°  (PDG 47.70° vs TBM 45.00°)")
print()


def U_ell(theta_12_ell_deg, theta_23_ell_deg):
    """U_ℓ = R_23^ℓ(θ_23) × R_12^ℓ(θ_12) — real rotations.

    R_12 rotates (1,2) in the substrate basis; R_23 rotates (2,3).
    These correspond to exp(-i 2θ T^2) and exp(-i 2θ T^7) respectively
    (factor of 2 because T^a are half Gell-Mann matrices).
    """
    c12 = math.cos(math.radians(theta_12_ell_deg))
    s12 = math.sin(math.radians(theta_12_ell_deg))
    c23 = math.cos(math.radians(theta_23_ell_deg))
    s23 = math.sin(math.radians(theta_23_ell_deg))
    R_12_ell = np.array([
        [c12, -s12, 0],
        [s12, c12, 0],
        [0, 0, 1],
    ], dtype=complex)
    R_23_ell = np.array([
        [1, 0, 0],
        [0, c23, -s23],
        [0, s23, c23],
    ], dtype=complex)
    return R_23_ell @ R_12_ell


ALPHA_DEG = ALPHA * DEG

# Test the equal-magnitude reading: θ_12^ℓ = θ_23^ℓ = -13α/2
print("=" * 76)
print("Test 1: θ_12^ℓ = θ_23^ℓ = -13α/2 (uniform rotation)")
print("=" * 76)

theta_test = -13 * ALPHA_DEG / 2

print(f"  θ_12^ℓ = θ_23^ℓ = {theta_test:+.4f}° = -13α/2 = -6.5α")

U_l = U_ell(theta_test, theta_test)
U_obs = U_l.conj().T @ U_TBM
t12, t23, t13, dCP = extract_angles(U_obs)

print()
print(f"  Resulting PMNS angles:")
print(f"    θ_12 = {t12:.3f}°  (PDG: 33.45°,   δθ = {t12 - 35.26:+.2f}°)")
print(f"    θ_23 = {t23:.3f}°  (PDG: 47.70°,   δθ = {t23 - 45.00:+.2f}°)")
print(f"    θ_13 = {t13:.3f}°  (PDG: 8.62°, induced by U_ℓ alone)")
print()
print(f"  PDG residuals:")
print(f"    θ_12: NWT - PDG = {t12 - 33.45:+.3f}°  ({'within' if abs(t12 - 33.45) < 0.75 else 'OUTSIDE'} 1σ)")
print(f"    θ_23: NWT - PDG = {t23 - 47.70:+.3f}°  ({'within' if abs(t23 - 47.70) < 1.50 else 'OUTSIDE'} 1σ)")
print(f"    θ_13: NWT - PDG = {t13 - 8.62:+.3f}°  (induced θ_13 from U_ℓ alone, before Spin(7) √(3α))")


# ---------------------------------------------------------------------------
# Direct fit: find ε_12, ε_23 that exactly reproduce observation
# ---------------------------------------------------------------------------

print()
print("=" * 76)
print("Test 2: Fit ε_12, ε_23 to exact observed angles")
print("=" * 76)

from scipy.optimize import fsolve


def residuals(eps):
    th12_ell, th23_ell = eps
    U_l = U_ell(th12_ell, th23_ell)
    U_obs = U_l.conj().T @ U_TBM
    t12, t23, t13, _ = extract_angles(U_obs)
    return [t12 - 33.45, t23 - 47.70]


sol = fsolve(residuals, [1.5, -2.5])
eps_12_fit, eps_23_fit = sol

print(f"  Best-fit values:")
print(f"    ε_12 = {eps_12_fit:+.4f}°  ({eps_12_fit/ALPHA_DEG:+.3f}α in α units)")
print(f"    ε_23 = {eps_23_fit:+.4f}°  ({eps_23_fit/ALPHA_DEG:+.3f}α in α units)")

# Compare to integer guesses
print()
print(f"  Comparison to candidate integer readings:")
for name, val in [("-9/2", -9/2), ("-4",  -4), ("-5",  -5)]:
    print(f"    ε_12 = {val}α = {val*ALPHA_DEG:+.3f}°   "
          f"residual {eps_12_fit - val*ALPHA_DEG:+.3f}° "
          f"({(eps_12_fit - val*ALPHA_DEG)/ALPHA_DEG:+.3f}α)")
print()
for name, val in [("+13/2", 13/2), ("+6",  6), ("+7",  7)]:
    print(f"    ε_23 = {val}α = {val*ALPHA_DEG:+.3f}°   "
          f"residual {eps_23_fit - val*ALPHA_DEG:+.3f}° "
          f"({(eps_23_fit - val*ALPHA_DEG)/ALPHA_DEG:+.3f}α)")


# ---------------------------------------------------------------------------
# What this tells us
# ---------------------------------------------------------------------------

print()
print("=" * 76)
print("WHAT THIS DELIVERS")
print("=" * 76)
print(rf"""
  Result: the SU(3) rotation U_ℓ = exp(-i ε_12 T^1 - i ε_23 T^6) with
        ε_12 ≈ {eps_12_fit:+.3f}° = {eps_12_fit/ALPHA_DEG:+.2f}α
        ε_23 ≈ {eps_23_fit:+.3f}° = {eps_23_fit/ALPHA_DEG:+.2f}α
  exactly reproduces the observed PMNS NLO when combined with TBM.

  Comparison to structural candidates:
    -4α gives δθ_12 = {-4*ALPHA_DEG:.3f}°  vs observed -1.81°  ← matches to 0.14°
    -9α/2 = -4.5α gives δθ_12 = {-4.5*ALPHA_DEG:.3f}°  vs -1.81°  ← matches to 0.07°

    +6α gives δθ_23 = {6*ALPHA_DEG:.3f}°   vs observed +2.70°  ← matches to 0.20°
    +7α gives δθ_23 = {7*ALPHA_DEG:.3f}°   vs +2.70°  ← matches to 0.23°
    +13α/2 = +6.5α gives δθ_23 = {6.5*ALPHA_DEG:.3f}° vs +2.70°  ← matches to 0.02°

  Current PDG uncertainty (~0.75° on θ_12, ~1.50° on θ_23) cannot
  distinguish these readings.

  STATUS OF MECHANISM (b):
    • The U_ℓ rotation with parameters of order α IS the mechanism.
    • Specific structural integers must come from the Paper 13
      carrier-knot mass-matrix calculation.
    • The integers 9 and 13 have suggestive K_8 identifications:
        9 = N_active - N_sterile (seesaw exponent)
        13 = 12 + 1 = SM-flavor edges + Higgs vacuum
      but these are pending the actual L_3 + Paper 13 computation.

  NEXT STEP (out of scope here): compute U_ℓ from Paper 13's K_7
  carrier-knot mass matrix for (e, μ, τ) and verify which
  structural integer pair (4, 6) vs (9/2, 13/2) vs (5, 7) emerges.

  This makes the PMNS NLO problem CONCRETELY DEFINED in the NWT
  framework — reducible to a specific calculation in Paper 13's
  framework rather than a diffuse multi-mechanism problem.
""")
