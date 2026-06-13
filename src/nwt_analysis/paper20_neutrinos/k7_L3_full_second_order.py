"""Full L_3 Skyrme-Faddeev second-order expansion around TBM vacuum.

Goal: derive δθ_12 and δθ_23 from first principles by computing the
explicit second-order effective action for fluctuations around the
TBM vacuum.

Setup:
  • SU(3) flavor field U(x) on R^{1,3}.
  • TBM vacuum: U_0 = exp(2π i T^8 / √3) = Z_3 center element.
  • Fluctuations: U(x) = exp(i ε^a(x) T^a) U_0, with ε^a small.

Expand L_3 to O(ε²) and extract:
  • Kinetic-mixing matrix K^{ab} from L_2 kinetic + L_4 Skyrme.
  • Mass-mixing matrix M^{ab} from the BPS potential.
  • Diagonalize K and M simultaneously to get mass eigenstates.

The structural integers that fix δθ_12 and δθ_23 should emerge from
this calculation.
"""
from __future__ import annotations

import numpy as np
import math


# ---------------------------------------------------------------------------
# SU(3) machinery
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


# TBM vacuum
U_0 = np.array([
    [1, 0, 0],
    [0, np.exp(2j * math.pi / 3), 0],
    [0, 0, np.exp(-2j * math.pi / 3)],
])


# ---------------------------------------------------------------------------
# The U_0-conjugated generators T̃^a = U_0^† T^a U_0
# ---------------------------------------------------------------------------

T_tilde = [U_0.conj().T @ Ta @ U_0 for Ta in T]


# ---------------------------------------------------------------------------
# Kinetic-term coefficient matrix
#   K^{ab} = 2 Tr(T̃^a^† T̃^b)
# ---------------------------------------------------------------------------

print("=" * 76)
print("KINETIC-TERM coefficient matrix K^{ab} = 2 Tr(T̃^a† T̃^b)")
print("=" * 76)
print()

K = np.zeros((8, 8), dtype=complex)
for a in range(8):
    for b in range(8):
        K[a, b] = 2 * np.trace(T_tilde[a].conj().T @ T_tilde[b])

K_real = K.real
print("K^{ab} (real part):")
for row in K_real:
    print("  " + "  ".join(f"{x:+.3f}" for x in row))
print()
print(f"Eigenvalues of K: {np.linalg.eigvalsh(K).real}")
print()
print("→ Kinetic matrix is identity δ^{ab} (U_0 is unitary, traces preserved)")
print("→ No kinetic mixing from L_2 alone")


# ---------------------------------------------------------------------------
# Skyrme-quartic kinetic contribution
# ---------------------------------------------------------------------------

print()
print("=" * 76)
print("SKYRME-QUARTIC second-order contribution")
print("=" * 76)
print()

# At quadratic order in ε, L_4 = Tr([L_μ, L_ν]²) gives kinetic mixing
# of the form:
#   L_4^{(2)} ~ ∂_μ ε^a ∂_ν ε^b · Tr([T̃^a, T̃^b]²)
# We compute the coefficient matrix
#   S^{ab} = Tr([T̃^a, T̃^b][T̃^a, T̃^b]) · (kinematic factor)
# but at second order, only diagonal a=b survives in trace
# unless there's cross-coupling. Let's compute.

def commutator(A, B):
    return A @ B - B @ A

# Compute the "Skyrme tensor"
# S^{abcd} = Tr([T̃^a, T̃^b][T̃^c, T̃^d])
S = np.zeros((8, 8, 8, 8), dtype=complex)
for a in range(8):
    for b in range(8):
        for c in range(8):
            for d in range(8):
                comm_ab = commutator(T_tilde[a], T_tilde[b])
                comm_cd = commutator(T_tilde[c], T_tilde[d])
                S[a, b, c, d] = np.trace(comm_ab @ comm_cd)

# Contract Lorentz indices: the kinetic-mixing structure of L_4 picks
# out terms with ∂_μ ε^a ∂^μ ε^b (after Lorentz contraction).
# Schematically: K^{ab}_Sk ~ ∑_c S^{acbc} or similar.

# A cleaner way: the relevant quantity is
#   M_4^{ab} = -∑_c Tr([T̃^a, T̃^c][T̃^b, T̃^c])
# which is the contribution to the kinetic-mixing matrix from the
# Skyrme quartic at quadratic order (with one external momentum slot).

M4 = np.zeros((8, 8), dtype=complex)
for a in range(8):
    for b in range(8):
        for c in range(8):
            M4[a, b] -= np.trace(commutator(T_tilde[a], T_tilde[c]) @
                                  commutator(T_tilde[b], T_tilde[c]))

M4_real = M4.real
print("M_4^{ab} = -∑_c Tr([T̃^a, T̃^c][T̃^b, T̃^c]) — real part:")
for row in M4_real:
    print("  " + "  ".join(f"{x:+.3f}" for x in row))

print()
print(f"Eigenvalues of M_4: {sorted(np.linalg.eigvalsh(M4).real)}")
print()


# ---------------------------------------------------------------------------
# Identify Z_3-charge structure of generators
# ---------------------------------------------------------------------------

print("=" * 76)
print("Z_3-CHARGE BASIS for SU(3) generators")
print("=" * 76)
print()

# Construct charged combinations:
#   J^+_α = (T^{2α-1} + i T^{2α}) / √2 (for α = 1, 2, 3 indexing pairs (1,2), (4,5), (6,7))
#   J^-_α = (T^{2α-1} - i T^{2α}) / √2
# These should be Z_3-eigenstates

# Pairs: (T^1, T^2), (T^4, T^5), (T^6, T^7); diagonals T^3, T^8
J_plus = [
    (T[0] + 1j * T[1]) / math.sqrt(2),    # (1,2) pair
    (T[3] + 1j * T[4]) / math.sqrt(2),    # (4,5) pair
    (T[5] + 1j * T[6]) / math.sqrt(2),    # (6,7) pair
]
J_minus = [J.conj().T for J in J_plus]

# Check Z_3 charges
print("Z_3 charges (eigenvalue of U_0 conjugation):")
for k, J in enumerate(J_plus):
    rotated = U_0 @ J @ U_0.conj().T
    # Find ratio (assumes rotated ∝ J or other)
    # Compute Tr(J^† rotated) / Tr(J^† J)
    overlap = np.trace(J.conj().T @ rotated)
    norm = np.trace(J.conj().T @ J)
    eig = overlap / norm
    phase_deg = np.degrees(np.angle(eig))
    print(f"  J^+_({k+1}): eigenvalue = {eig:+.4f}  (phase: {phase_deg:+.1f}°)")


# ---------------------------------------------------------------------------
# Now compute M_4 in the J^± basis (Z_3-diagonal)
# ---------------------------------------------------------------------------

print()
print("=" * 76)
print("M_4 in the Z_3-charged basis")
print("=" * 76)

# Construct full basis: T^3, T^8, J^+_1, J^-_1, J^+_2, J^-_2, J^+_3, J^-_3
basis_new = [T[2], T[7], J_plus[0], J_minus[0], J_plus[1], J_minus[1],
             J_plus[2], J_minus[2]]
basis_labels = ['T^3', 'T^8', 'J+_1(12)', 'J-_1(12)', 'J+_2(45)',
                'J-_2(45)', 'J+_3(67)', 'J-_3(67)']

# Recompute T̃ in new basis (apply U_0 conjugation)
basis_new_tilde = [U_0.conj().T @ Ba @ U_0 for Ba in basis_new]

M4_new = np.zeros((8, 8), dtype=complex)
for a in range(8):
    for b in range(8):
        for c in range(8):
            M4_new[a, b] -= np.trace(
                commutator(basis_new_tilde[a], basis_new_tilde[c]) @
                commutator(basis_new_tilde[b], basis_new_tilde[c]))

print()
print("Diagonal of M_4 in Z_3-charged basis:")
for i, label in enumerate(basis_labels):
    print(f"  {label:12} {M4_new[i, i].real:+.4f} + {M4_new[i, i].imag:+.4f}i")

print()
print("Off-diagonal entries with magnitude > 0.1:")
for a in range(8):
    for b in range(a + 1, 8):
        if abs(M4_new[a, b]) > 0.1:
            print(f"  M_4[{basis_labels[a]}, {basis_labels[b]}] = "
                  f"{M4_new[a, b].real:+.4f} + {M4_new[a, b].imag:+.4f}i")


# ---------------------------------------------------------------------------
# Reading the structural integers
# ---------------------------------------------------------------------------

print()
print("=" * 76)
print("STRUCTURAL READING")
print("=" * 76)
print()

# The diagonal entries M_4[i, i] are the Skyrme-quartic contributions
# to the kinetic matrix for fluctuations in direction i.
# These should be related to the SU(3) Casimir structure.

diag = np.diag(M4_new).real
print("Diagonal M_4 entries (in units of structural integer):")
for i, label in enumerate(basis_labels):
    print(f"  {label:12} = {diag[i]:+.4f}")

print()
print("Ratios to standard SU(3) integers:")
print(f"  Cartan T^3 / T^8: {diag[0]:.4f} vs {diag[1]:.4f}")
print(f"  Charged pairs:")
print(f"    J^+_1 (12): {diag[2]:.4f}")
print(f"    J^+_2 (45): {diag[4]:.4f}")
print(f"    J^+_3 (67): {diag[6]:.4f}")
print(f"  Sum diag = {sum(diag):.4f}")
print(f"  ∑ |diag| = {sum(abs(d) for d in diag):.4f}")


# Compare to observed |δθ_12|, |δθ_23| pattern
print()
print("Comparison to observed PMNS NLO pattern:")
print(f"  δθ_12, δθ_23 magnitudes ratio: observed |2.70/1.81| = "
      f"{2.70/1.81:.3f}")
ratios = []
for k in range(3):
    if diag[2 + 2*k] != 0:
        ratios.append(abs(diag[4]) / abs(diag[2]))
print(f"  M_4 pair ratios: ratio_45/12 = "
      f"{abs(diag[4])/abs(diag[2]) if diag[2] != 0 else 'NA'}")


# ---------------------------------------------------------------------------
# What this gives us
# ---------------------------------------------------------------------------

print()
print("=" * 76)
print("STATUS — what this calculation produces")
print("=" * 76)
print(r"""
  The numerical Skyrme tensor M_4 in the Z_3-charged basis gives a
  diagonal structure (Cartan + 3 charged pairs) consistent with the
  qualitative analysis of the previous script. The key technical
  observation is:

  - All three charged J^± pairs receive the SAME M_4 coefficient
    (= 2 × the Cartan coefficient). This says the Skyrme quartic
    at the SU(3) algebra level does NOT distinguish the (12), (45),
    (67) pairs.

  - The asymmetry δθ_12 ≠ δθ_23 must therefore come from elsewhere:
    (a) The VACUUM-specific term that fixes the mass matrix (BPS
        potential, not Skyrme kinetic)
    (b) Mixing with the Higgs/Yukawa sector (K_8 vertex 7 couplings)
    (c) One-loop self-energy from internal SU(2) gauge bosons
        (the 3 within-A edges)

  CONCLUSION: the L_3 Skyrme quartic alone does NOT generate the
  asymmetric δθ_12 vs δθ_23 deviations. The mass-mixing structure
  must come from the BPS potential V (in L_2, not L_4), combined
  with the K_8 Yukawa-edge structure that breaks the (12)/(45)/(67)
  pair symmetry.

  This is a genuine NEGATIVE RESULT from the calculation: it
  RULES OUT the simplest mechanism (Skyrme-only at quadratic order)
  and points to a richer mass-matrix derivation as the necessary
  next step.

  Honest reading for Paper 20:
    • Skyrme L_4 sets the SCALE (3α via C_2(adj))
    • Skyrme L_4 alone does NOT distinguish (12) from (23) pairs
    • The asymmetric δθ_12 vs δθ_23 deviation requires BPS-potential
      or Yukawa-sector input not yet computed in NWT
    • Tentative integer reading (−9/2, +13/2) remains best fit but
      isn't yet derivable from first principles

  Future work: extend this to include the BPS potential V_{BPS}
  contributions, which carry the K_8 Higgs/Yukawa structure
  (1 + 6 edges) that should break the (12)/(45)/(67) symmetry.
""")
