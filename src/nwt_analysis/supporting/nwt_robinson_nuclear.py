#!/usr/bin/env python3
"""Robinson Charge-Overlap Nuclear Binding Model.

Robinson's nuclear model ("The Common Sense Universe", Ch. 5) treats nucleon
binding as electrostatic overlap of harmonic charge rings.  Each nucleon's
internal photon oscillation generates f/3, f/9, f/27 harmonics that create
concentric rings of alternating charge.  When nucleons sit side-by-side,
opposite-charge rings from neighbouring nucleons overlap → attraction → binding.

In NWT, Robinson's harmonics are the pure toroidal overtone modes on a k=3
torus (analyses.py Section 7).  This script tests whether the charge-overlap
mechanism produces correct binding energies.

Energy scale: k_e · e² / (1 fm) = 1.44 MeV — already nuclear.

Usage:
    python nwt_robinson_nuclear.py
"""
import os
import numpy as np
from scipy.special import ellipk
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ---- Physical constants ----
c0 = 2.99792458e8           # m/s
hbar = 1.054571817e-34       # J·s
e_charge = 1.602176634e-19   # C
m_p = 1.67262192369e-27      # kg
alpha = 1 / 137.035999084
epsilon_0 = 8.8541878128e-12
k_e = 1 / (4 * np.pi * epsilon_0)  # Coulomb constant
MeV = 1.602176634e-13        # J per MeV
fm = 1e-15                    # m per fm

hbar_c_MeV_fm = 197.3269804  # MeV·fm
k_e_MeV_fm = k_e * e_charge**2 / (fm * MeV)  # k_e·e²/(fm) in MeV ≈ 1.44

out_dir = os.path.dirname(os.path.abspath(__file__))

# ---- Robinson's nucleon harmonic data (p216b) ----
# Each harmonic: {name, n (harmonic number), E_MeV, Q_e (charge in units of e)}
# Radii: r_n = hbar_c / E_n  in fm

PROTON_HARMONICS = [
    {'name': 'P/3',  'n': 3,  'E_MeV': 312.9, 'Q_e': +2.0/3},
    {'name': 'P/9',  'n': 9,  'E_MeV': 104.3, 'Q_e': +1.0/3},
]

NEUTRON_HARMONICS = [
    {'name': 'N/3',  'n': 3,  'E_MeV': 313.2, 'Q_e': -1.0/3},
    {'name': 'N/9',  'n': 9,  'E_MeV': 104.5, 'Q_e': +2.0/3},
    {'name': 'N/27', 'n': 27, 'E_MeV':  34.8, 'Q_e': -1.0/3},
]

# Experimental binding energies (total, in MeV)
EXPERIMENTAL_BE = {
    '2H':   2.224,
    '3H':   8.482,
    '3He':  7.718,
    '4He': 28.296,
    '6Li': 31.995,
    '7Li': 39.245,
    '12C': 92.162,
    '16O': 127.619,
    '56Fe': 492.254,
}


# ---- 1. Charge profiles ----

def harmonic_ring_radius(E_MeV):
    """r_n = hbar_c / E_n in fm."""
    return hbar_c_MeV_fm / E_MeV


def nucleon_charge_density(r_fm, harmonics, sigma_fm=0.2):
    """Radial charge density rho(r) in units of e/fm.

    Sum of Gaussian rings centred at each harmonic radius:
        Q_n · exp(-(r - r_n)² / (2σ²)) / (√(2π) σ)
    """
    rho = np.zeros_like(np.atleast_1d(r_fm), dtype=float)
    for h in harmonics:
        r_n = harmonic_ring_radius(h['E_MeV'])
        Q = h['Q_e']
        rho += Q * np.exp(-0.5 * ((r_fm - r_n) / sigma_fm)**2) / (np.sqrt(2 * np.pi) * sigma_fm)
    return rho


def verify_charge_normalisation(harmonics, sigma_fm=0.2, label=''):
    """Check that charge density integrates to the expected total."""
    Q_total_expected = sum(h['Q_e'] for h in harmonics)
    r = np.linspace(0, 20, 5000)
    rho = nucleon_charge_density(r, harmonics, sigma_fm)
    dr = r[1] - r[0]
    Q_integrated = np.sum(rho) * dr
    print(f"  {label:8s}  Q_expected = {Q_total_expected:+.3f}e  "
          f"Q_integrated = {Q_integrated:+.4f}e  (Δ = {abs(Q_integrated - Q_total_expected):.1e})")
    return Q_integrated


# ---- 2. Ring-ring Coulomb integral ----

def ring_ring_coulomb(r1_fm, r2_fm, d_fm, n_theta=200):
    """Coulomb energy between two unit-charge coplanar rings at separation d.

    Ring 1: radius r1, centred at origin, in xy-plane.
    Ring 2: radius r2, centred at (d, 0, 0), in xy-plane.

    Uses the elliptic-integral reduction.  For two coplanar rings the
    double angular integral reduces to a single integral over theta2
    with the inner integral giving a complete elliptic integral K(k).

    Returns energy in MeV (for unit charges e).
    """
    if r1_fm < 1e-12 and r2_fm < 1e-12:
        # Both point charges
        if d_fm < 1e-12:
            return 0.0
        return k_e_MeV_fm / d_fm

    if r1_fm < 1e-12:
        # Ring 1 is a point at origin, ring 2 is a ring centred at d
        # E = k_e e² / (2π) ∫ dθ / |r_point - r_ring(θ)|
        # |r|² = r2² + d² - 2 r2 d cos(θ)
        def integrand(theta):
            dist2 = r2_fm**2 + d_fm**2 - 2 * r2_fm * d_fm * np.cos(theta)
            return 1.0 / np.sqrt(max(dist2, 1e-30))
        val, _ = quad(integrand, 0, 2 * np.pi, limit=200)
        return k_e_MeV_fm * val / (2 * np.pi)

    if r2_fm < 1e-12:
        return ring_ring_coulomb(r2_fm, r1_fm, d_fm, n_theta)

    # General case: numerical integration over theta2
    # For each theta2, the distance between point on ring1 at (r1, 0)
    # and point on ring2 at (d + r2 cos θ2, r2 sin θ2) averaged over
    # ring1 angle theta1 gives an elliptic integral.
    #
    # |P1(θ1) - P2(θ2)|² = (r1 cos θ1 - d - r2 cos θ2)² + (r1 sin θ1 - r2 sin θ2)²
    #   = r1² + (d + r2 cos θ2)² + (r2 sin θ2)² - 2 r1 [(d + r2 cos θ2) cos θ1 + r2 sin θ2 sin θ1]
    #
    # Let A² = r1² + (d + r2 cos θ2)² + r2² sin²θ2
    #     B  = r1 √[(d + r2 cos θ2)² + r2² sin²θ2]
    # Then inner integral ∫ dθ1 / √(A² - 2B cos(θ1 - φ))
    #   = (2/√(A² + 2B)) × 2 K(k)   where k² = 4B/(A² + 2B) ... wait,
    # let me use the standard form more carefully.
    #
    # Actually let's just do direct numerical double integration for robustness.

    theta1 = np.linspace(0, 2 * np.pi, n_theta, endpoint=False)
    dtheta1 = 2 * np.pi / n_theta
    theta2 = np.linspace(0, 2 * np.pi, n_theta, endpoint=False)
    dtheta2 = 2 * np.pi / n_theta

    # Ring 1 points: (r1 cos θ1, r1 sin θ1, 0)
    x1 = r1_fm * np.cos(theta1)
    y1 = r1_fm * np.sin(theta1)

    # Ring 2 points: (d + r2 cos θ2, r2 sin θ2, 0)
    x2 = d_fm + r2_fm * np.cos(theta2)
    y2 = r2_fm * np.sin(theta2)

    # Double sum: E = (k_e e²) / (2π)² × ΣΣ dθ1 dθ2 / |r1 - r2|
    total = 0.0
    for j in range(n_theta):
        dx = x1 - x2[j]
        dy = y1 - y2[j]
        dist = np.sqrt(dx**2 + dy**2)
        dist = np.maximum(dist, 1e-12)  # avoid division by zero
        total += np.sum(1.0 / dist)

    total *= dtheta1 * dtheta2
    energy = k_e_MeV_fm * total / (2 * np.pi)**2
    return energy


def verify_point_charge_limit():
    """Ring-ring Coulomb with r1=r2≈0, d>0 should → k_e e²/d."""
    d = 1.0  # fm
    E_computed = ring_ring_coulomb(0.001, 0.001, d, n_theta=100)
    E_expected = k_e_MeV_fm / d
    print(f"  Point-charge limit: E_computed = {E_computed:.6f} MeV, "
          f"E_expected = {E_expected:.6f} MeV, "
          f"ratio = {E_computed/E_expected:.6f}")


# ---- 3. Nucleon overlap energy ----

def nucleon_overlap_energy(harmonics_1, harmonics_2, d_fm, sigma_fm=0.2):
    """Total electrostatic energy between two nucleons at separation d.

    Sums over all pairs of harmonic rings, weighting by charge fractions.
    The finite ring width σ is accounted for by convolving the ring-ring
    Coulomb with Gaussian profiles — in practice we use a small number of
    sample radii around each ring centre.

    Returns E in MeV (negative = attractive).
    """
    # For each harmonic pair, compute the ring-ring Coulomb weighted by charges
    # To account for σ: sample each ring at a few radii and weight by Gaussian
    n_sigma_pts = 5  # points per ring for σ convolution
    sigma_range = 2.0  # sample out to ±2σ

    total_energy = 0.0
    for h1 in harmonics_1:
        r1_centre = harmonic_ring_radius(h1['E_MeV'])
        Q1 = h1['Q_e']
        for h2 in harmonics_2:
            r2_centre = harmonic_ring_radius(h2['E_MeV'])
            Q2 = h2['Q_e']

            if sigma_fm < 0.01:
                # Delta-function limit: just use ring centres
                E_pair = Q1 * Q2 * ring_ring_coulomb(r1_centre, r2_centre, d_fm)
                total_energy += E_pair
            else:
                # Gaussian convolution: sample radii
                r1_samples = np.linspace(
                    max(0.01, r1_centre - sigma_range * sigma_fm),
                    r1_centre + sigma_range * sigma_fm,
                    n_sigma_pts)
                r2_samples = np.linspace(
                    max(0.01, r2_centre - sigma_range * sigma_fm),
                    r2_centre + sigma_range * sigma_fm,
                    n_sigma_pts)
                dr1 = r1_samples[1] - r1_samples[0]
                dr2 = r2_samples[1] - r2_samples[0]

                # Gaussian weights
                w1 = np.exp(-0.5 * ((r1_samples - r1_centre) / sigma_fm)**2)
                w1 /= np.sum(w1)
                w2 = np.exp(-0.5 * ((r2_samples - r2_centre) / sigma_fm)**2)
                w2 /= np.sum(w2)

                E_conv = 0.0
                for i, r1 in enumerate(r1_samples):
                    for j, r2 in enumerate(r2_samples):
                        E_conv += w1[i] * w2[j] * ring_ring_coulomb(
                            r1, r2, d_fm, n_theta=100)
                total_energy += Q1 * Q2 * E_conv

    return total_energy


# ---- 4. Kinetic energy correction ----

def kinetic_energy_correction(d_fm, mu_MeV=469.5):
    """Zero-point kinetic energy of relative motion.

    E_kin = ℏc² / (2 μ d²)  in natural units = (ℏc)² / (2 μ c² d²)
    where μ is the reduced mass in MeV/c².

    For p-n system: μ = m_p m_n / (m_p + m_n) ≈ m_p/2 ≈ 469.5 MeV/c².
    """
    return hbar_c_MeV_fm**2 / (2 * mu_MeV * d_fm**2)


# ---- 5. Deuteron binding ----

def deuteron_binding(sigma_fm=0.2, d_range=(0.2, 8.0), n_points=60, verbose=True):
    """Scan p-n separation, find E_total = E_overlap + E_kin minimum.

    Compare to measured 2.224 MeV.
    """
    d_vals = np.linspace(d_range[0], d_range[1], n_points)
    E_overlap = np.zeros(n_points)
    E_kin = np.zeros(n_points)
    E_total = np.zeros(n_points)

    for i, d in enumerate(d_vals):
        E_overlap[i] = nucleon_overlap_energy(
            PROTON_HARMONICS, NEUTRON_HARMONICS, d, sigma_fm)
        E_kin[i] = kinetic_energy_correction(d)
        E_total[i] = E_overlap[i] + E_kin[i]

    idx_min = np.argmin(E_total)
    d_min = d_vals[idx_min]
    E_min = E_total[idx_min]

    if verbose:
        # Show pair-by-pair breakdown at a few separations
        print(f"\n  Deuteron pair-by-pair breakdown (σ = {sigma_fm:.2f} fm):")
        print(f"    {'d (fm)':>8s}  {'E_overlap':>10s}  {'E_kin':>10s}  {'E_total':>10s}")
        print(f"    " + "─" * 44)
        for d_show in [0.5, 1.0, 1.5, 2.0, 3.0, 5.0]:
            idx = np.argmin(np.abs(d_vals - d_show))
            print(f"    {d_vals[idx]:8.2f}  {E_overlap[idx]:+10.4f}  "
                  f"{E_kin[idx]:+10.4f}  {E_total[idx]:+10.4f}")

        # Show individual ring-pair contributions at d = 2 fm
        d_diag = 2.0
        print(f"\n  Ring-pair contributions at d = {d_diag:.1f} fm:")
        for h1 in PROTON_HARMONICS:
            r1 = harmonic_ring_radius(h1['E_MeV'])
            for h2 in NEUTRON_HARMONICS:
                r2 = harmonic_ring_radius(h2['E_MeV'])
                E_bare = ring_ring_coulomb(r1, r2, d_diag)
                E_weighted = h1['Q_e'] * h2['Q_e'] * E_bare
                sign = '+' if E_weighted > 0 else ''
                print(f"    {h1['name']:>4s}({h1['Q_e']:+.2f}e, r={r1:.2f}) × "
                      f"{h2['name']:>4s}({h2['Q_e']:+.2f}e, r={r2:.2f}): "
                      f"E_bare={E_bare:.4f}  E_net={sign}{E_weighted:.4f} MeV")

        print(f"\n  Deuteron (σ = {sigma_fm:.2f} fm):")
        print(f"    Equilibrium separation: {d_min:.3f} fm")
        print(f"    E_overlap at min:       {E_overlap[idx_min]:+.4f} MeV")
        print(f"    E_kinetic at min:       {E_kin[idx_min]:+.4f} MeV")
        print(f"    E_total at min:         {E_min:+.4f} MeV")
        print(f"    Measured BE(²H):        -2.224 MeV")
        if E_min < 0:
            print(f"    Ratio pred/expt:        {E_min / -2.224:.3f}")
        else:
            print(f"    *** No bound state found (E_min > 0) ***")

    return {
        'd_vals': d_vals,
        'E_overlap': E_overlap,
        'E_kin': E_kin,
        'E_total': E_total,
        'd_min': d_min,
        'E_min': E_min,
    }


# ---- 6. Alpha particle ----

def alpha_binding(sigma_fm=0.2, d_eq_fm=None, verbose=True):
    """⁴He: 2p + 2n in tetrahedral arrangement.

    6 pairwise distances, all equal for a regular tetrahedron.
    4 p-n pairs + 1 p-p pair + 1 n-n pair.
    """
    # If no equilibrium distance given, find from deuteron
    if d_eq_fm is None:
        deut = deuteron_binding(sigma_fm, verbose=False)
        d_eq_fm = deut['d_min']

    # Pairwise energies (all at same separation for tetrahedron)
    E_pn = nucleon_overlap_energy(PROTON_HARMONICS, NEUTRON_HARMONICS, d_eq_fm, sigma_fm)
    E_pp = nucleon_overlap_energy(PROTON_HARMONICS, PROTON_HARMONICS, d_eq_fm, sigma_fm)
    E_nn = nucleon_overlap_energy(NEUTRON_HARMONICS, NEUTRON_HARMONICS, d_eq_fm, sigma_fm)

    # 4 p-n + 1 p-p + 1 n-n
    E_total = 4 * E_pn + E_pp + E_nn

    # Kinetic energy: 4 nucleons in a volume of size ~ d_eq
    # Approximate as 4 independent particles each with ZPE ~ ℏc/(d)
    E_kin_total = 4 * kinetic_energy_correction(d_eq_fm, mu_MeV=938.3 * 3/4)

    E_bound = E_total + E_kin_total

    if verbose:
        print(f"\n  ⁴He alpha particle (σ = {sigma_fm:.2f} fm, d = {d_eq_fm:.2f} fm):")
        print(f"    E(p-n) × 4:   {4*E_pn:+.4f} MeV")
        print(f"    E(p-p) × 1:   {E_pp:+.4f} MeV")
        print(f"    E(n-n) × 1:   {E_nn:+.4f} MeV")
        print(f"    E_overlap:     {E_total:+.4f} MeV")
        print(f"    E_kinetic:     {E_kin_total:+.4f} MeV")
        print(f"    E_total:       {E_bound:+.4f} MeV")
        print(f"    Measured:      -28.296 MeV")
        if E_bound < 0:
            print(f"    Ratio:         {E_bound / -28.296:.3f}")

    return {
        'E_pn': E_pn, 'E_pp': E_pp, 'E_nn': E_nn,
        'E_overlap': E_total, 'E_kin': E_kin_total,
        'E_total': E_bound, 'd_fm': d_eq_fm,
    }


# ---- 7. Light nuclei survey ----

# Bond counts and geometry for light nuclei
# (Z, N, n_pn_bonds, n_pp_bonds, n_nn_bonds, geometry_note)
LIGHT_NUCLEI = {
    '2H':  (1, 1,  1, 0, 0, 'dimer'),
    '3H':  (1, 2,  2, 0, 1, 'triangle'),
    '3He': (2, 1,  2, 1, 0, 'triangle'),
    '4He': (2, 2,  4, 1, 1, 'tetrahedron'),
    '6Li': (3, 3,  9, 3, 3, 'octahedron'),
    '7Li': (3, 4, 12, 3, 6, 'pentagonal bipyramid'),
    '12C': (6, 6, 36, 15, 15, '3 alphas'),
    '16O': (8, 8, 64, 28, 28, '4 alphas'),
    '56Fe': (26, 30, 780, 325, 435, 'bulk'),
}


def light_nuclei_survey(sigma_fm=0.2, d_eq_fm=None, verbose=True):
    """Compute binding energies for light nuclei using bond counting.

    Uses Robinson's layer geometry: nearest-neighbour bonds only.
    """
    if d_eq_fm is None:
        deut = deuteron_binding(sigma_fm, verbose=False)
        d_eq_fm = deut['d_min']

    E_pn = nucleon_overlap_energy(PROTON_HARMONICS, NEUTRON_HARMONICS, d_eq_fm, sigma_fm)
    E_pp = nucleon_overlap_energy(PROTON_HARMONICS, PROTON_HARMONICS, d_eq_fm, sigma_fm)
    E_nn = nucleon_overlap_energy(NEUTRON_HARMONICS, NEUTRON_HARMONICS, d_eq_fm, sigma_fm)

    if verbose:
        print(f"\n  Light nuclei survey (σ = {sigma_fm:.2f} fm, d = {d_eq_fm:.2f} fm):")
        print(f"    Pair energies: E(p-n) = {E_pn:+.4f}, E(p-p) = {E_pp:+.4f}, E(n-n) = {E_nn:+.4f} MeV")
        print(f"\n    {'Nucleus':>8s}  {'Z':>3s} {'N':>3s}  {'Bonds':>8s}  "
              f"{'E_pred':>10s}  {'E_expt':>10s}  {'Ratio':>8s}  {'B/A pred':>10s}  {'B/A expt':>10s}")
        print(f"    " + "─" * 80)

    results = {}
    for name, (Z, N, n_pn, n_pp, n_nn, geom) in LIGHT_NUCLEI.items():
        A = Z + N
        E_pred = n_pn * E_pn + n_pp * E_pp + n_nn * E_nn
        E_expt = -EXPERIMENTAL_BE.get(name, 0)
        BpA_pred = -E_pred / A if E_pred < 0 else 0
        BpA_expt = EXPERIMENTAL_BE.get(name, 0) / A

        if verbose and name in EXPERIMENTAL_BE:
            ratio = E_pred / E_expt if E_expt != 0 else 0
            print(f"    {name:>8s}  {Z:3d} {N:3d}  {n_pn+n_pp+n_nn:8d}  "
                  f"{E_pred:+10.3f}  {E_expt:+10.3f}  {ratio:8.3f}  "
                  f"{BpA_pred:10.3f}  {BpA_expt:10.3f}")

        results[name] = {
            'Z': Z, 'N': N, 'A': A,
            'E_pred': E_pred, 'E_expt': E_expt,
            'BpA_pred': BpA_pred, 'BpA_expt': BpA_expt,
        }

    return results


# ---- 8. Parameter sensitivity ----

def sigma_scan(sigma_range=(0.05, 0.8), n_points=16, verbose=True):
    """Scan ring width σ, show deuteron BE vs σ."""
    sigmas = np.linspace(sigma_range[0], sigma_range[1], n_points)
    results = []

    if verbose:
        print(f"\n  σ scan for deuteron:")
        print(f"    {'σ (fm)':>8s}  {'d_eq (fm)':>10s}  {'E_min (MeV)':>12s}  {'vs -2.224':>10s}")
        print(f"    " + "─" * 44)

    for s in sigmas:
        deut = deuteron_binding(s, verbose=False)
        ratio = deut['E_min'] / -2.224 if deut['E_min'] < 0 else 0
        results.append((s, deut['d_min'], deut['E_min'], ratio))
        if verbose:
            print(f"    {s:8.3f}  {deut['d_min']:10.3f}  {deut['E_min']:+12.4f}  {ratio:10.3f}")

    return results


def charge_scan(verbose=True):
    """Scan proton/neutron charge assignments around Robinson's values."""
    if verbose:
        print(f"\n  Charge assignment scan (deuteron at σ = 0.20 fm):")
        print(f"    {'Label':>20s}  {'P/3':>6s} {'P/9':>6s}  {'N/3':>6s} {'N/9':>6s} {'N/27':>6s}  {'E_min':>10s}")
        print(f"    " + "─" * 70)

    configs = [
        ('Robinson (default)', [+2/3, +1/3], [-1/3, +2/3, -1/3]),
        ('Quark-like v2',      [+2/3, +1/3], [-1/3, +1/3,  0.0]),
        ('Equal thirds',       [+1/3, +2/3], [-1/3, +1/3,  0.0]),
        ('Strong inner',       [+1.0, +0.0], [-1/3, +2/3, -1/3]),
        ('Inverted proton',    [+1/3, +2/3], [-2/3, +1/3, +1/3]),
    ]

    results = []
    for label, p_charges, n_charges in configs:
        p_harm = [dict(h) for h in PROTON_HARMONICS]
        n_harm = [dict(h) for h in NEUTRON_HARMONICS]
        for i, q in enumerate(p_charges):
            p_harm[i]['Q_e'] = q
        for i, q in enumerate(n_charges):
            n_harm[i]['Q_e'] = q

        # Quick deuteron scan
        d_vals = np.linspace(0.5, 5.0, 30)
        E_best = np.inf
        d_best = 0
        for d in d_vals:
            E = nucleon_overlap_energy(p_harm, n_harm, d, 0.2) + kinetic_energy_correction(d)
            if E < E_best:
                E_best = E
                d_best = d

        results.append((label, p_charges, n_charges, d_best, E_best))
        if verbose:
            pq = ' '.join(f'{q:+.2f}' for q in p_charges)
            nq = ' '.join(f'{q:+.2f}' for q in n_charges)
            print(f"    {label:>20s}  {pq}  {nq}  {E_best:+10.4f}")

    return results


# ---- 9. Display & Plot ----

def make_figure(deut_result, alpha_result, nuclei_results, sigma_results):
    """Multi-panel figure: charge profiles, E(d), BE comparison, σ scan."""
    fig = plt.figure(figsize=(18, 14), facecolor='#0a0a0a')
    gs = GridSpec(2, 3, hspace=0.35, wspace=0.35)

    for ax_pos in [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]:
        ax = fig.add_subplot(gs[ax_pos[0], ax_pos[1]])
        ax.set_facecolor('#111111')
        ax.tick_params(colors='white')
        for spine in ax.spines.values():
            spine.set_color('#444444')
        ax.xaxis.label.set_color('white')
        ax.yaxis.label.set_color('white')
        ax.title.set_color('white')

    # Panel 1: Charge density profiles
    ax1 = fig.axes[0]
    r = np.linspace(0, 8, 500)
    rho_p = nucleon_charge_density(r, PROTON_HARMONICS, sigma_fm=0.2)
    rho_n = nucleon_charge_density(r, NEUTRON_HARMONICS, sigma_fm=0.2)
    ax1.fill_between(r, rho_p, 0, where=rho_p > 0, alpha=0.3, color='#ff4444')
    ax1.fill_between(r, rho_p, 0, where=rho_p < 0, alpha=0.3, color='#4444ff')
    ax1.plot(r, rho_p, '-', color='#ff4444', linewidth=2, label='Proton')
    ax1.fill_between(r, rho_n, 0, where=rho_n > 0, alpha=0.2, color='#00cc88')
    ax1.fill_between(r, rho_n, 0, where=rho_n < 0, alpha=0.2, color='#8800cc')
    ax1.plot(r, rho_n, '-', color='#00cc88', linewidth=2, label='Neutron')
    # Mark ring positions
    for h in PROTON_HARMONICS:
        rn = harmonic_ring_radius(h['E_MeV'])
        ax1.axvline(rn, color='#ff4444', alpha=0.3, linestyle=':')
        ax1.text(rn, ax1.get_ylim()[1] * 0.9 if hasattr(ax1, '_ylim_set') else 1.5,
                 h['name'], color='#ff4444', fontsize=7, ha='center')
    for h in NEUTRON_HARMONICS:
        rn = harmonic_ring_radius(h['E_MeV'])
        ax1.axvline(rn, color='#00cc88', alpha=0.3, linestyle=':')
    ax1.axhline(0, color='white', alpha=0.3, linewidth=0.5)
    ax1.set_xlabel('r (fm)')
    ax1.set_ylabel('ρ(r) (e/fm)')
    ax1.set_title('Robinson Charge Profiles (σ = 0.2 fm)', fontsize=11, color='white')
    ax1.legend(fontsize=8, facecolor='#222222', edgecolor='#444444', labelcolor='white')
    ax1.grid(True, alpha=0.2, color='white')

    # Panel 2: Deuteron E(d) curve
    ax2 = fig.axes[1]
    d = deut_result['d_vals']
    ax2.plot(d, deut_result['E_overlap'], '--', color='#00ccff', label='E_overlap', linewidth=1.5)
    ax2.plot(d, deut_result['E_kin'], '--', color='#ffaa00', label='E_kinetic', linewidth=1.5)
    ax2.plot(d, deut_result['E_total'], '-', color='#ff44ff', linewidth=2.5, label='E_total')
    ax2.axhline(-2.224, color='#00ff44', linestyle=':', alpha=0.7, label='Expt −2.224 MeV')
    ax2.axhline(0, color='white', alpha=0.3, linewidth=0.5)
    d_min = deut_result['d_min']
    E_min = deut_result['E_min']
    ax2.plot(d_min, E_min, 'o', color='#ff44ff', markersize=8, zorder=5)
    ax2.annotate(f'd={d_min:.2f} fm\nE={E_min:.3f} MeV',
                 xy=(d_min, E_min), xytext=(d_min + 0.5, E_min + 0.5),
                 color='white', fontsize=8,
                 arrowprops=dict(arrowstyle='->', color='white', lw=0.8))
    ax2.set_xlabel('p-n separation d (fm)')
    ax2.set_ylabel('Energy (MeV)')
    ax2.set_title('Deuteron Potential Curve', fontsize=11, color='white')
    ax2.legend(fontsize=8, facecolor='#222222', edgecolor='#444444', labelcolor='white')
    ax2.grid(True, alpha=0.2, color='white')
    ax2.set_xlim(0.3, 5.0)

    # Panel 3: Binding energy comparison (bar chart)
    ax3 = fig.axes[2]
    nuclei_names = [n for n in LIGHT_NUCLEI if n in EXPERIMENTAL_BE and n in nuclei_results]
    x = np.arange(len(nuclei_names))
    BE_pred = [-nuclei_results[n]['E_pred'] for n in nuclei_names]
    BE_expt = [EXPERIMENTAL_BE[n] for n in nuclei_names]
    width = 0.35
    ax3.bar(x - width/2, BE_expt, width, color='#00ff44', alpha=0.7, label='Experiment')
    ax3.bar(x + width/2, BE_pred, width, color='#ff44ff', alpha=0.7, label='Robinson model')
    ax3.set_xticks(x)
    ax3.set_xticklabels(nuclei_names, color='white', fontsize=8)
    ax3.set_ylabel('Total BE (MeV)')
    ax3.set_title('Binding Energy: Model vs Experiment', fontsize=11, color='white')
    ax3.legend(fontsize=8, facecolor='#222222', edgecolor='#444444', labelcolor='white')
    ax3.grid(True, alpha=0.2, color='white', axis='y')

    # Panel 4: B/A curve
    ax4 = fig.axes[3]
    A_pred = [nuclei_results[n]['A'] for n in nuclei_names]
    BpA_pred = [nuclei_results[n]['BpA_pred'] for n in nuclei_names]
    BpA_expt = [nuclei_results[n]['BpA_expt'] for n in nuclei_names]
    ax4.plot(A_pred, BpA_expt, 's-', color='#00ff44', markersize=6, label='Experiment')
    ax4.plot(A_pred, BpA_pred, 'o-', color='#ff44ff', markersize=6, label='Robinson model')
    ax4.set_xlabel('Mass number A')
    ax4.set_ylabel('B/A (MeV)')
    ax4.set_title('Binding Energy per Nucleon', fontsize=11, color='white')
    ax4.legend(fontsize=8, facecolor='#222222', edgecolor='#444444', labelcolor='white')
    ax4.grid(True, alpha=0.2, color='white')

    # Panel 5: σ scan
    ax5 = fig.axes[4]
    sigmas = [r[0] for r in sigma_results]
    E_mins = [r[2] for r in sigma_results]
    ax5.plot(sigmas, E_mins, 'o-', color='#00ccff', markersize=5, linewidth=2)
    ax5.axhline(-2.224, color='#00ff44', linestyle=':', alpha=0.7, label='Expt −2.224 MeV')
    ax5.axhline(0, color='white', alpha=0.3, linewidth=0.5)
    ax5.set_xlabel('Ring width σ (fm)')
    ax5.set_ylabel('Deuteron E_min (MeV)')
    ax5.set_title('σ Sensitivity: Deuteron BE', fontsize=11, color='white')
    ax5.legend(fontsize=8, facecolor='#222222', edgecolor='#444444', labelcolor='white')
    ax5.grid(True, alpha=0.2, color='white')

    # Panel 6: Overlapping charge profiles at finite separation
    ax6 = fig.axes[5]
    d_show = 2.0  # fm separation to illustrate
    r_plot = np.linspace(-2, 8, 600)
    rho_p_shifted = nucleon_charge_density(r_plot, PROTON_HARMONICS, 0.2)
    rho_n_shifted = nucleon_charge_density(r_plot - d_show, NEUTRON_HARMONICS, 0.2)
    ax6.fill_between(r_plot, rho_p_shifted, 0, where=rho_p_shifted > 0,
                     alpha=0.25, color='#ff4444')
    ax6.fill_between(r_plot, rho_p_shifted, 0, where=rho_p_shifted < 0,
                     alpha=0.25, color='#4444ff')
    ax6.plot(r_plot, rho_p_shifted, '-', color='#ff4444', linewidth=1.5, label='Proton @ 0')
    ax6.fill_between(r_plot, rho_n_shifted, 0, where=rho_n_shifted > 0,
                     alpha=0.2, color='#00cc88')
    ax6.fill_between(r_plot, rho_n_shifted, 0, where=rho_n_shifted < 0,
                     alpha=0.2, color='#8800cc')
    ax6.plot(r_plot, rho_n_shifted, '-', color='#00cc88', linewidth=1.5,
             label=f'Neutron @ d={d_show} fm')
    # Highlight overlap region
    overlap = rho_p_shifted * rho_n_shifted
    ax6.fill_between(r_plot, overlap * 5, 0, where=overlap < 0,
                     alpha=0.3, color='#ffff00', label='Attractive overlap (×5)')
    ax6.axhline(0, color='white', alpha=0.3, linewidth=0.5)
    ax6.set_xlabel('r (fm)')
    ax6.set_ylabel('ρ(r) (e/fm)')
    ax6.set_title(f'Charge Overlap at d = {d_show} fm', fontsize=11, color='white')
    ax6.legend(fontsize=7, facecolor='#222222', edgecolor='#444444', labelcolor='white')
    ax6.grid(True, alpha=0.2, color='white')

    fig.suptitle('Robinson Charge-Overlap Nuclear Binding — NWT Analysis',
                 fontsize=14, fontweight='bold', color='white', y=0.99)

    outpath = os.path.join(out_dir, 'nwt_robinson_nuclear.png')
    fig.savefig(outpath, dpi=180, bbox_inches='tight', facecolor='#0a0a0a')
    print(f"\n  Saved: {outpath}")
    return outpath


# ---- 10. Main ----

def main():
    print("=" * 70)
    print("  Robinson Charge-Overlap Nuclear Binding Model")
    print("  Electrostatic overlap of harmonic charge rings")
    print("  Energy scale: k_e·e²/fm = {:.4f} MeV".format(k_e_MeV_fm))
    print("=" * 70)

    # ---- Harmonic ring data ----
    print("\n--- Robinson's Nucleon Harmonics ---")
    print(f"  {'Ring':>6s}  {'E (MeV)':>10s}  {'r (fm)':>10s}  {'Q (e)':>8s}")
    print(f"  " + "─" * 40)
    for h in PROTON_HARMONICS:
        r = harmonic_ring_radius(h['E_MeV'])
        print(f"  {h['name']:>6s}  {h['E_MeV']:10.1f}  {r:10.3f}  {h['Q_e']:+8.3f}")
    for h in NEUTRON_HARMONICS:
        r = harmonic_ring_radius(h['E_MeV'])
        print(f"  {h['name']:>6s}  {h['E_MeV']:10.1f}  {r:10.3f}  {h['Q_e']:+8.3f}")

    # ---- Verify charge normalisation ----
    print("\n--- Charge Normalisation ---")
    verify_charge_normalisation(PROTON_HARMONICS, 0.2, 'Proton')
    verify_charge_normalisation(NEUTRON_HARMONICS, 0.2, 'Neutron')

    # ---- Verify point-charge limit ----
    print("\n--- Point-Charge Limit Verification ---")
    verify_point_charge_limit()

    # ---- Deuteron ----
    print("\n" + "=" * 70)
    print("  DEUTERON BINDING")
    print("=" * 70)
    deut = deuteron_binding(sigma_fm=0.2)

    # ---- Alpha particle ----
    print("\n" + "=" * 70)
    print("  ALPHA PARTICLE")
    print("=" * 70)
    alpha_res = alpha_binding(sigma_fm=0.2)

    # ---- Light nuclei survey ----
    print("\n" + "=" * 70)
    print("  LIGHT NUCLEI SURVEY")
    print("=" * 70)
    nuclei = light_nuclei_survey(sigma_fm=0.2)

    # ---- σ sensitivity ----
    print("\n" + "=" * 70)
    print("  PARAMETER SENSITIVITY")
    print("=" * 70)
    sigma_res = sigma_scan()

    # ---- Charge assignment scan ----
    charge_scan()

    # ---- Figure ----
    print("\n--- Generating figure ---")
    make_figure(deut, alpha_res, nuclei, sigma_res)

    # ---- Verdict ----
    print("\n" + "=" * 70)
    print("  VERDICT")
    print("=" * 70)
    E_deut = deut['E_min']
    E_alpha = alpha_res['E_total']
    print(f"\n  Deuteron:  E_pred = {E_deut:+.4f} MeV  vs  E_expt = -2.224 MeV")
    if E_deut < 0:
        print(f"             Ratio = {E_deut / -2.224:.3f}")
    else:
        print(f"             *** NOT BOUND — model produces repulsion ***")
    print(f"  Alpha:     E_pred = {E_alpha:+.4f} MeV  vs  E_expt = -28.296 MeV")
    if E_alpha < 0:
        print(f"             Ratio = {E_alpha / -28.296:.3f}")

    print(f"\n  Robinson mechanism assessment:")
    print(f"    Energy scale k_e·e²/fm = {k_e_MeV_fm:.4f} MeV — correct ballpark")

    if E_deut < 0 and abs(E_deut / -2.224 - 1) < 1:
        print(f"    Deuteron BE within factor of 2 — mechanism is viable")
    elif E_deut < 0:
        print(f"    Deuteron BE off by factor {-2.224 / abs(E_deut):.1f} — needs tuning")
    else:
        print(f"    Deuteron not bound — charge assignments or σ need adjustment")

    # Find best σ
    best_sigma = None
    best_ratio = 0
    for s, d, E, ratio in sigma_res:
        if E < 0 and abs(ratio - 1) < abs(best_ratio - 1):
            best_ratio = ratio
            best_sigma = s
    if best_sigma is not None:
        print(f"\n  Best σ for deuteron: {best_sigma:.3f} fm (ratio = {best_ratio:.3f})")
    else:
        print(f"\n  No σ in scan range produces bound deuteron")

    print("\n" + "=" * 70)


if __name__ == '__main__':
    main()
