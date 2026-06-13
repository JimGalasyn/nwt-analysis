#!/usr/bin/env python3
"""
Spin(7) fusion ring computation via Freudenthal weight recursion
+ Brauer/Klimyk tensor-product algorithm.

From-scratch implementation in pure Python with rationals (Fraction).
No external Lie-algebra libraries required.

Pipeline:
  1. Enumerate all dominant weights μ ≤ λ for each rep.
  2. Compute weight multiplicities m_λ(μ) via Freudenthal recursion.
  3. Verify dim(V_λ) = Σ_μ m_λ(μ) · |Weyl orbit(μ)| matches Weyl
     dimension formula.
  4. Decompose V_a ⊗ V_b via Brauer / Klimyk: for each weight w of V_b,
     reflect (hwt(a) + w + ρ) to dominant chamber, sign-track, subtract ρ.
  5. Verify against Phase 6.3's V⊗V, V⊗Adj, Adj⊗Adj.
  6. Compute all unordered pairs in the 6-channel fusion ring.
  7. Enumerate all unordered triples (α,β,γ) with 1 ∈ α⊗β⊗γ.

Channels (Adj ⊗ Adj decomposition, B_3 / Spin(7)):
  {1, Adj, Sym²V, Λ³V, 168, (2,1,1)_189}
  with hwts {(0,0,0), (1,1,0), (2,0,0), (1,1,1), (2,2,0), (2,1,1)}.

CORRECTION (2026-04-25):  earlier project memory listed Adj⊗Adj as
the 8-component sum {1, V, Adj, Sym², Λ³V, 77, 105, 168}.  This is
WRONG.  Direct Freudenthal + Brauer (and hand-verification of the
V-multiplicity by exhaustive sum over all 21 Adj-weights) gives
6 components: {1, Adj, Sym²V, Λ³V, 168, (2,1,1)_189}.  The reps 77
and 105 are real B_3 irreps but they appear in V⊗Sym²V, not in
Adj⊗Adj.  The 189-dim rep (2,1,1) takes their place.

We also retain V, 77, 105 in EXTENDED_REPS so the Brauer algorithm
can name them when they appear in larger products (e.g. V⊗Sym2 =
V + 77 + 105).
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from itertools import combinations_with_replacement, product, permutations
from typing import Dict, List, Tuple, Set, Iterable

Vec = Tuple[Fraction, Fraction, Fraction]


# =====================================================================
# 1. B_3 root system data (orthogonal basis, rationals)
# =====================================================================

SIMPLE_ROOTS: List[Vec] = [
    (Fraction(1), Fraction(-1), Fraction(0)),     # α_1 = e_1 - e_2  (long)
    (Fraction(0), Fraction(1),  Fraction(-1)),    # α_2 = e_2 - e_3  (long)
    (Fraction(0), Fraction(0),  Fraction(1)),     # α_3 = e_3        (short)
]

POSITIVE_ROOTS: List[Vec] = [
    (Fraction(1),  Fraction(-1), Fraction(0)),
    (Fraction(1),  Fraction(0),  Fraction(-1)),
    (Fraction(0),  Fraction(1),  Fraction(-1)),
    (Fraction(1),  Fraction(1),  Fraction(0)),
    (Fraction(1),  Fraction(0),  Fraction(1)),
    (Fraction(0),  Fraction(1),  Fraction(1)),
    (Fraction(1),  Fraction(0),  Fraction(0)),
    (Fraction(0),  Fraction(1),  Fraction(0)),
    (Fraction(0),  Fraction(0),  Fraction(1)),
]

RHO: Vec = (Fraction(5, 2), Fraction(3, 2), Fraction(1, 2))


def vadd(a: Vec, b: Vec) -> Vec:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def vsub(a: Vec, b: Vec) -> Vec:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def vscale(c: Fraction, a: Vec) -> Vec:
    return (c * a[0], c * a[1], c * a[2])


def vdot(a: Vec, b: Vec) -> Fraction:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


# =====================================================================
# 2. Weyl group (signed permutations of B_3) and dominance
# =====================================================================

def reflect_to_dominant(v: Vec) -> Tuple[Vec, int]:
    """Apply Weyl reflections to bring v into the dominant chamber
    (v_1 ≥ v_2 ≥ v_3 ≥ 0).  Returns (v_dominant, sign) where sign
    is the parity of the number of reflections used.
    """
    a, b, c = v
    sign = 1
    if a < 0:
        a = -a; sign = -sign
    if b < 0:
        b = -b; sign = -sign
    if c < 0:
        c = -c; sign = -sign
    # Bubble-sort descending (each swap is one Weyl reflection)
    if a < b:
        a, b = b, a; sign = -sign
    if b < c:
        b, c = c, b; sign = -sign
    if a < b:
        a, b = b, a; sign = -sign
    return (a, b, c), sign


def is_dominant(v: Vec) -> bool:
    return v[0] >= v[1] >= v[2] >= 0


def is_below_hwt(hwt: Vec, mu: Vec) -> bool:
    """True iff hwt - mu is a non-negative integer combination of
    simple roots.  For B_3 with simple roots α_1=(1,-1,0), α_2=(0,1,-1),
    α_3=(0,0,1):
        diff = hwt - mu = c_1 α_1 + c_2 α_2 + c_3 α_3
        ⇒ c_1 = diff_0, c_2 = diff_0 + diff_1, c_3 = diff_0 + diff_1 + diff_2.
    """
    d0 = hwt[0] - mu[0]
    d1 = hwt[1] - mu[1]
    d2 = hwt[2] - mu[2]
    c1 = d0
    c2 = d0 + d1
    c3 = d0 + d1 + d2
    for c in (c1, c2, c3):
        if c < 0 or c.denominator != 1:
            return False
    return True


def weyl_orbit(v: Vec) -> Set[Vec]:
    """Generate the Weyl orbit of v (signed permutations)."""
    seen: Set[Vec] = set()
    for perm in permutations(range(3)):
        for s in product([1, -1], repeat=3):
            new = (Fraction(s[0]) * v[perm[0]],
                   Fraction(s[1]) * v[perm[1]],
                   Fraction(s[2]) * v[perm[2]])
            seen.add(new)
    return seen


def weyl_orbit_size(v: Vec) -> int:
    return len(weyl_orbit(v))


# =====================================================================
# 3. Enumerate all dominant weights ≤ hwt
# =====================================================================

def all_dominant_below(hwt: Vec, c_max: int = 30) -> Set[Vec]:
    """Enumerate dominant weights μ with hwt - μ a non-negative integer
    combination of simple roots, c_i ≤ c_max.
    """
    found: Set[Vec] = set()
    for c1 in range(c_max + 1):
        for c2 in range(c_max + 1):
            for c3 in range(c_max + 1):
                # μ = hwt - c1*α_1 - c2*α_2 - c3*α_3
                #    = (hwt_0 - c1, hwt_1 + c1 - c2, hwt_2 + c2 - c3)
                mu = (hwt[0] - Fraction(c1),
                      hwt[1] + Fraction(c1 - c2),
                      hwt[2] + Fraction(c2 - c3))
                if mu[0] >= mu[1] >= mu[2] >= Fraction(0):
                    found.add(mu)
    return found


# =====================================================================
# 4. Freudenthal recursion for weight multiplicities
# =====================================================================
#
#   m_λ(μ) · [|λ+ρ|² - |μ+ρ|²]
#       = 2 · Σ_{α>0} Σ_{k≥1} ⟨μ + kα, α⟩ · m_λ(μ + kα).
#
# We process dominant weights in order of decreasing height (μ, ρ̌)
# where ρ̌ = (3, 2, 1) is the dual Weyl vector for B_3.  Multiplicities
# are Weyl-invariant, so m_λ(w) = m_λ(reflect_to_dominant(w)).  For
# any non-dominant w in the support, height(reflect(w)) > height(w),
# so the dominant rep is already computed.
# =====================================================================

RHO_CHECK: Vec = (Fraction(3), Fraction(2), Fraction(1))


def height(mu: Vec) -> Fraction:
    """Height in the ρ̌ inner product."""
    return vdot(mu, RHO_CHECK)


def freudenthal(hwt: Vec) -> Dict[Vec, Fraction]:
    """Compute m_λ(μ) for all dominant μ ≤ λ via Freudenthal recursion.
    Returns a dict mapping dominant weight -> multiplicity.
    """
    dom = all_dominant_below(hwt)
    sorted_dom = sorted(dom, key=lambda mu: -height(mu))

    mult: Dict[Vec, Fraction] = {hwt: Fraction(1)}
    lpr = vadd(hwt, RHO)
    lpr_sq = vdot(lpr, lpr)

    for mu in sorted_dom:
        if mu == hwt:
            continue
        mpr = vadd(mu, RHO)
        denom = lpr_sq - vdot(mpr, mpr)
        if denom == 0:
            # μ in same Weyl orbit as λ — but μ is dominant ≠ λ, so
            # this shouldn't happen in practice.
            mult[mu] = Fraction(0)
            continue

        rhs = Fraction(0)
        for alpha in POSITIVE_ROOTS:
            for k in range(1, 50):
                w = vadd(mu, vscale(Fraction(k), alpha))
                dom_w, _ = reflect_to_dominant(w)
                if dom_w not in dom:
                    break
                m_w = mult.get(dom_w, Fraction(0))
                if m_w == 0:
                    # If 0 (genuinely) we keep going — could ascend
                    # past it.  But if dom_w not yet in dict, it's a
                    # bug.
                    if dom_w not in mult:
                        raise RuntimeError(
                            f"Freudenthal ordering bug: μ={mu}, α={alpha}, "
                            f"k={k}, dom_w={dom_w} not yet computed")
                    continue
                rhs += vdot(w, alpha) * m_w

        mult[mu] = (Fraction(2) * rhs) / denom

    return mult


# =====================================================================
# 5. Weyl dimension formula (verification)
# =====================================================================

def weyl_dim(hwt: Vec) -> Fraction:
    num = Fraction(1)
    den = Fraction(1)
    lpr = vadd(hwt, RHO)
    for alpha in POSITIVE_ROOTS:
        num *= vdot(lpr, alpha)
        den *= vdot(RHO, alpha)
    return num / den


def total_dim_from_mults(mult: Dict[Vec, Fraction]) -> Fraction:
    """Sum of (multiplicity × Weyl orbit size) over dominant weights."""
    total = Fraction(0)
    for mu, m in mult.items():
        total += m * weyl_orbit_size(mu)
    return total


# =====================================================================
# 6. Reps in the Adj⊗Adj fusion algebra (8 channels)
# =====================================================================

# Adj⊗Adj fusion channels (the actual 6 components for B_3 = Spin(7))
REPS: Dict[str, Tuple[Vec, int]] = {
    '1':     ((Fraction(0), Fraction(0), Fraction(0)),   1),
    'Adj':   ((Fraction(1), Fraction(1), Fraction(0)),  21),
    'Sym2':  ((Fraction(2), Fraction(0), Fraction(0)),  27),
    'L3V':   ((Fraction(1), Fraction(1), Fraction(1)),  35),
    '168':   ((Fraction(2), Fraction(2), Fraction(0)), 168),
    '189':   ((Fraction(2), Fraction(1), Fraction(1)), 189),   # (2,1,1) — replaces V+77+105
}

# Extended reps that may appear in higher products (V⊗Sym2 etc.)
EXTENDED_REPS: Dict[str, Tuple[Vec, int]] = {
    **REPS,
    'V':     ((Fraction(1), Fraction(0), Fraction(0)),   7),
    '77':    ((Fraction(3), Fraction(0), Fraction(0)),  77),
    '105':   ((Fraction(2), Fraction(1), Fraction(0)), 105),
    '8s':    ((Fraction(1, 2), Fraction(1, 2), Fraction(1, 2)), 8),  # spinor (not used here)
}

# Cache: name -> {dominant μ : m_λ(μ)}
_DOMINANT_MULTS: Dict[str, Dict[Vec, Fraction]] = {}


def get_dominant_mults(name: str) -> Dict[Vec, Fraction]:
    if name not in _DOMINANT_MULTS:
        hwt, _ = EXTENDED_REPS[name]
        _DOMINANT_MULTS[name] = freudenthal(hwt)
    return _DOMINANT_MULTS[name]


def expand_to_full_weights(mult_dom: Dict[Vec, Fraction]) -> Dict[Vec, Fraction]:
    """Expand dominant-weight multiplicities to full weight set via Weyl orbits."""
    out: Dict[Vec, Fraction] = defaultdict(lambda: Fraction(0))
    for mu_dom, m in mult_dom.items():
        for mu in weyl_orbit(mu_dom):
            out[mu] += m
    return dict(out)


# =====================================================================
# 7. Brauer / Klimyk tensor-product algorithm
# =====================================================================
#
# V_λ ⊗ V_μ = ⊕_ν N^ν_{λμ} V_ν, with
#
#   N^ν_{λμ} = Σ_{w ∈ wts(V_μ)} sign(σ) · m_μ(w)
#               where (λ + w + ρ) reflects via σ ∈ Weyl to (ν + ρ).
#
# Skip contributions where (λ + w + ρ) lies on a Weyl wall (reflection
# would have a fixed point), which corresponds to ν not being strictly
# dominant after subtracting ρ.
# =====================================================================

def tensor_decompose(name_a: str, name_b: str) -> Dict[Vec, Fraction]:
    hwt_a, _ = EXTENDED_REPS[name_a]
    mult_b_dom = get_dominant_mults(name_b)
    mult_b_full = expand_to_full_weights(mult_b_dom)

    decomp: Dict[Vec, Fraction] = defaultdict(lambda: Fraction(0))

    for w, m_w in mult_b_full.items():
        cand_plus_rho = vadd(vadd(hwt_a, w), RHO)
        dom, sign = reflect_to_dominant(cand_plus_rho)
        nu = vsub(dom, RHO)
        # ν must be dominant integral (already guaranteed integral
        # since hwt_a, w, ρ-shift all match).  Wall contributions
        # automatically vanish: dom on a Weyl wall ⇒ dom has equal
        # coords or zero coord ⇒ ν fails dominance after subtracting ρ.
        if nu[0] >= nu[1] >= nu[2] >= 0:
            decomp[nu] += sign * m_w

    return {k: v for k, v in decomp.items() if v != 0}


# =====================================================================
# 8. Map highest weights back to channel names
# =====================================================================

NAME_BY_HWT: Dict[Vec, str] = {hwt: name for name, (hwt, _) in EXTENDED_REPS.items()}


def hwt_label(hwt: Vec) -> str:
    if hwt in NAME_BY_HWT:
        return NAME_BY_HWT[hwt]
    return f"({int(hwt[0])},{int(hwt[1])},{int(hwt[2])})"


def hwt_dim(hwt: Vec) -> int:
    if hwt in NAME_BY_HWT:
        return EXTENDED_REPS[NAME_BY_HWT[hwt]][1]
    return int(weyl_dim(hwt))


# =====================================================================
# 9. Verification
# =====================================================================

def section(title: str) -> None:
    print()
    print('=' * 78)
    print(' ' + title)
    print('=' * 78)


def verify_dimensions() -> bool:
    section('Verification 1: dim(V_λ) from Freudenthal vs Weyl formula')
    ok = True
    for name, (hwt, dim) in EXTENDED_REPS.items():
        m = get_dominant_mults(name)
        d_freud = total_dim_from_mults(m)
        d_weyl = weyl_dim(hwt)
        match = (d_freud == dim and d_weyl == dim)
        ok = ok and match
        print(f"  {name:>5}  hwt={hwt_label(hwt):<6} "
              f"Weyl-dim={int(d_weyl):>4}  Freud-dim={int(d_freud):>4}  "
              f"expected={dim:>4}  {'OK' if match else 'FAIL'}")
    return ok


def verify_known_products() -> bool:
    section('Verification 2: known products (analytic + cross-checks)')

    # Analytical expectations (verified independently).
    # Note: V is in EXTENDED_REPS (not in fusion-channel REPS) — used for
    # cross-checks with the standard Phase 6 V-decompositions.
    expected = {
        # V × V = 1 + Adj + Sym²V                       (so(7) standard)
        ('V', 'V'):     ['1', 'Adj', 'Sym2'],
        # V × Adj = V + Λ³V + (2,1,0)_105               (Phase 6, verified)
        ('V', 'Adj'):   ['V', 'L3V', '105'],
        # V × Sym² = V + 77 + 105                       (Sym³V = V + 77 in so(7))
        ('V', 'Sym2'):  ['V', '77', '105'],
        # Adj × Adj — CORRECTED (was wrong in Phase 6.3 memory).
        # The 8-component {1,V,Adj,Sym2,L3V,77,105,168} list cited in
        # earlier memory is actually wrong: the 6-component decomposition
        # below has the same total dim 441 = 1+21+27+35+168+189.
        ('Adj', 'Adj'): ['1', 'Adj', 'Sym2', 'L3V', '168', '189'],
    }

    ok = True
    for (a, b), expected_components in expected.items():
        decomp = tensor_decompose(a, b)
        d_in = EXTENDED_REPS[a][1] * EXTENDED_REPS[b][1]
        d_out = sum(hwt_dim(k) * v for k, v in decomp.items())
        match_dim = (d_out == d_in)
        names = [hwt_label(k) for k, v in decomp.items() if v > 0]
        match_names = sorted(names) == sorted(expected_components)
        match = match_dim and match_names
        ok = ok and match

        comp_str = ' + '.join(
            f"{int(v)}·{hwt_label(k)}" if v != 1 else hwt_label(k)
            for k, v in sorted(decomp.items(), key=lambda x: -hwt_dim(x[0]))
        )
        print(f"  {a} ⊗ {b}:")
        print(f"    computed: {comp_str}")
        print(f"    expected: {' + '.join(expected_components)}")
        print(f"    dim {EXTENDED_REPS[a][1]} × {EXTENDED_REPS[b][1]} = {d_in}; "
              f"sum out = {int(d_out)}  {'OK' if match else 'FAIL'}")
    return ok


# =====================================================================
# 10. Full fusion ring on 8 channels + singlet-triple enumeration
# =====================================================================

def compute_full_fusion() -> Dict[Tuple[str, str], Dict[str, int]]:
    """Compute α ⊗ β for all unordered pairs of fusion channels.
    Returns dict mapping (name pair) -> {component name: multiplicity}.

    Components outside the 6 channels are listed by hwt-tuple label.
    """
    section('Full fusion-ring computation: unordered pairs of channels')

    channel_names = list(REPS.keys())
    fusion: Dict[Tuple[str, str], Dict[str, int]] = {}

    for a, b in combinations_with_replacement(channel_names, 2):
        decomp = tensor_decompose(a, b)
        d_in = REPS[a][1] * REPS[b][1]
        d_out = sum(hwt_dim(k) * v for k, v in decomp.items())
        decomp_named: Dict[str, int] = defaultdict(int)
        in_channel = 0
        out_channel = 0
        for hwt, m in decomp.items():
            label = hwt_label(hwt)
            decomp_named[label] += int(m)
            if label in REPS:
                in_channel += hwt_dim(hwt) * int(m)
            else:
                out_channel += hwt_dim(hwt) * int(m)
        fusion[(a, b)] = dict(decomp_named)

        # Print
        comp_str = ' + '.join(
            f"{m}·{lbl}" if m != 1 else lbl
            for lbl, m in sorted(
                decomp_named.items(),
                key=lambda x: -hwt_dim(REPS[x[0]][0]) if x[0] in REPS else -1
            )
        )
        check = 'OK' if d_out == d_in else 'FAIL'
        print(f"  {a:>4} ⊗ {b:<4} = {comp_str}")
        print(f"           dim {REPS[a][1]:>3} × {REPS[b][1]:>3} = {d_in:>5}  "
              f"sum {int(d_out):>5}  in-ch {in_channel:>5}  out-ch {out_channel:>5}  "
              f"{check}")
    return fusion


def enumerate_singlet_triples(
    fusion: Dict[Tuple[str, str], Dict[str, int]]
) -> List[Tuple[str, str, str]]:
    """Find all unordered triples (a, b, c) of channels with 1 ∈ a⊗b⊗c.

    Self-conjugacy:  for these reps, V_λ* = V_λ, so 1 ∈ a⊗b⊗c iff c
    appears in a⊗b.
    """
    section('Singlet-containing triples (1 ∈ α⊗β⊗γ)')

    channels = list(REPS.keys())
    triples = []

    for a, b, c in combinations_with_replacement(channels, 3):
        # Check c in a⊗b
        # We need the unordered pair (a,b) lookup
        if (a, b) in fusion:
            decomp = fusion[(a, b)]
        elif (b, a) in fusion:
            decomp = fusion[(b, a)]
        else:
            continue
        if decomp.get(c, 0) > 0:
            triples.append((a, b, c))

    print(f"\n  Found {len(triples)} singlet-containing triples "
          f"(out of {len(list(combinations_with_replacement(channels, 3)))} total):")
    for t in triples:
        print(f"    {t}")

    return triples


# =====================================================================
# 11. Main
# =====================================================================

def main() -> None:
    section('Spin(7) fusion ring via Freudenthal + Brauer')

    print(f"""
Computing weight multiplicities, dimensions, and tensor-product
decompositions for the 6 channels {{1, Adj, Sym²V, Λ³V, 168, (2,1,1)_189}}
of Adj⊗Adj in B_3 = Spin(7).

(Note: V, 77, 105 are real B_3 reps but do NOT appear in Adj⊗Adj.
The earlier Phase 6.3 memory listed an 8-component decomposition that
was incorrect; the actual decomposition has 6 components.)

All arithmetic in exact rationals (fractions.Fraction).
""")

    ok_dim = verify_dimensions()
    if not ok_dim:
        print("\n  *** Dimension check FAILED — aborting.")
        return

    ok_prod = verify_known_products()
    if not ok_prod:
        print("\n  *** Known-product check FAILED — aborting.")
        return

    fusion = compute_full_fusion()
    triples = enumerate_singlet_triples(fusion)

    section('Summary')
    n_pairs = len(list(combinations_with_replacement(REPS.keys(), 2)))
    n_triples_total = len(list(combinations_with_replacement(REPS.keys(), 3)))
    print(f"""
  Reps verified (Freudenthal + Weyl agree):       {len(EXTENDED_REPS)} / {len(EXTENDED_REPS)}
  Adj⊗Adj decomposition (corrected):              6 components
                                                  {{1, Adj, Sym²V, Λ³V, 168, 189}}
  Known products verified:                        OK
  Total unordered pairs in fusion ring:           {n_pairs}
  Total unordered triples enumerated:             {n_triples_total}
  Singlet-containing triples found:               {len(triples)}

These triples constrain the V_6(Adj⁶) vertex factors at each of the
7 vertices of K_7 in the Spin(7)_k=32 RT amplitude.

NEXT STEP:  for each of the {len(triples)} triples, look up Racah-Wigner 6j
symbols of B_3 (Spin(7)) at level k=32 and assemble the K_7
amplitude in the Schur-Racah basis.
""")


if __name__ == '__main__':
    main()
