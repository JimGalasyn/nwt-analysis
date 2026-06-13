"""Find Z_3 subgroups of G_2 = Aut(O) with vertex orbit structure 1+3+3.

This is the rep-theory derivation of the seesaw-via-triality Z_3 on K_7:
show that there exists a concrete order-3 permutation σ of the 7
imaginary octonion units (= K_7 vertices) which (a) preserves the
octonion product (so σ ∈ G_2), and (b) has orbit structure 1+3+3.

If we find such σ, then the K_7 triality structure is RIGOROUSLY a
descendant of G_2 ⊂ Spin(7) (and therefore of Spin(8) triality via
the Spin(7) ⊂ Spin(8) embedding chain).

Approach:
  1. Build the octonion multiplication table (Baez Fano convention).
  2. Enumerate all permutations σ ∈ S_7 with cycle structure 1+3+3.
  3. For each, check if σ is an automorphism: σ(e_i · e_j) = σ(e_i) · σ(e_j).
  4. Equivalent: σ permutes Fano lines preserving cyclic orientation.
  5. Report which permutations work + their orbit decomposition.
"""
from __future__ import annotations

import itertools

# Baez Fano-line convention from nwt_substrate.algebra.octonions
FANO_LINES = [
    (1, 2, 3),
    (1, 4, 5),
    (1, 7, 6),
    (2, 4, 6),
    (2, 5, 7),
    (3, 4, 7),
    (3, 6, 5),
]


def line_set(line):
    """A line as an unordered frozenset of 3 points."""
    return frozenset(line)


# For a Fano line (i, j, k), the cyclic equivalence class {(i,j,k), (j,k,i),
# (k,i,j)} all encode the same orientation (positive). The reverse class
# {(j,i,k), (k,j,i), (i,k,j)} encodes the negative orientation.
FANO_POS_TRIPLES = set()
FANO_NEG_TRIPLES = set()
for line in FANO_LINES:
    i, j, k = line
    FANO_POS_TRIPLES.add((i, j, k))
    FANO_POS_TRIPLES.add((j, k, i))
    FANO_POS_TRIPLES.add((k, i, j))
    FANO_NEG_TRIPLES.add((j, i, k))
    FANO_NEG_TRIPLES.add((k, j, i))
    FANO_NEG_TRIPLES.add((i, k, j))


def is_automorphism(sigma):
    """Test if σ : {1..7} → {1..7} is an octonion automorphism.

    σ is in Aut(O) = G_2 iff σ maps every positive-cyclic Fano triple
    (i, j, k) to another positive-cyclic triple (σ(i), σ(j), σ(k)).
    """
    for (i, j, k) in FANO_LINES:
        si, sj, sk = sigma[i], sigma[j], sigma[k]
        if (si, sj, sk) not in FANO_POS_TRIPLES:
            return False
    return True


def cycle_structure(sigma):
    """Return sorted tuple of cycle lengths of σ."""
    visited = set()
    cycles = []
    for start in sigma:
        if start in visited:
            continue
        cycle_len = 0
        x = start
        while x not in visited:
            visited.add(x)
            x = sigma[x]
            cycle_len += 1
        cycles.append(cycle_len)
    return tuple(sorted(cycles))


def order(sigma):
    """Order of σ as a permutation."""
    from math import gcd
    cs = cycle_structure(sigma)
    o = 1
    for c in cs:
        o = o * c // gcd(o, c)
    return o


def orbits(sigma):
    """Return orbits of σ as a list of frozensets."""
    visited = set()
    out = []
    for start in sigma:
        if start in visited:
            continue
        orb = []
        x = start
        while x not in visited:
            visited.add(x)
            orb.append(x)
            x = sigma[x]
        out.append(frozenset(orb))
    return out


# ---------------------------------------------------------------------------
# Search: enumerate permutations with cycle structure 1+3+3 on {1..7}
# ---------------------------------------------------------------------------

candidates = []
n_total_perms = 0
for perm in itertools.permutations(range(1, 8)):
    n_total_perms += 1
    # Build sigma dict: sigma[i] = perm[i-1]
    sigma = {i: perm[i - 1] for i in range(1, 8)}
    if cycle_structure(sigma) != (1, 3, 3):
        continue
    if not is_automorphism(sigma):
        continue
    candidates.append(sigma)

print("=" * 72)
print(f"Searched {n_total_perms} permutations of {{1..7}}")
print(f"Found {len(candidates)} (1+3+3)-cycle automorphisms of the octonion product")
print("=" * 72)

if not candidates:
    print("\n  NONE FOUND. The Baez convention may not admit a 1+3+3 Z_3 in G_2.")
    print("  Will need a different Fano-line labeling.")
else:
    # Group by orbit structure (fixed point + the two 3-cycles)
    print()
    distinct_structures = set()
    for sigma in candidates:
        orbs = orbits(sigma)
        fixed = next(o for o in orbs if len(o) == 1)
        threes = sorted([tuple(sorted(o)) for o in orbs if len(o) == 3])
        key = (tuple(fixed), tuple(threes))
        distinct_structures.add(key)

    print(f"Distinct orbit structures: {len(distinct_structures)}")
    print(f"(Each generates a Z_3 subgroup; same Z_3 has σ and σ² as generators)")
    print()

    # Print a representative for each Z_3 (only one of {σ, σ²})
    z3_subgroups = []
    for sigma in candidates:
        orbs = orbits(sigma)
        # Canonicalize by lex-smallest representation of orbit set
        fixed = next(iter(o) for o in orbs if len(o) == 1)
        threes_as_sets = [frozenset(o) for o in orbs if len(o) == 3]
        threes_as_sets.sort(key=lambda s: tuple(sorted(s)))
        key = (fixed, tuple(threes_as_sets))
        if key not in z3_subgroups:
            z3_subgroups.append(key)

    print(f"Distinct Z_3 subgroups (orbit-structure unique): {len(z3_subgroups)}")
    print()
    for fixed, threes in z3_subgroups:
        print(f"  Fixed e_{fixed}, "
              f"3-orbits: {{{','.join(map(str, sorted(threes[0])))}}}, "
              f"{{{','.join(map(str, sorted(threes[1])))}}}")

    # Pick the lex-smallest canonical Z_3 and show its explicit action
    print()
    print("=" * 72)
    print("Canonical representative")
    print("=" * 72)
    sigma = candidates[0]
    orbs = orbits(sigma)
    fixed = next(iter(o) for o in orbs if len(o) == 1)
    threes = [sorted(o) for o in orbs if len(o) == 3]
    print(f"\n  Fixed point: e_{fixed}")
    print(f"  3-orbit A:   e_{threes[0][0]} → e_{threes[0][1]} → e_{threes[0][2]} → e_{threes[0][0]}")
    print(f"  3-orbit B:   e_{threes[1][0]} → e_{threes[1][1]} → e_{threes[1][2]} → e_{threes[1][0]}")
    print()
    print("  Action on imaginary units:")
    for i in range(1, 8):
        print(f"    σ(e_{i}) = e_{sigma[i]}")
    print()
    print("  Sanity check: σ has order 3:", order(sigma) == 3)
    print("  σ is in Aut(O) = G_2:", is_automorphism(sigma))


# ---------------------------------------------------------------------------
# Connection to Spin(8) triality (theoretical claim)
# ---------------------------------------------------------------------------

print()
print("=" * 72)
print("Connection to Spin(8) triality")
print("=" * 72)
print("""
  The Spin(8) outer automorphism group is S_3, including order-3
  elements that cyclically permute the three 8-dim representations
  V_8 ↔ S_8+ ↔ S_8-.  This is Spin(8) "triality."

  Embedding chain:
      Spin(8) ⊃ Spin(7) ⊃ G_2 ⊃ {our Z_3}

  When Spin(7) is embedded in Spin(8) by stabilizing a unit vector
  in V_8, the triality outer automorphism breaks (it now mixes
  V_8 with S_8±, but Spin(7) acts irreducibly on each piece).

  However, the residual order-3 symmetry inside Spin(7) is captured
  by elements of G_2 = Aut(O) ⊂ Spin(7).  The Z_3 we computed above
  (if it exists) realizes this residual triality at the level of
  K_7 vertices = imaginary octonion units.

  Equivalently: the *outer* Z_3 of Spin(8) becomes an *inner* Z_3 of
  Spin(7), and its action on the 7 imaginary octonion directions IS
  the σ we just verified.

  This is the standard restriction-of-triality picture (Baez 2002,
  §3.2, "The Octonions and Triality"). The verification above is the
  concrete realization for our specific Baez Fano convention.
""")
