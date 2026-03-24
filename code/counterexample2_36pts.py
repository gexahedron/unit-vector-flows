#!/usr/bin/env python3
"""
Counterexample 2 to Jain's second conjecture (S^2 nz5-flow conjecture).

36-point subset of S^2 constructed via the "square root arithmetic" method
with parameters v1=1, v2=3, w=2.  All coordinates lie in
Q(sqrt(3), sqrt(sqrt(3)-1)), a degree-4 extension of Q with denominator 2.

This point set admits a nowhere-zero 6-flow (values in {±1,...,±5})
but NOT a nowhere-zero 5-flow (values in {±1,...,±4}).

Reference: "Graph Puzzles II.1: Counterexamples to Jain's second unit
vector flows conjecture" by Nikolay Ulyanov, 2026.
"""

from __future__ import annotations

import math
from math import sqrt
from typing import Dict, List, Tuple

# ─── Algebraic constants ─────────────────────────────────────────────────────
# The 7 distinct absolute coordinate values:
#   0, (2-√3)/2, (√3-1)/2, 1/2, √(√3-1), √3/2, 1

_s3 = sqrt(3)          # √3
_sv = sqrt(_s3 - 1)    # √(√3-1)
_a  = (2 - _s3) / 2   # (2-√3)/2  ≈ 0.1340
_b  = (_s3 - 1) / 2   # (√3-1)/2  ≈ 0.3660
_c  = 1 / 2            # 1/2
_d  = _sv              # √(√3-1)   ≈ 0.8556
_e  = _s3 / 2          # √3/2      ≈ 0.8660
_f  = 1.0              # 1

# ─── 36 vertices on S^2 ─────────────────────────────────────────────────────
# Coordinates are exact elements of Q(√3, √(√3-1)) with denominator 2.
# Antipodal pairs: vertex i is antipodal to vertex 35-i.
#
# Notation for the 7 distinct |coord| values:
#   a = (2-√3)/2,  b = (√3-1)/2,  c = 1/2,
#   d = √(√3-1),   e = √3/2,      f = 1

POINTS: List[Tuple[float, float, float]] = [
    (-_f,  0,   0  ),  #  0: (-f,  0,  0)
    (-_e,  0,  -_c ),  #  1: (-e,  0, -c)
    (-_e,  0,  +_c ),  #  2: (-e,  0, +c)
    (-_c,  0,  -_e ),  #  3: (-c,  0, -e)
    (-_c,  0,  +_e ),  #  4: (-c,  0, +e)
    (-_c, -_d, -_a ),  #  5: (-c, -d, -a)
    (-_c, -_d, +_a ),  #  6: (-c, -d, +a)
    (-_c, +_d, -_a ),  #  7: (-c, +d, -a)
    (-_c, +_d, +_a ),  #  8: (-c, +d, +a)
    (-_b, -_d, -_b ),  #  9: (-b, -d, -b)
    (-_b, -_d, +_b ),  # 10: (-b, -d, +b)
    (-_b, +_d, -_b ),  # 11: (-b, +d, -b)
    (-_b, +_d, +_b ),  # 12: (-b, +d, +b)
    (-_a, -_d, -_c ),  # 13: (-a, -d, -c)
    (-_a, -_d, +_c ),  # 14: (-a, -d, +c)
    (-_a, +_d, -_c ),  # 15: (-a, +d, -c)
    (-_a, +_d, +_c ),  # 16: (-a, +d, +c)
    ( 0,   0,  -_f ),  # 17: ( 0,  0, -f)
    ( 0,   0,  +_f ),  # 18: ( 0,  0, +f)
    (+_a, -_d, -_c ),  # 19: (+a, -d, -c)
    (+_a, -_d, +_c ),  # 20: (+a, -d, +c)
    (+_a, +_d, -_c ),  # 21: (+a, +d, -c)
    (+_a, +_d, +_c ),  # 22: (+a, +d, +c)
    (+_b, -_d, -_b ),  # 23: (+b, -d, -b)
    (+_b, -_d, +_b ),  # 24: (+b, -d, +b)
    (+_b, +_d, -_b ),  # 25: (+b, +d, -b)
    (+_b, +_d, +_b ),  # 26: (+b, +d, +b)
    (+_c, -_d, -_a ),  # 27: (+c, -d, -a)
    (+_c, -_d, +_a ),  # 28: (+c, -d, +a)
    (+_c, +_d, -_a ),  # 29: (+c, +d, -a)
    (+_c, +_d, +_a ),  # 30: (+c, +d, +a)
    (+_c,  0,  -_e ),  # 31: (+c,  0, -e)
    (+_c,  0,  +_e ),  # 32: (+c,  0, +e)
    (+_e,  0,  -_c ),  # 33: (+e,  0, -c)
    (+_e,  0,  +_c ),  # 34: (+e,  0, +c)
    (+_f,  0,   0  ),  # 35: (+f,  0,  0)
]

# ─── 13 great-circle triples ─────────────────────────────────────────────────
# Each triple (a, b, c) satisfies POINTS[a] + POINTS[b] + POINTS[c] = (0,0,0),
# i.e. the three points are equidistant (120° apart) on a great circle.

TRIPLES: List[Tuple[int, int, int]] = [
    ( 0, 31, 32),
    ( 2, 17, 34),
    ( 5,  8, 35),
    ( 5, 11, 34),
    ( 6,  7, 35),
    ( 7,  9, 34),
    ( 8, 10, 33),
    ( 9, 15, 32),
    (10, 16, 31),
    (11, 13, 32),
    (12, 14, 31),
    (14, 17, 22),
    (16, 17, 20),
]

# ─── 18 antipodal pairs ──────────────────────────────────────────────────────
# Vertex i is antipodal to vertex 35-i.

ANTIPODE: Dict[int, int] = {i: 35 - i for i in range(36)}


# ─── Verification ────────────────────────────────────────────────────────────

def verify_on_sphere(tol: float = 1e-10) -> bool:
    """Verify all 36 points lie on the unit sphere."""
    for i, (x, y, z) in enumerate(POINTS):
        r2 = x*x + y*y + z*z
        if abs(r2 - 1.0) > tol:
            print(f"Point {i} has |r|^2 = {r2:.15f}")
            return False
    return True


def verify_antipodal(tol: float = 1e-10) -> bool:
    """Verify antipodal pairs: POINTS[i] + POINTS[35-i] = 0."""
    for i in range(18):
        j = ANTIPODE[i]
        pi, pj = POINTS[i], POINTS[j]
        err = (pi[0]+pj[0])**2 + (pi[1]+pj[1])**2 + (pi[2]+pj[2])**2
        if err > tol:
            print(f"Points {i} and {j} are not antipodal (err={err:.2e})")
            return False
    return True


def verify_triples(tol: float = 1e-10) -> bool:
    """Verify each triple sums to zero (equidistant on a great circle)."""
    for a, b, c in TRIPLES:
        pa, pb, pc = POINTS[a], POINTS[b], POINTS[c]
        sx = pa[0] + pb[0] + pc[0]
        sy = pa[1] + pb[1] + pc[1]
        sz = pa[2] + pb[2] + pc[2]
        if sx*sx + sy*sy + sz*sz > tol:
            print(f"Triple ({a},{b},{c}): sum = ({sx:.2e},{sy:.2e},{sz:.2e})")
            return False
        # Check coplanarity with origin (det = 0)
        det = (pa[0]*(pb[1]*pc[2] - pb[2]*pc[1])
             - pa[1]*(pb[0]*pc[2] - pb[2]*pc[0])
             + pa[2]*(pb[0]*pc[1] - pb[1]*pc[0]))
        if abs(det) > tol:
            print(f"Triple ({a},{b},{c}): det = {det:.2e} (not coplanar with origin)")
            return False
    return True


def print_coordinate_summary() -> None:
    """Print the 7 distinct absolute coordinate values."""
    abs_coords = sorted(set(round(abs(c), 10) for p in POINTS for c in p))
    names = {
        round(0.0,   10): '0',
        round(_a,    10): '(2-√3)/2',
        round(_b,    10): '(√3-1)/2',
        round(_c,    10): '1/2',
        round(_d,    10): '√(√3-1)',
        round(_e,    10): '√3/2',
        round(1.0,   10): '1',
    }
    print(f"  Distinct |coord| values ({len(abs_coords)}):")
    for v in abs_coords:
        print(f"    {v:.10f}  =  {names.get(round(v,10), '?')}")


def main() -> None:
    print("=" * 60)
    print("Counterexample 2: 36-point subset of S^2")
    print("Parameters: v1=1, v2=3, w=2")
    print("=" * 60)
    print(f"\n  Points:         {len(POINTS)}")
    print(f"  Triples:        {len(TRIPLES)}")
    print(f"  Antipodal pairs:{len(ANTIPODE) // 2}")

    assert verify_on_sphere(),  "FAIL: points not on sphere"
    print("\n  ✓ All 36 points lie on S^2")

    assert verify_antipodal(),  "FAIL: antipodal pairs incorrect"
    print("  ✓ All 18 antipodal pairs verified")

    assert verify_triples(),    "FAIL: triples invalid"
    print("  ✓ All 13 great-circle triples verified (sum to zero)")

    print()
    print_coordinate_summary()
    print()
    print("  Algebraic field: Q(√3, √(√3-1))  [degree 4 over Q, denominator 2]")
    print()
    print("  To verify nz5-flow impossibility and nz6-flow existence, run:")
    print("    python3 sat_verify_36pts.py")


if __name__ == "__main__":
    main()
