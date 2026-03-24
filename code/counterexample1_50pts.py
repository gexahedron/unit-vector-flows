#!/usr/bin/env python3
"""
Counterexample 1 to Jain's second conjecture (S^2 nz5-flow).

50-point subset of S^2 obtained by expanding the icosidodecahedron.
This point set admits a nowhere-zero 6-flow (values in {±1,...,±5})
but NOT a nowhere-zero 5-flow (values in {±1,...,±4}).

The unsatisfiability of the nz5-flow labeling is verified via SAT solving.

Reference: "Graph Puzzles II.1: Counterexamples to Jain's second unit
vector flows conjecture" by Nikolay Ulyanov.
"""

from __future__ import annotations

import math
from math import sqrt
from typing import Dict, List, Tuple

# Golden ratio
PHI = (1 + sqrt(5)) / 2

# Scaling radius: R = 2*phi
R = 2 * PHI

# Auxiliary values used in coordinates
# x = 2 / 5^(1/4), y = x * phi
X_VAL = 2 / (5 ** 0.25)
Y_VAL = X_VAL * PHI

# ─── 50 vertices (scaled coordinates, sphere radius = R = 2φ) ───────────────

COORDS: Dict[int, Tuple[str, str, str]] = {
    0:  ("y+phi-1",     "-(x-1)",       "-phi"),
    1:  ("-(y-1)",      "x+phi",        "phi-1"),
    2:  ("-phi",        "-(phi+1)",     "1"),
    3:  ("phi+1",       "-1",           "-phi"),
    4:  ("0",           "0",            "2*phi"),
    5:  ("-(phi+1)",    "1",            "-phi"),
    6:  ("phi",         "phi+1",        "1"),
    7:  ("y-phi+1",     "-(x+1)",       "-phi"),
    8:  ("-(y+1)",      "-(phi-x)",     "phi-1"),
    9:  ("-(y-phi+1)",  "x+1",          "-phi"),
    10: ("phi",         "phi+1",        "-1"),
    11: ("y-1",         "-(x+phi)",     "-(phi-1)"),
    12: ("-(y+phi-1)",  "x-1",          "phi"),
    13: ("phi",         "-(phi+1)",     "-1"),
    14: ("-2*phi",      "0",            "0"),
    15: ("phi+1",       "-1",           "phi"),
    16: ("-1",          "-phi",         "-(phi+1)"),
    17: ("-phi",        "phi+1",        "1"),
    18: ("y+phi-1",     "-(x-1)",       "phi"),
    19: ("0",           "0",            "-2*phi"),
    20: ("y-phi+1",     "-(x+1)",       "phi"),
    21: ("-(y-phi+1)",  "x+1",          "phi"),
    22: ("y-1",         "-(x+phi)",     "phi-1"),
    23: ("-(y+phi-1)",  "x-1",          "-phi"),
    24: ("-(phi+1)",    "1",            "phi"),
    25: ("-1",          "-phi",         "phi+1"),
    26: ("-phi",        "phi+1",        "-1"),
    27: ("1",           "-phi",         "phi+1"),
    28: ("-(phi+1)",    "-1",           "-phi"),
    29: ("phi+1",       "1",            "-phi"),
    30: ("-1",          "phi",          "phi+1"),
    31: ("-phi",        "-(phi+1)",     "-1"),
    32: ("2*phi",       "0",            "0"),
    33: ("y+1",         "phi-x",        "phi-1"),
    34: ("-y",          "x",            "2"),
    35: ("1",           "phi",          "phi+1"),
    36: ("0",           "2*phi",        "0"),
    37: ("1",           "phi",          "-(phi+1)"),
    38: ("0",           "-2*phi",       "0"),
    39: ("1",           "-phi",         "-(phi+1)"),
    40: ("-(phi+1)",    "-1",           "phi"),
    41: ("y+1",         "phi-x",        "-(phi-1)"),
    42: ("-y",          "x",            "-2"),
    43: ("-(y+1)",      "-(phi-x)",     "-(phi-1)"),
    44: ("-(y-1)",      "x+phi",        "-(phi-1)"),
    45: ("y",           "-x",           "2"),
    46: ("y",           "-x",           "-2"),
    47: ("-1",          "phi",          "-(phi+1)"),
    48: ("phi",         "-(phi+1)",     "1"),
    49: ("phi+1",       "1",            "phi"),
}

# ─── 40 great-circle triples ────────────────────────────────────────────────

TRIPLES: List[Tuple[int, int, int]] = [
    (41, 25, 42),
    (18, 19, 12),
    (49, 19, 40),
    (22, 37, 34),
    (10, 11, 12),
    (0,  4,  23),
    (29, 4,  28),
    (48, 37, 24),
    (33, 9,  2),
    (20, 19, 21),
    (6,  39, 40),
    (11, 35, 42),
    (7,  4,  9),
    (6,  7,  8),
    (29, 30, 31),
    (41, 21, 31),
    (35, 38, 47),
    (32, 26, 2),
    (0,  1,  2),
    (10, 27, 28),
    (45, 16, 1),
    (10, 48, 14),
    (33, 16, 34),
    (6,  22, 23),
    (3,  25, 26),
    (18, 44, 31),
    (15, 16, 17),
    (10, 20, 43),
    (46, 35, 43),
    (6,  13, 14),
    (45, 37, 8),
    (37, 38, 30),
    (27, 36, 16),
    (13, 35, 5),
    (49, 47, 2),
    (15, 19, 24),
    (3,  4,  5),
    (32, 17, 31),
    (39, 36, 25),
    (46, 25, 44),
]


def eval_coord(expr: str) -> float:
    """Evaluate a symbolic coordinate expression to a float."""
    phi = PHI
    x = X_VAL
    y = Y_VAL
    expr = expr.strip()
    # Handle negation of parenthesized expressions
    result = eval(expr, {"__builtins__": {}},
                  {"phi": phi, "x": x, "y": y})
    return float(result)


def get_numeric_coords() -> List[Tuple[float, float, float]]:
    """Return the 50 points as numeric (x, y, z) tuples on S^2."""
    points = []
    for i in range(len(COORDS)):
        sx, sy, sz = COORDS[i]
        px = eval_coord(sx) / R
        py = eval_coord(sy) / R
        pz = eval_coord(sz) / R
        points.append((px, py, pz))
    return points


def verify_on_sphere(points: List[Tuple[float, float, float]],
                     tol: float = 1e-10) -> bool:
    """Verify all points lie on the unit sphere."""
    for i, (x, y, z) in enumerate(points):
        r2 = x*x + y*y + z*z
        if abs(r2 - 1.0) > tol:
            print(f"Point {i} has |r|^2 = {r2}, off by {abs(r2-1.0)}")
            return False
    return True


def verify_antipodal(points: List[Tuple[float, float, float]],
                     tol: float = 1e-10) -> Dict[int, int]:
    """Find and verify antipodal pairs."""
    n = len(points)
    antipode: Dict[int, int] = {}
    for i in range(n):
        for j in range(i + 1, n):
            dx = points[i][0] + points[j][0]
            dy = points[i][1] + points[j][1]
            dz = points[i][2] + points[j][2]
            if dx*dx + dy*dy + dz*dz < tol:
                antipode[i] = j
                antipode[j] = i
    return antipode


def verify_triples(points: List[Tuple[float, float, float]],
                   triples: List[Tuple[int, int, int]],
                   tol: float = 1e-8) -> bool:
    """Verify each triple forms an equilateral triangle on a great circle."""
    for a, b, c in triples:
        pa, pb, pc = points[a], points[b], points[c]
        # Check sum is zero (equidistant on great circle)
        sx = pa[0] + pb[0] + pc[0]
        sy = pa[1] + pb[1] + pc[1]
        sz = pa[2] + pb[2] + pc[2]
        if sx*sx + sy*sy + sz*sz > tol:
            print(f"Triple ({a},{b},{c}): sum = ({sx},{sy},{sz})")
            return False
        # Check all three are on the same great circle (coplanar with origin)
        # det(pa, pb, pc) = 0
        det = (pa[0] * (pb[1]*pc[2] - pb[2]*pc[1])
             - pa[1] * (pb[0]*pc[2] - pb[2]*pc[0])
             + pa[2] * (pb[0]*pc[1] - pb[1]*pc[0]))
        if abs(det) > tol:
            print(f"Triple ({a},{b},{c}): det = {det} (not coplanar)")
            return False
    return True


def main() -> None:
    print("=" * 60)
    print("Counterexample 1: 50-point expansion of icosidodecahedron")
    print("=" * 60)

    points = get_numeric_coords()
    print(f"Number of points: {len(points)}")
    print(f"Number of triples: {len(TRIPLES)}")

    # Verify geometry
    assert verify_on_sphere(points), "Some points are not on S^2!"
    print("✓ All 50 points lie on S^2")

    antipode = verify_antipodal(points)
    n_pairs = sum(1 for v, w in antipode.items() if v < w)
    print(f"✓ Found {n_pairs} antipodal pairs")

    assert verify_triples(points, TRIPLES), "Some triples are invalid!"
    print("✓ All 40 triples are valid great-circle equilateral triangles")

    # Print some sample coordinates
    print("\nSample points (first 5):")
    for i in range(5):
        sx, sy, sz = COORDS[i]
        px, py, pz = points[i]
        print(f"  v{i}: ({sx}, {sy}, {sz})")
        print(f"       = ({px:.10f}, {py:.10f}, {pz:.10f})")

    print("\nSample triples (first 5):")
    for a, b, c in TRIPLES[:5]:
        print(f"  ({a}, {b}, {c})")

    print("\nTo verify UNSAT (no nz5-flow labeling), run:")
    print("  python3 sat_verify_50pts.py")


if __name__ == "__main__":
    main()
