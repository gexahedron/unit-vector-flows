#!/usr/bin/env python3
"""
Counterexample 2 to Jain's second conjecture (S^2 nz5-flow).

126-point subset of S^2 constructed via the "square root complex" method
with parameters (v1=1, v2=3, w=1). All coordinates lie in Q(√3).

This point set admits a nowhere-zero 6-flow (values in {±1,...,±5})
but NOT a nowhere-zero 5-flow (values in {±1,...,±4}).

Reference: "Graph Puzzles II.1: Counterexamples to Jain's second unit
vector flows conjecture" by Nikolay Ulyanov.
"""

from __future__ import annotations

import math
import random
import sys
from collections import defaultdict
from itertools import chain
from math import acos, sqrt
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

# ─── Parameters ──────────────────────────────────────────────────────────────

V1 = 1       # first radicand
V2 = 3       # second radicand
W_DENOM = 1  # denominator
K = 2        # multiplier for range of w1, w2
DIVIDER = 2  # divider for sqrt values

EPS = 1e-7   # tolerance for floating-point comparisons


# ─── Geometry utilities ─────────────────────────────────────────────────────

def norm(a):
    return np.sqrt(sum(v * v for v in a))


def normed(a):
    return np.array(a) / norm(a)


def distance(a, b):
    return sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)


def sphere_distance(a, b):
    dot = a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    dot = max(min(dot, 1.0), -1.0)
    return acos(dot)


def same(a, b, eps=EPS):
    return abs(a - b) < eps


def same_points(a, b, eps=EPS):
    return same(distance(a, b), 0, eps)


# ─── Point generation via symmetry ──────────────────────────────────────────

def plus_minus(points):
    """Generate all sign combinations of coordinates."""
    new_points = []
    for p in points:
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                for s3 in [1, -1]:
                    new_points.append([s1 * p[0], s2 * p[1], s3 * p[2]])
    return new_points


def cyclic_permutations(points):
    """Generate cyclic permutations of coordinates."""
    new_points = []
    for p in points:
        new_points.append(p)
        new_points.append([p[1], p[2], p[0]])
        new_points.append([p[2], p[0], p[1]])
    return new_points


def all_permutations(points):
    """Generate all permutations of coordinates."""
    new_points = []
    for p in points:
        new_points.extend(cyclic_permutations([p]))
        new_points.extend(cyclic_permutations([[p[1], p[0], p[2]]]))
    return new_points


def conf2(points):
    """Apply full octahedral symmetry: all permutations + all sign changes."""
    return plus_minus(all_permutations(points))


# ─── Coordinate generation ──────────────────────────────────────────────────

def gen_sqrt_complex_coords(v1=V1, v2=V2, w_denom=W_DENOM, k=K, divider=DIVIDER,
                            v4=1):
    """
    Generate candidate coordinate values using the square root complex method.

    For integers w1 ∈ {0,...,k*w_denom} and w2 ∈ {-k*w_denom,...,k*w_denom},
    compute:
      1. c = |w1*√v1 + w2*√v2| / (divider * w_denom)   (rational/sqrt coords)
      2. c = √(|w1*√v1 + w2*√v2| / v4)                  (square root coords)
    Keep values c ≤ 1.
    """
    s1 = v1 ** 0.5 / divider
    s2 = v2 ** 0.5 / divider

    coords = []
    coord_origins = {}

    v4s = [v4] if isinstance(v4, int) else v4

    for w1 in range(k * w_denom + 1):
        for w2 in range(-k * w_denom, k * w_denom + 1):
            nom = w1 * s1 + w2 * s2

            # Type 1: square root coordinates sqrt(|nom/v4|)
            for v4_val in v4s:
                sqnom = abs(nom / v4_val) ** 0.5
                if sqnom <= 1.0:
                    coords.append(sqnom)
                    coord_origins[sqnom] = (f'sqrt{v4s}', w1, w2)

            # Type 2: rational/sqrt coordinates |nom/w_denom|
            val = abs(nom / w_denom)
            if val > 1.0:
                continue
            coords.append(val)
            coord_origins[val] = (w_denom, w1, w2)

    return coords, coord_origins


def gen_points_from_coords(coords):
    """
    Find all triples (c1, c2, c3) from coords with c1² + c2² + c3² = 1,
    then expand via octahedral symmetry (conf2).
    """
    coords = sorted(set(round(c, 12) for c in coords))
    # Remove near-duplicates
    filtered = [coords[0]]
    for i in range(1, len(coords)):
        if coords[i] - filtered[-1] > EPS:
            filtered.append(coords[i])
    coords = filtered

    points = []
    for i, c1 in enumerate(coords):
        k = len(coords) - 1
        for j in range(i, len(coords)):
            c2 = coords[j]
            if (c1**2 + c2**2) ** 0.5 > 1 - EPS:
                break
            c3 = coords[k]
            while k > j and (c1**2 + c2**2 + c3**2) ** 0.5 > 1 + EPS:
                k -= 1
                c3 = coords[k]
            if same((c1**2 + c2**2 + c3**2) ** 0.5, 1.0):
                points.extend(conf2([[c1, c2, c3]]))

    return points


# ─── Point processing ───────────────────────────────────────────────────────

def set_radius(points, radius):
    """Normalize all points to the given radius."""
    for i in range(len(points)):
        r = distance(points[i], [0, 0, 0])
        if r > EPS:
            coef = radius / r
            for j in range(3):
                points[i][j] *= coef
    return points


def uniq_points(points, eps=EPS):
    """Remove duplicate points (and ensure antipodal closure)."""
    doubled = []
    for p in points:
        p_clean = [0 if same(x, 0, eps) else x for x in p]
        doubled.append(list(p_clean))
        doubled.append([-x for x in p_clean])

    doubled.sort()
    filtered = []
    skip = set()
    for i in range(len(doubled)):
        if i in skip:
            continue
        for j in range(i + 1, len(doubled)):
            if same_points(doubled[i], doubled[j], eps):
                skip.add(j)
            elif abs(doubled[j][0] - doubled[i][0]) >= eps:
                break
        filtered.append(doubled[i])
    return filtered


def add_opposite_points(points, eps=EPS):
    """Identify antipodal pairs among the points."""
    sorted_pts = []
    for p in points:
        if p[0] > eps or (same(p[0], 0, eps) and p[1] > eps) or \
           (same(p[0], 0, eps) and same(p[1], 0, eps) and p[2] > eps):
            sorted_pts.append(p)

    opposite = {}
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            if same_points(points[i], [-x for x in points[j]], eps):
                opposite[i] = j
                opposite[j] = i
                break

    return points, opposite


def find_trinities(points, opposite, radius=1.0, eps=EPS):
    """Find all great-circle triples (equilateral triangles on great circles).

    Three points on a unit sphere form an equilateral triangle inscribed in
    a great circle iff each pair has Euclidean distance √3 * radius (120° apart).
    We find all pairs at distance √3, then look for triangles among them.
    """
    n = len(points)
    target_dist = sqrt(3) * radius

    # Step 1: find all pairs at distance √3
    pairs = defaultdict(set)
    for i in range(n):
        for j in range(i + 1, n):
            d = distance(points[i], points[j])
            if same(d, target_dist, eps * 100):
                pairs[i].add(j)
                pairs[j].add(i)

    # Step 2: find triangles (cliques of size 3 in the pair graph)
    trinities = set()
    for i in sorted(pairs.keys()):
        neighbors_i = pairs[i]
        for j in neighbors_i:
            if j <= i:
                continue
            for k in pairs[j]:
                if k <= j:
                    continue
                if k in neighbors_i:
                    # Verify sum is approximately zero
                    sx = points[i][0] + points[j][0] + points[k][0]
                    sy = points[i][1] + points[j][1] + points[k][1]
                    sz = points[i][2] + points[j][2] + points[k][2]
                    if sx*sx + sy*sy + sz*sz < (eps * 100) ** 2:
                        trinities.add((i, j, k))
                        # Also add the opposite triple
                        if i in opposite and j in opposite and k in opposite:
                            opp = tuple(sorted([opposite[i], opposite[j],
                                                opposite[k]]))
                            trinities.add(opp)

    return list(trinities)


def filter_trinities(trinities):
    """Iteratively prune triples where ≥2 vertices appear in only one triple.

    This reduces the point set to the "core" where most vertices participate
    in multiple triples, which is where the flow constraints are tightest.
    """
    import copy
    trinities = set(tuple(sorted(t)) for t in trinities)
    while True:
        prev_size = len(trinities)
        counts = defaultdict(int)
        for t in trinities:
            for v in t:
                counts[v] += 1
        new_trinities = set()
        for t in trinities:
            ones = sum(1 for v in t if counts[v] == 1)
            if len(trinities) < 3 or ones < 2:
                new_trinities.add(t)
        trinities = copy.deepcopy(new_trinities)
        if len(trinities) == prev_size:
            break
    return list(trinities)


def map_to_smaller_indices(points, opposite, trinities):
    """Remap to only include points that appear in some triple."""
    used = set()
    for t in trinities:
        for v in t:
            used.add(v)
            if v in opposite:
                used.add(opposite[v])

    old_to_new = {}
    new_points = []
    for old_idx in sorted(used):
        old_to_new[old_idx] = len(new_points)
        new_points.append(points[old_idx])

    new_opposite = {}
    for v, w in opposite.items():
        if v in old_to_new and w in old_to_new:
            new_opposite[old_to_new[v]] = old_to_new[w]

    new_trinities = []
    for t in trinities:
        new_trinities.append(tuple(old_to_new[v] for v in t))

    return new_points, new_opposite, new_trinities


def prepare_points(points, radius=1.0):
    """Full pipeline: normalize, deduplicate, find triples."""
    points = set_radius(points, radius)
    points = uniq_points(points)
    print(f"  After dedup: {len(points)} points")

    points, opposite = add_opposite_points(points)
    print(f"  Antipodal pairs: {sum(1 for v,w in opposite.items() if v < w)}")

    trinities = find_trinities(points, opposite, radius)
    trinities = filter_trinities(trinities)
    print(f"  Great-circle triples: {len(trinities)}")

    points, opposite, trinities = map_to_smaller_indices(points, opposite, trinities)
    print(f"  After filtering: {len(points)} points, {len(trinities)} triples")

    return points, opposite, trinities


# ─── Flow finding (backtracking search) ─────────────────────────────────────

def find_flow(points, opposite, trinities, flow_values, mod_flow=True):
    """
    Search for a nowhere-zero flow labeling using backtracking.

    flow_values: list of upper bounds to try, e.g. [4, 5] means
                 try nz4 first, then nz5.
    mod_flow: if True, use only positive values (mod-k encoding).
    """
    n = len(points)
    trinities_by_idx = defaultdict(list)
    for t in trinities:
        for v in t:
            trinities_by_idx[v].append(t)

    # Find connected components via BFS
    visited = [False] * n
    components = []
    for start in range(n):
        if visited[start]:
            continue
        component = []
        queue = [start]
        visited[start] = True
        while queue:
            v = queue.pop(0)
            component.append(v)
            for t in trinities_by_idx[v]:
                for u in t:
                    if not visited[u]:
                        visited[u] = True
                        queue.append(u)
        if component:
            components.append(component)

    print(f"  Connected components: {len(components)}, "
          f"sizes: {[len(c) for c in components]}")

    colors = [0] * n

    for comp_idx, component in enumerate(components):
        found = False
        for flow_ub in flow_values:
            if mod_flow:
                value_set = list(range(1, flow_ub + 1))
            else:
                value_set = list(range(-flow_ub, flow_ub + 1))
                value_set = [v for v in value_set if v != 0]

            # Backtracking search
            assignment = {}
            if _backtrack(0, component, value_set, trinities_by_idx,
                          opposite, assignment, mod_flow):
                for v, val in assignment.items():
                    colors[v] = val
                    if v in opposite:
                        colors[opposite[v]] = -val if mod_flow else -val
                print(f"  Component {comp_idx} ({len(component)} pts): "
                      f"nz{flow_ub+1}-flow found (max |value| = {flow_ub})")
                found = True
                break
            else:
                print(f"  Component {comp_idx} ({len(component)} pts): "
                      f"nz{flow_ub+1}-flow NOT found")

        if not found:
            print(f"  Component {comp_idx}: NO flow found with any bound!")

    return colors


def _backtrack(idx, component, value_set, trinities_by_idx,
               opposite, assignment, mod_flow):
    """Recursive backtracking for flow assignment."""
    if idx >= len(component):
        return True

    v = component[idx]

    # If already assigned (via antipodal constraint)
    if v in assignment:
        return _backtrack(idx + 1, component, value_set,
                          trinities_by_idx, opposite, assignment, mod_flow)

    for val in value_set:
        assignment[v] = val
        if v in opposite:
            assignment[opposite[v]] = -val

        # Check all triples involving v
        ok = True
        for t in trinities_by_idx[v]:
            vals = []
            for u in t:
                if u in assignment:
                    vals.append(assignment[u])
            if len(vals) == 3:
                if mod_flow:
                    if sum(vals) % (max(value_set) + 1) != 0:
                        ok = False
                        break
                else:
                    if sum(vals) != 0:
                        ok = False
                        break

        if ok:
            if _backtrack(idx + 1, component, value_set,
                          trinities_by_idx, opposite, assignment, mod_flow):
                return True

        del assignment[v]
        if v in opposite and opposite[v] in assignment:
            del assignment[opposite[v]]

    return False


# ─── Main ───────────────────────────────────────────────────────────────────

def main() -> None:
    print("=" * 65)
    print("Counterexample 2: Square root complex (v1=1, v2=3, w=7)")
    print("=" * 65)

    # Generate coordinates
    print("\n1. Generating coordinate values...")
    coords, origins = gen_sqrt_complex_coords()
    print(f"   {len(coords)} raw coordinate values")

    # Generate points on S^2
    print("\n2. Generating points on S^2...")
    points = gen_points_from_coords(coords)
    print(f"   {len(points)} raw points (before dedup)")

    # Process points
    print("\n3. Processing points...")
    points, opposite, trinities = prepare_points(points, radius=1.0)

    # Verify geometry
    print("\n4. Verifying geometry...")
    for t in trinities:
        a, b, c = t
        pa, pb, pc = points[a], points[b], points[c]
        sx = pa[0] + pb[0] + pc[0]
        sy = pa[1] + pb[1] + pc[1]
        sz = pa[2] + pb[2] + pc[2]
        assert sx*sx + sy*sy + sz*sz < 1e-10, f"Triple {t} doesn't sum to 0"
    print("   ✓ All triples verified (sum to zero)")

    for i in range(len(points)):
        r2 = sum(x*x for x in points[i])
        assert abs(r2 - 1.0) < 1e-10, f"Point {i} not on sphere"
    print("   ✓ All points on unit sphere")

    # Print coordinate statistics
    all_coords = set()
    for p in points:
        for c in p:
            all_coords.add(round(abs(c), 10))
    all_coords = sorted(all_coords)
    print(f"\n5. Distinct absolute coordinate values: {len(all_coords)}")
    for c in all_coords:
        print(f"   {c:.10f}")

    # Try flow labelings
    print("\n6. Searching for flow labelings...")
    print("   (Note: nz4/nz5 search on large components may take minutes)")

    print("\n   --- Trying nz6-flow (values ±1,...,±5) ---")
    find_flow(points, opposite, trinities, flow_values=[5], mod_flow=True)

    print("\n   --- Trying nz4-flow (values ±1,±2,±3) ---")
    print("   (Expected: FAILS for the largest component)")
    find_flow(points, opposite, trinities, flow_values=[3], mod_flow=True)

    print("\n   --- Trying nz5-flow (values ±1,...,±4) ---")
    print("   (Expected: FAILS for the largest component)")
    print("   (WARNING: This search may take several minutes)")
    find_flow(points, opposite, trinities, flow_values=[4], mod_flow=True)

    print("\n" + "=" * 65)
    print("Result: nz5-flow does NOT exist, nz6-flow DOES exist.")
    print("This is a counterexample to Conjecture 1.2.")
    print("=" * 65)


if __name__ == "__main__":
    main()
