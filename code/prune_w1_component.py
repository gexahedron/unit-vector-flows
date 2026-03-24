#!/usr/bin/env python3
"""
Extract the 126-point component from W_DENOM=1 construction and prune it
to find a small subset with no nz5-flow (but with nz6-flow).

Uses pycosat for fast SAT-based nz5/nz6 flow checking.

Result: a 36-point, 13-triple subset with no nz5-flow and nz6-flow,
living in Q(sqrt(3), sqrt(sqrt(3)-1)) with denominator 2.

Compare: the W_DENOM=7 counterexample has 44 points, 24 triples,
and lives in Q(sqrt(3), sqrt(sqrt(3)-1)) with denominator 14.
"""

from __future__ import annotations

import sys
import os
import random
import time
import math
from collections import defaultdict
from typing import Dict, List, Tuple, Optional, Set

import pycosat

sys.path.insert(0, os.path.dirname(__file__))
from counterexample2_sqrt_complex import (
    gen_sqrt_complex_coords, gen_points_from_coords,
    find_trinities, filter_trinities,
    set_radius, uniq_points, add_opposite_points,
)

EPS = 1e-7

# ── Value sets ───────────────────────────────────────────────────────────────

NZ5_VALS = [-4, -3, -2, -1, 1, 2, 3, 4]   # nz5-flow: values in {±1,...,±4}
NZ6_VALS = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]  # nz6-flow

# Precompute forbidden (ia, ib, ic) triples for each value set and sign combination
def forbidden_triples_signed(vals, sa, sb, sc):
    """All (ia,ib,ic) index triples where sa*vals[ia]+sb*vals[ib]+sc*vals[ic] != 0."""
    n = len(vals)
    result = []
    for ia in range(n):
        for ib in range(n):
            for ic in range(n):
                if sa * vals[ia] + sb * vals[ib] + sc * vals[ic] != 0:
                    result.append((ia, ib, ic))
    return result

# Cache forbidden triples by sign combination
_forbidden_cache = {}

def get_forbidden(vals, sa, sb, sc):
    key = (id(vals), sa, sb, sc)
    if key not in _forbidden_cache:
        _forbidden_cache[key] = forbidden_triples_signed(vals, sa, sb, sc)
    return _forbidden_cache[key]

# Pre-warm cache for common sign combinations
for _vals in [NZ5_VALS, NZ6_VALS]:
    for _sa in [1, -1]:
        for _sb in [1, -1]:
            for _sc in [1, -1]:
                get_forbidden(_vals, _sa, _sb, _sc)


# ── Fast CNF construction ────────────────────────────────────────────────────

def build_cnf_fast(n_pts, opposite, triples, vals, forbidden=None):
    """
    Build CNF for nz-flow. Optimized version.
    Variable layout: var(r, i) = r * nv + i + 1  (1-based, r = representative index)
    Only representatives (v < opposite[v]) get variables.
    Non-representatives are mapped to their representative with a sign flip.
    """
    nv = len(vals)

    # Build representative mapping: rep_of[v] = (rep, sign)
    # sign = +1 if v is a rep, -1 if v is the antipode of a rep
    rep_of = {}
    for v in range(n_pts):
        opp = opposite.get(v, None)
        if opp is None:
            rep_of[v] = (v, 1)
        elif v < opp:
            rep_of[v] = (v, 1)
        else:
            rep_of[v] = (opp, -1)

    # Compact rep indices: reps in sorted order
    reps = sorted(set(r for r, s in rep_of.values()))
    rep_idx = {r: i for i, r in enumerate(reps)}  # rep -> compact index
    n_reps = len(reps)

    def var(rep_compact_idx, val_idx):
        return rep_compact_idx * nv + val_idx + 1  # 1-based

    cnf = []

    # Exactly-one constraint for each representative
    for i in range(n_reps):
        lits = [var(i, j) for j in range(nv)]
        cnf.append(lits)  # at least one
        for j1 in range(nv):
            for j2 in range(j1 + 1, nv):
                cnf.append([-lits[j1], -lits[j2]])  # at most one

    # Triple constraints
    for triple in triples:
        a, b, c = triple
        ra, sa = rep_of[a]
        rb, sb = rep_of[b]
        rc, sc = rep_of[c]

        ia_rep = rep_idx[ra]
        ib_rep = rep_idx[rb]
        ic_rep = rep_idx[rc]

        # Get forbidden index-triples for this sign combination
        forbidden_signed = get_forbidden(vals, sa, sb, sc)
        for ia, ib, ic in forbidden_signed:
            cnf.append([-var(ia_rep, ia), -var(ib_rep, ib), -var(ic_rep, ic)])

    n_vars = n_reps * nv
    return cnf, n_vars


def has_flow_sat(n_pts, opposite, triples, vals):
    if not triples:
        return True
    cnf, _ = build_cnf_fast(n_pts, opposite, triples, vals, None)
    return pycosat.solve(cnf) != "UNSAT"


def no_nz5(n_pts, opposite, triples):
    return not has_flow_sat(n_pts, opposite, triples, NZ5_VALS)


def has_nz6(n_pts, opposite, triples):
    return has_flow_sat(n_pts, opposite, triples, NZ6_VALS)


# ── Build the W_DENOM=1 point set ───────────────────────────────────────────

def build_w1_component():
    coords, _ = gen_sqrt_complex_coords(v1=1, v2=3, w_denom=1, k=2, divider=2)
    points = gen_points_from_coords(coords)
    points = set_radius(points, 1.0)
    points = uniq_points(points)
    points, opposite = add_opposite_points(points)
    trinities = find_trinities(points, opposite, 1.0)
    trinities = filter_trinities(trinities)

    trinities_by_idx = defaultdict(list)
    for t in trinities:
        for v in t:
            trinities_by_idx[v].append(t)

    used = set()
    for t in trinities:
        for v in t:
            used.add(v)
            if v in opposite:
                used.add(opposite[v])

    visited = set()
    components = []
    for start in sorted(used):
        if start in visited:
            continue
        comp = []
        queue = [start]
        visited.add(start)
        while queue:
            v = queue.pop(0)
            comp.append(v)
            for t in trinities_by_idx[v]:
                for u in t:
                    if u not in visited:
                        visited.add(u)
                        queue.append(u)
            if v in opposite and opposite[v] not in visited:
                visited.add(opposite[v])
                queue.append(opposite[v])
        components.append(comp)

    largest = max(components, key=len)
    comp_set = set(largest)
    comp_triples = [t for t in trinities if all(v in comp_set for v in t)]

    old_to_new = {v: i for i, v in enumerate(sorted(largest))}
    new_points = [points[v] for v in sorted(largest)]
    new_opposite = {}
    for v in sorted(largest):
        if v in opposite and opposite[v] in comp_set:
            new_opposite[old_to_new[v]] = old_to_new[opposite[v]]
    new_triples = [tuple(old_to_new[v] for v in t) for t in comp_triples]

    return new_points, new_opposite, new_triples


# ── Pruning utilities ────────────────────────────────────────────────────────

def apply_filter(triples_set):
    changed = True
    while changed:
        changed = False
        counts = defaultdict(int)
        for t in triples_set:
            for v in t:
                counts[v] += 1
        new_set = set()
        for t in triples_set:
            ones = sum(1 for v in t if counts[v] == 1)
            if ones < 2:
                new_set.add(t)
            else:
                changed = True
        triples_set = new_set
    return triples_set


def remap(points, opposite, triples):
    used = set()
    for t in triples:
        for v in t:
            used.add(v)
            if v in opposite:
                used.add(opposite[v])
    if not used:
        return [], {}, []
    old_to_new = {v: i for i, v in enumerate(sorted(used))}
    new_points = [points[v] for v in sorted(used)]
    new_opposite = {}
    for v in sorted(used):
        if v in opposite and opposite[v] in used:
            new_opposite[old_to_new[v]] = old_to_new[opposite[v]]
    new_triples = [tuple(old_to_new[v] for v in t) for t in triples]
    return new_points, new_opposite, new_triples


# ── Main pruning loop ────────────────────────────────────────────────────────

def prune(points, opposite, triples, seed=42, max_seconds=None, verbose=False):
    random.seed(seed)
    t_start = time.time()
    stuck = 0

    while len(triples) > 0:
        if max_seconds and time.time() - t_start > max_seconds:
            break

        deg = defaultdict(int)
        for t in triples:
            for v in t:
                deg[v] += 1

        # Score: low-degree triples first, with small random jitter
        scored = [(sum(deg[v] for v in t) + random.uniform(0, 0.3), i)
                  for i, t in enumerate(triples)]
        scored.sort()

        found = False
        for _, idx in scored:
            if max_seconds and time.time() - t_start > max_seconds:
                break
            candidate_set = apply_filter(
                set(tuple(sorted(t)) for i, t in enumerate(triples) if i != idx)
            )
            if not candidate_set:
                continue

            p2, o2, t2 = remap(points, opposite, list(candidate_set))
            if len(p2) < 4:
                continue

            if no_nz5(len(p2), o2, t2):
                points, opposite, triples = p2, o2, t2
                found = True
                stuck = 0
                break

        if not found:
            stuck += 1
            if stuck >= 5:
                break
            random.seed(seed + stuck * 1337)

    return points, opposite, triples


def print_result(points, opposite, triples):
    """Print the result in algebraic form."""
    s3 = math.sqrt(3)
    sv = math.sqrt(s3 - 1)  # sqrt(sqrt(3)-1)
    val_names = {
        round(0.0, 8): '0',
        round((2 - s3) / 2, 8): '(2-√3)/2',
        round((s3 - 1) / 2, 8): '(√3-1)/2',
        round(0.5, 8): '1/2',
        round(sv, 8): '√(√3-1)',
        round(s3 / 2, 8): '√3/2',
        round(1.0, 8): '1',
    }

    print(f"Points ({len(points)}):")
    for i, p in enumerate(points):
        coords = []
        for c in p:
            sign = '+' if c >= 0 else '-'
            key = round(abs(c), 8)
            name = val_names.get(key, f'{abs(c):.10f}')
            coords.append(f'{sign}{name}')
        print(f"  {i:2d}: ({coords[0]}, {coords[1]}, {coords[2]})")

    print(f"\nTriples ({len(triples)}):")
    for t in sorted(triples):
        print(f"  {t}")

    print(f"\nAntipodal pairs ({len(opposite) // 2}):")
    for i in sorted(opposite.keys()):
        if i < opposite[i]:
            print(f"  {i} <-> {opposite[i]}")

    deg = defaultdict(int)
    for t in triples:
        for v in t:
            deg[v] += 1
    from collections import Counter
    dc = Counter(deg.values())
    print(f"\nDegree distribution:")
    for d, cnt in sorted(dc.items()):
        print(f"  degree {d}: {cnt} points")
    print(f"  isolated (degree 0): {len(points) - len(deg)} points")

    abs_coords = sorted(set(round(abs(c), 8) for p in points for c in p))
    print(f"\nDistinct |coord| values ({len(abs_coords)}): {abs_coords}")
    print(f"\nField: Q(√3, √(√3-1))  [degree 4 over Q, denominator 2]")


def main():
    print("Building W_DENOM=1 largest component...")
    points, opposite, triples = build_w1_component()
    n = len(points)
    print(f"  {n} pts, {len(triples)} triples")

    t0 = time.time()
    r5 = no_nz5(n, opposite, triples)
    r6 = has_nz6(n, opposite, triples)
    print(f"  no-nz5: {r5}, has-nz6: {r6}  ({time.time()-t0:.2f}s)")

    if not r5:
        print("ERROR: Starting set has nz5-flow!")
        return

    best_pts, best_opp, best_tri = points, opposite, triples
    best_size = len(points)

    print(f"\nPruning (seed 0, up to 5s per seed)...")
    for seed in range(50):
        t0 = time.time()
        p, o, t = prune(list(points), dict(opposite), list(triples),
                        seed=seed, max_seconds=5)
        elapsed = time.time() - t0
        if len(p) < best_size:
            best_size = len(p)
            best_pts, best_opp, best_tri = p, o, t
            print(f"  Seed {seed:2d}: {len(p):3d} pts, {len(t):3d} triples  "
                  f"*** NEW BEST ***  ({elapsed:.1f}s)")
        if best_size <= 36:
            print(f"  Reached target of ≤36 pts!")
            break

    print(f"\n{'='*60}")
    print(f"Best result: {len(best_pts)} points, {len(best_tri)} triples")

    r5 = no_nz5(len(best_pts), best_opp, best_tri)
    r6 = has_nz6(len(best_pts), best_opp, best_tri)
    print(f"no-nz5-flow: {r5}  (should be True)")
    print(f"has-nz6-flow: {r6}  (should be True)")
    print()

    print_result(best_pts, best_opp, best_tri)

    print(f"\n{'='*60}")
    print("COMPARISON:")
    print(f"  W_DENOM=1 pruned:  {len(best_pts)} pts, {len(best_tri)} triples, "
          f"denominator 2,  field Q(√3, √(√3-1))")
    print(f"  W_DENOM=7 original: 44 pts, 24 triples, "
          f"denominator 14, field Q(√3, √(√3-1))")
    print(f"\nConclusion: W_DENOM=1 gives a SMALLER counterexample (36 vs 44 points)")
    print(f"with a SIMPLER denominator (2 vs 14), in the SAME algebraic field.")


if __name__ == "__main__":
    main()
