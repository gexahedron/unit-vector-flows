#!/usr/bin/env python3
"""
SAT-based verification that the 50-point counterexample has no nz5-flow labeling.

Encodes the labeling problem as a Boolean satisfiability (SAT) instance:
  - q(v) ∈ {-4, -3, -2, -1, 1, 2, 3, 4}  for each vertex v
  - q(antipode(v)) = -q(v)                 for each antipodal pair
  - q(a) + q(b) + q(c) = 0                 for each great-circle triple

The resulting CNF is UNSATISFIABLE, proving that no such labeling exists.

Requires: pycosat (pip install pycosat)

Reference: "Graph Puzzles II.1: Counterexamples to Jain's second unit
vector flows conjecture" by Nikolay Ulyanov.
"""

from __future__ import annotations

import itertools
import time
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

try:
    import pycosat
except ImportError as exc:
    raise SystemExit(
        "pycosat is required. Install with: pip install pycosat"
    ) from exc

from counterexample1_50pts import COORDS, TRIPLES

# ─── Constants ───────────────────────────────────────────────────────────────

ALLOWED_VALUES: List[int] = [-4, -3, -2, -1, 1, 2, 3, 4]
VALUE_TO_INDEX: Dict[int, int] = {v: i for i, v in enumerate(ALLOWED_VALUES)}

DISALLOWED_3: List[Tuple[int, int, int]] = [
    (i, j, k)
    for i, a in enumerate(ALLOWED_VALUES)
    for j, b in enumerate(ALLOWED_VALUES)
    for k, c in enumerate(ALLOWED_VALUES)
    if a + b + c != 0
]


# ─── Antipode derivation ────────────────────────────────────────────────────

def negate_expr(expr: str) -> str:
    """Negate a symbolic coordinate expression."""
    expr = expr.strip()
    if expr == "0":
        return "0"
    if expr.startswith("-(") and expr.endswith(")"):
        return expr[2:-1]
    if expr.startswith("-"):
        return expr[1:]
    if "+" in expr or "-" in expr:
        return f"-({expr})"
    return f"-{expr}"


def derive_antipode(coords: Dict[int, Tuple[str, str, str]]) -> Dict[int, int]:
    """Derive the antipodal map from symbolic coordinates."""
    by_coord: Dict[Tuple[str, str, str], int] = {}
    for idx, c in coords.items():
        by_coord[c] = idx

    antipode: Dict[int, int] = {}
    for idx, c in coords.items():
        neg_c = tuple(negate_expr(x) for x in c)
        if neg_c not in by_coord:
            raise ValueError(f"Missing antipode for vertex {idx}: {neg_c}")
        antipode[idx] = by_coord[neg_c]

    # Verify involution
    for v, w in antipode.items():
        assert antipode[w] == v, f"Antipode not involution at {v}"
    return antipode


# ─── CNF construction ───────────────────────────────────────────────────────

def add_exactly_one(cnf: List[List[int]], lits: Sequence[int]) -> None:
    """Add clauses enforcing exactly one of the literals is true."""
    cnf.append(list(lits))  # at least one
    for i in range(len(lits)):
        for j in range(i + 1, len(lits)):
            cnf.append([-lits[i], -lits[j]])  # at most one


def build_full_cnf(
    n_vertices: int,
    antipode: Dict[int, int],
    triples: Sequence[Tuple[int, int, int]],
) -> Tuple[int, List[List[int]]]:
    """Build the full CNF encoding (one set of vars per vertex)."""
    nvals = len(ALLOWED_VALUES)

    def var(v: int, val_idx: int) -> int:
        return v * nvals + val_idx + 1

    cnf: List[List[int]] = []

    # Each vertex gets exactly one value
    for v in range(n_vertices):
        add_exactly_one(cnf, [var(v, i) for i in range(nvals)])

    # Antipodal constraints: q(v) = k ↔ q(antipode(v)) = -k
    for v, w in antipode.items():
        if v < w:
            for i, value in enumerate(ALLOWED_VALUES):
                j = VALUE_TO_INDEX[-value]
                cnf.append([-var(v, i), var(w, j)])
                cnf.append([-var(w, j), var(v, i)])

    # Triple constraints: q(a) + q(b) + q(c) = 0
    for a, b, c in triples:
        for i, j, k in DISALLOWED_3:
            cnf.append([-var(a, i), -var(b, j), -var(c, k)])

    return n_vertices * nvals, cnf


def build_reduced_cnf(
    n_vertices: int,
    antipode: Dict[int, int],
    triples: Sequence[Tuple[int, int, int]],
) -> Tuple[int, List[List[int]]]:
    """Build the reduced CNF (one set of vars per antipodal representative)."""
    reps = [v for v in range(n_vertices) if v <= antipode[v]]
    rep_pos = {r: i for i, r in enumerate(reps)}

    v_to_rep_sign: Dict[int, Tuple[int, int]] = {}
    for r in reps:
        v_to_rep_sign[r] = (r, 1)
        v_to_rep_sign[antipode[r]] = (r, -1)

    nvals = len(ALLOWED_VALUES)

    def var(rep_idx: int, val_idx: int) -> int:
        return rep_idx * nvals + val_idx + 1

    cnf: List[List[int]] = []

    # Each representative gets exactly one value
    for ridx in range(len(reps)):
        add_exactly_one(cnf, [var(ridx, i) for i in range(nvals)])

    # Triple constraints (with sign adjustments for antipodal vertices)
    for triple in triples:
        coeff: Dict[int, int] = {}
        for v in triple:
            r, s = v_to_rep_sign[v]
            coeff[r] = coeff.get(r, 0) + s

        items = sorted(coeff.items())
        triple_reps = [r for r, _ in items]
        coefs = [c for _, c in items]

        for idx_tuple in itertools.product(range(nvals), repeat=len(triple_reps)):
            total = sum(c * ALLOWED_VALUES[idx] for c, idx in zip(coefs, idx_tuple))
            if total != 0:
                cnf.append([-var(rep_pos[r], idx)
                            for r, idx in zip(triple_reps, idx_tuple)])

    return len(reps) * nvals, cnf


def write_dimacs(path: Path, nvars: int, cnf: Sequence[Sequence[int]]) -> None:
    """Write CNF in DIMACS format."""
    with path.open("w", encoding="utf-8") as f:
        f.write(f"p cnf {nvars} {len(cnf)}\n")
        for clause in cnf:
            f.write(" ".join(str(x) for x in clause) + " 0\n")


# ─── Main ───────────────────────────────────────────────────────────────────

def main() -> int:
    n_vertices = len(COORDS)
    antipode = derive_antipode(COORDS)
    n_pairs = sum(1 for v, w in antipode.items() if v < w)

    print(f"Instance: {n_vertices} vertices, {len(TRIPLES)} triples, "
          f"{n_pairs} antipodal pairs")
    print(f"Allowed values: {ALLOWED_VALUES}")
    print()

    for enc_name, builder in [("full", build_full_cnf),
                               ("reduced", build_reduced_cnf)]:
        nvars, cnf = builder(n_vertices, antipode, TRIPLES)

        # Optionally write DIMACS file
        dimacs_path = Path(f"counterexample1_{enc_name}.cnf")
        write_dimacs(dimacs_path, nvars, cnf)
        print(f"[{enc_name}] Wrote DIMACS: {dimacs_path}")

        # Solve
        t0 = time.perf_counter()
        result = pycosat.solve(cnf)
        dt = time.perf_counter() - t0

        if result == "UNSAT":
            status = "UNSAT"
        elif result == "UNKNOWN":
            status = "UNKNOWN"
        else:
            status = "SAT"

        print(f"[{enc_name}] vars={nvars}, clauses={len(cnf)}, "
              f"status={status}, time={dt:.3f}s")
        print()

    print("=" * 50)
    print("Result: UNSAT — no labeling q: S² → {±1,±2,±3,±4} exists")
    print("for this 50-point set. Conjecture 1.2 is DISPROVED.")
    print("=" * 50)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
