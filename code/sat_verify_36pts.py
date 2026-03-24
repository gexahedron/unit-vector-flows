#!/usr/bin/env python3
"""
SAT-based verification that the 36-point counterexample has no nz5-flow labeling
but does admit a nz6-flow labeling.

Encodes the labeling problem as a Boolean satisfiability (SAT) instance using
the *reduced* encoding (one variable block per antipodal representative):

  nz5-flow: q(v) ∈ {-4,-3,-2,-1,1,2,3,4}  →  18 reps × 8 values = 144 variables
  nz6-flow: q(v) ∈ {-5,-4,-3,-2,-1,1,2,3,4,5}  →  18 reps × 10 values = 180 variables

Constraints:
  - Exactly one value per representative (exactly-one clauses)
  - q(a) + q(b) + q(c) = 0 for each great-circle triple (forbidden-combo clauses)
    (antipodal constraint is implicit: representative carries sign ±1)

Results:
  nz5-flow: UNSAT  (144 variables, 6710 clauses)  — no labeling exists
  nz6-flow: SAT    (180 variables, 13048 clauses) — labeling exists

Requires: pycosat (pip install pycosat)

Reference: "Graph Puzzles II.1: Counterexamples to Jain's second unit
vector flows conjecture" by Nikolay Ulyanov, 2026.
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

from counterexample2_36pts import POINTS, TRIPLES, ANTIPODE


# ─── CNF construction ────────────────────────────────────────────────────────

def add_exactly_one(cnf: List[List[int]], lits: Sequence[int]) -> None:
    """Add clauses enforcing exactly one of the literals is true."""
    cnf.append(list(lits))
    for i in range(len(lits)):
        for j in range(i + 1, len(lits)):
            cnf.append([-lits[i], -lits[j]])


def build_reduced_cnf(
    n_vertices: int,
    antipode: Dict[int, int],
    triples: Sequence[Tuple[int, int, int]],
    allowed_values: List[int],
) -> Tuple[int, List[List[int]]]:
    """
    Build the reduced CNF encoding.

    One block of variables per antipodal representative (v <= antipode[v]).
    The antipodal constraint q(antipode(v)) = -q(v) is handled implicitly:
    each triple vertex is mapped to its representative with a sign, and the
    forbidden-combo check uses the signed value.

    Returns (nvars, cnf) where nvars = n_reps * len(allowed_values).
    """
    nvals = len(allowed_values)
    value_to_idx: Dict[int, int] = {v: i for i, v in enumerate(allowed_values)}

    # Representatives: vertices v with v <= antipode[v]
    reps: List[int] = [v for v in range(n_vertices) if v <= antipode[v]]
    rep_pos: Dict[int, int] = {r: i for i, r in enumerate(reps)}

    # Map every vertex to (representative, sign)
    v_to_rep_sign: Dict[int, Tuple[int, int]] = {}
    for r in reps:
        v_to_rep_sign[r] = (r, +1)
        v_to_rep_sign[antipode[r]] = (r, -1)

    def var(rep_idx: int, val_idx: int) -> int:
        """1-based SAT variable for representative rep_idx having value allowed_values[val_idx]."""
        return rep_idx * nvals + val_idx + 1

    cnf: List[List[int]] = []

    # Exactly-one constraint for each representative
    for ridx in range(len(reps)):
        add_exactly_one(cnf, [var(ridx, i) for i in range(nvals)])

    # Triple constraints: for each triple (a, b, c), forbid all value combos
    # where the signed sum is non-zero.
    for triple in triples:
        # Compute (representative, coefficient) for each vertex in the triple.
        # A vertex may appear with sign +1 or -1 depending on whether it is
        # the representative or its antipode.
        coeff: Dict[int, int] = {}
        for v in triple:
            r, s = v_to_rep_sign[v]
            coeff[r] = coeff.get(r, 0) + s

        items = sorted(coeff.items())          # deterministic order
        triple_reps = [r for r, _ in items]
        coefs       = [c for _, c in items]

        for idx_tuple in itertools.product(range(nvals), repeat=len(triple_reps)):
            total = sum(c * allowed_values[idx]
                        for c, idx in zip(coefs, idx_tuple))
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


def decode_solution(
    solution: List[int],
    n_vertices: int,
    antipode: Dict[int, int],
    allowed_values: List[int],
) -> Dict[int, int]:
    """Decode a pycosat solution into a labeling q: {0,...,n-1} → allowed_values."""
    nvals = len(allowed_values)
    reps = [v for v in range(n_vertices) if v <= antipode[v]]

    labeling: Dict[int, int] = {}
    for ridx, r in enumerate(reps):
        for val_idx, val in enumerate(allowed_values):
            lit = ridx * nvals + val_idx + 1
            if lit in solution:
                labeling[r] = val
                labeling[antipode[r]] = -val
                break

    return labeling


# ─── Main ────────────────────────────────────────────────────────────────────

def main() -> int:
    n_vertices = len(POINTS)
    antipode   = ANTIPODE
    n_pairs    = sum(1 for v, w in antipode.items() if v < w)

    print("=" * 60)
    print("SAT verification: 36-point counterexample to S² nz5-flow")
    print("=" * 60)
    print(f"\n  Vertices:       {n_vertices}")
    print(f"  Triples:        {len(TRIPLES)}")
    print(f"  Antipodal pairs:{n_pairs}")
    print(f"  Encoding:       reduced (one block per antipodal rep)")
    print()

    configs = [
        ("nz5", list(range(-4, 0)) + list(range(1, 5)),  "UNSAT"),
        ("nz6", list(range(-5, 0)) + list(range(1, 6)),  "SAT"),
    ]

    all_ok = True

    for label, allowed, expected in configs:
        nvars, cnf = build_reduced_cnf(n_vertices, antipode, TRIPLES, allowed)

        dimacs_path = Path(f"counterexample2_{label}.cnf")
        write_dimacs(dimacs_path, nvars, cnf)

        t0     = time.perf_counter()
        result = pycosat.solve(cnf)
        dt     = time.perf_counter() - t0

        if result == "UNSAT":
            status = "UNSAT"
        elif result == "UNKNOWN":
            status = "UNKNOWN"
        else:
            status = "SAT"

        ok = (status == expected)
        mark = "✓" if ok else "✗"
        all_ok = all_ok and ok

        print(f"  [{label}] vars={nvars}, clauses={len(cnf)}, "
              f"result={status}, time={dt:.3f}s  {mark}")

        if status == "SAT" and isinstance(result, list):
            labeling = decode_solution(result, n_vertices, antipode, allowed)
            # Spot-check: verify all triples sum to zero
            triple_ok = all(
                labeling[a] + labeling[b] + labeling[c] == 0
                for a, b, c in TRIPLES
            )
            antipodal_ok = all(
                labeling[v] + labeling[antipode[v]] == 0
                for v in range(n_vertices)
            )
            if triple_ok and antipodal_ok:
                print(f"         (labeling verified: all triples sum to 0, "
                      f"antipodal constraint satisfied)")
            else:
                print(f"         WARNING: labeling verification FAILED")
                all_ok = False

    print()
    print("=" * 60)
    if all_ok:
        print("Result: nz5-flow UNSAT, nz6-flow SAT.")
        print("The 36-point set is a counterexample to Jain's S² nz5-flow")
        print("conjecture (Conjecture 1.2 is DISPROVED).")
    else:
        print("Result: UNEXPECTED — check output above.")
    print("=" * 60)

    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
