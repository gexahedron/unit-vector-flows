# Counterexamples to Jain's Second Unit Vector Flows Conjecture

Code and data accompanying the paper:

> **Graph Puzzles II.1: Counterexamples to Jain's second unit vector flows conjecture**  
> Nikolay Ulyanov, 2026  
> [`gp2-1.pdf`](gp2-1.pdf)

---

## Background

K. Jain posed two conjectures on [Open Problem Garden](http://www.openproblemgarden.org/op/unit_vector_flows) website about unit vector flows on graphs. The **second conjecture** ($S^2$ nz5-flow conjecture) states:

> There exists a map $q : S^2 \to \{-4,-3,-2,-1,1,2,3,4\}$ such that antipodal points sum to 0, and any three points equidistant on a great circle also sum to 0.

Together with Jain's first conjecture (every bridgeless graph has a unit vector flow), this would imply **Tutte's nowhere-zero 5-flow conjecture**.

**This repository disproves the second conjecture** by exhibiting two finite subsets of $S^2$ for which no such labeling exists — each requires values $\pm 5$ (a nowhere-zero 6-flow).

---

## Counterexample 1: 50-point expansion of the icosidodecahedron

**50 points, 40 great-circle triples, 25 antipodal pairs.**

Starting from the 30 vertices of the icosidodecahedron (which correspond to the Petersen graph and admit a nz5-flow), we expand by intersecting pairs of small circles at height $h = 1/2$ above each vertex. The resulting 50-point set with 40 great-circle triples requires a nz6-flow.

See it in (motion)[https://www.youtube.com/watch?v=W0n6gg3YEm0]!

Coordinates lie in $\mathbb{Q}(\varphi, x)$ where $\varphi = \frac{1+\sqrt{5}}{2}$ is the golden ratio and $x = 2/5^{1/4}$.

On the antipodal quotient, the 10 "old" triple-pairs form the **Petersen graph** and the 10 "new" triple-pairs form a **10-vertex Möbius ladder**. An explicit nz6-flow labeling is shown in [`figures/counterexample1_petersen_moebius.pdf`](figures/counterexample1_petersen_moebius.pdf).

**SAT verification (reduced encoding, 25 reps × 8 values):**
- nz5-flow (UNSAT): **200 variables, 19765 clauses**

**Files:**
- [`code/counterexample1_50pts.py`](code/counterexample1_50pts.py) — point coordinates and triples
- [`code/sat_verify_50pts.py`](code/sat_verify_50pts.py) — SAT encoding and verification
- [`code/draw_counterexample1_structure.py`](code/draw_counterexample1_structure.py) — generates the Petersen+Möbius DOT figure
- [`dot/counterexample1_petersen_moebius.dot`](dot/counterexample1_petersen_moebius.dot) — DOT source for the figure

**Run:**
```bash
cd code
python3 counterexample1_50pts.py        # verify geometry
python3 sat_verify_50pts.py             # verify UNSAT (requires pycosat)
python3 draw_counterexample1_structure.py
neato -Tpdf ../dot/counterexample1_petersen_moebius.dot \
      -o ../dot/counterexample1_petersen_moebius.pdf
```

---

## Counterexample 2: 36-point construction via square root arithmetic

**36 points (18 antipodal pairs), 13 great-circle triples.**

We fix parameters $v_1 = 1$, $v_2 = 3$, $w = 2$ and generate candidate coordinate values:

$$c_{\text{rat}} = \frac{|w_1\sqrt{v_1} + w_2\sqrt{v_2}|}{2}, \qquad c_{\text{sqrt}} = \sqrt{\frac{|w_1\sqrt{v_1} + w_2\sqrt{v_2}|}{2}}$$

for integers $w_1 \in \{0,1,2\}$, $w_2 \in \{-2,\ldots,2\}$, keeping only $c \le 1$. We find all points on $S^2$ whose coordinates come from this set, identify great-circle triples, and prune iteratively.

The 7 distinct absolute coordinate values are:
$$0,\quad \frac{2-\sqrt{3}}{2},\quad \frac{\sqrt{3}-1}{2},\quad \frac{1}{2},\quad \sqrt{\sqrt{3}-1},\quad \frac{\sqrt{3}}{2},\quad 1$$

All coordinates lie in $\mathbb{Q}(\sqrt{3},\sqrt{\sqrt{3}-1})$, a degree-4 extension of $\mathbb{Q}$, with denominator 2.

**Construction pipeline:**
1. Generate 126-point connected component with 108 triples
2. Prune iteratively (remove triples where ≥2 of 3 vertices appear in only one triple)
3. Converge to **36 points, 13 triples**

**SAT verification (reduced encoding, 18 reps):**
- nz5-flow (UNSAT): **144 variables, 6710 clauses**
- nz6-flow (SAT): **180 variables, 13048 clauses**

**Files:**
- [`code/counterexample2_36pts.py`](code/counterexample2_36pts.py) — hardcoded 36-point data (exact symbolic coordinates, triples, antipodal pairs)
- [`code/sat_verify_36pts.py`](code/sat_verify_36pts.py) — SAT encoding and verification (nz5 UNSAT + nz6 SAT)
- [`code/prune_w1_component.py`](code/prune_w1_component.py) — full construction and pruning pipeline (reproduces the 36-point set from scratch)
- [`code/counterexample2_sqrt_complex.py`](code/counterexample2_sqrt_complex.py) — general square root complex construction (larger parameter $w=7$, 402-point set)

**Run:**
```bash
cd code
python3 counterexample2_36pts.py        # verify geometry (sphere, antipodal, triples)
python3 sat_verify_36pts.py             # verify nz5 UNSAT + nz6 SAT (requires pycosat)
python3 prune_w1_component.py           # reproduce 36-point set from scratch (requires pycosat)
```

---

## Icosidodecahedron and nz5-flow (non-counterexample)

As a warm-up, the paper also shows that the 30-point icosidodecahedron **does** admit a nz5-flow. Its 20 great-circle triples (10 antipodal pairs) correspond to the Petersen graph, which has a nz5-flow. The explicit labeling is shown in [`figures/petersen_uvf.pdf`](figures/petersen_uvf.pdf).

---

## Repository structure

```
unit-vector-flows/
├── gp2-1.pdf                  # Compiled paper
├── code/
│   ├── counterexample1_50pts.py          # CE1: 50-point data
│   ├── sat_verify_50pts.py               # CE1: SAT verification
│   ├── draw_counterexample1_structure.py # CE1: figure generation
│   ├── draw_sphere_counterexample1.py    # CE1: sphere visualization
│   ├── counterexample2_36pts.py          # CE2: 36-point data (hardcoded)
│   ├── sat_verify_36pts.py               # CE2: SAT verification (36-pt)
│   ├── counterexample2_sqrt_complex.py   # CE2: general sqrt construction
│   ├── prune_w1_component.py             # CE2: 36-point pruning pipeline
│   ├── generate_scaled_dot.py            # DOT generation utility
│   └── requirements.txt
├── dot/
│   ├── counterexample1_petersen_moebius.dot   # CE1 figure source
│   ├── counterexample2_36pts_active.dot       # CE2 hypergraph (active vertices)
│   ├── counterexample2_36pts_clean.dot        # CE2 hypergraph (all vertices)
│   └── *.pdf                                  # Rendered figures
├── figures/
│   ├── counterexample1_petersen_moebius.pdf   # CE1 main figure
│   └── petersen_uvf.pdf                       # Icosidodecahedron nz5-flow
└── lean/                      # Lean 4 formal verification (CE1)
    ├── lakefile.toml           # Lake build config (Mathlib v4.26.0-rc2)
    ├── lean-toolchain          # Lean version pin (leanprover/lean4:v4.26.0-rc2)
    ├── S2FlowsInLean.lean      # Root import file
    ├── S2FlowsInLean/
    │   └── SecondConjecture50BVDecide.lean  # bv_decide proof (auto-generated)
    └── scripts/
        └── generate_second50_bvdecide.py    # Generates the BVDecide Lean file
```

---

## Requirements

- Python 3.8+
- NumPy
- pycosat (SAT solver): `pip install pycosat`
- Graphviz (`neato`) for rendering DOT files

```bash
pip install -r code/requirements.txt
```

---

## Formal verification (Lean 4)

The nz5-flow impossibility for **Counterexample 1** (50-point set) has been formally verified in **Lean 4** using the `bv_decide` tactic from Mathlib.

The proof encodes each of the 50 antipodal representatives as a `BitVec 8` variable, asserts the domain constraint ($q(v) \in \{-4,\ldots,-1,1,\ldots,4\}$), the antipodal constraint ($q(\bar v) = -q(v)$), and the triple-sum constraint ($q(a)+q(b)+q(c)=0$) as a single Boolean formula, then calls `bv_decide` to certify UNSAT.

**Files:**
- [`lean/S2FlowsInLean/SecondConjecture50BVDecide.lean`](lean/S2FlowsInLean/SecondConjecture50BVDecide.lean) — the Lean 4 proof (auto-generated)
- [`lean/scripts/generate_second50_bvdecide.py`](lean/scripts/generate_second50_bvdecide.py) — Python script that generates the Lean file from the 50-point data
- [`lean/S2FlowsInLean.lean`](lean/S2FlowsInLean.lean) — root import file
- [`lean/lakefile.toml`](lean/lakefile.toml) — Lake build config (Mathlib `v4.26.0-rc2`)
- [`lean/lean-toolchain`](lean/lean-toolchain) — Lean version pin (`leanprover/lean4:v4.26.0-rc2`)

**To regenerate and check the proof:**
```bash
cd lean
# Regenerate the Lean file from the 50-point data:
python3 scripts/generate_second50_bvdecide.py
# Build with Lake (downloads Mathlib on first run, may take a while):
lake build
```

---

## License

MIT. See [`LICENSE`](LICENSE).
