#!/usr/bin/env python3
"""
Derive the graph structure behind Counterexample 1 and render a DOT figure.

The 50-point construction has 25 antipodal point-pairs and 20 antipodal
triple-pairs. The 10 triple-pairs supported only on icosidodecahedron points
form a Petersen graph. The 10 mixed triple-pairs form a Möbius ladder on
10 vertices.

The script also finds one integer nz6-flow labeling q : S^2 -> {±1,...,±5} by
backtracking and annotates each edge-orbit with its representative point and q.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

from counterexample1_50pts import COORDS, TRIPLES

ALLOWED_FLOW_VALUES = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]
PETERSEN_UVF_POSITIONS = {
    0: (-10.00000, 2.00000),
    1: (-11.73210, -1.00000),
    2: (-8.26795, -1.00000),
    3: (-10.00000, 0.00000),
    4: (-9.31596, -1.87939),
    5: (-10.68400, -1.87939),
    6: (-8.71443, 1.53209),
    7: (-8.03038, 0.347296),
    8: (-11.28560, 1.53209),
    9: (-11.96960, 0.347296),
}
MIXED_POSITIONS = {
    0: (-2.80000, 1.99872),
    1: (-1.62428, 1.61661),
    2: (-0.89763, 0.61725),
    3: (-0.89763, -0.61725),
    4: (-1.62428, -1.61661),
    5: (-2.80000, -1.99872),
    6: (-3.97572, -1.61661),
    7: (-4.70237, -0.61725),
    8: (-4.70237, 0.61725),
    9: (-3.97572, 1.61661),
}


@dataclass(frozen=True)
class TriplePair:
    pair_id: int
    triple: Tuple[int, int, int]
    antipodal_triple: Tuple[int, int, int]
    pair_reps: Tuple[int, int, int]
    new_pair_count: int


@dataclass(frozen=True)
class EdgeOrbit:
    rep: int
    antipode: int
    coord: Tuple[str, str, str]
    kind: str
    old_endpoints: Tuple[int, ...]
    mixed_endpoints: Tuple[int, ...]


def negate_expr(expr: str) -> str:
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


def build_antipode_map() -> Dict[int, int]:
    by_coord = {tuple(coord): idx for idx, coord in COORDS.items()}
    antipode: Dict[int, int] = {}
    for idx, coord in COORDS.items():
        neg_coord = tuple(negate_expr(part) for part in coord)
        antipode[idx] = by_coord[neg_coord]
    return antipode


def is_new_point(idx: int) -> bool:
    return any(token in "".join(COORDS[idx]) for token in ("x", "y"))


def pair_rep(idx: int, antipode: Dict[int, int]) -> int:
    return min(idx, antipode[idx])


def triple_pair_key(
    triple: Sequence[int], antipode: Dict[int, int]
) -> Tuple[Tuple[int, int, int], Tuple[int, int, int]]:
    triple_sorted = tuple(sorted(triple))
    anti_sorted = tuple(sorted(antipode[v] for v in triple))
    return min((triple_sorted, anti_sorted), (anti_sorted, triple_sorted))


def build_triple_pairs(antipode: Dict[int, int]) -> List[TriplePair]:
    pair_id_by_key: Dict[
        Tuple[Tuple[int, int, int], Tuple[int, int, int]], int
    ] = {}
    triple_pairs: List[TriplePair] = []
    for triple in TRIPLES:
        key = triple_pair_key(triple, antipode)
        if key in pair_id_by_key:
            continue
        triple_sorted, anti_sorted = key
        reps = tuple(pair_rep(v, antipode) for v in triple_sorted)
        triple_pairs.append(
            TriplePair(
                pair_id=len(triple_pairs),
                triple=triple_sorted,
                antipodal_triple=anti_sorted,
                pair_reps=reps,
                new_pair_count=sum(1 for v in reps if is_new_point(v)),
            )
        )
        pair_id_by_key[key] = len(triple_pairs) - 1
    return triple_pairs


def build_edge_orbits(
    antipode: Dict[int, int], triple_pairs: Sequence[TriplePair]
) -> List[EdgeOrbit]:
    incident_triple_pairs: Dict[int, List[int]] = defaultdict(list)
    triple_pair_by_key = {
        triple_pair_key(tp.triple, antipode): tp.pair_id for tp in triple_pairs
    }

    for triple in TRIPLES:
        key = triple_pair_key(triple, antipode)
        tp_id = triple_pair_by_key[key]
        for vertex in triple:
            incident_triple_pairs[pair_rep(vertex, antipode)].append(tp_id)

    edge_orbits: List[EdgeOrbit] = []
    for rep in sorted(v for v in COORDS if v < antipode[v]):
        tp_ids = tuple(sorted(set(incident_triple_pairs[rep])))
        old_endpoints = tuple(tp for tp in tp_ids if triple_pairs[tp].new_pair_count == 0)
        mixed_endpoints = tuple(
            tp for tp in tp_ids if triple_pairs[tp].new_pair_count == 2
        )
        if is_new_point(rep):
            kind = "new_only"
        elif len(old_endpoints) == 2 and len(mixed_endpoints) == 2:
            kind = "shared"
        else:
            kind = "old_only"
        edge_orbits.append(
            EdgeOrbit(
                rep=rep,
                antipode=antipode[rep],
                coord=COORDS[rep],
                kind=kind,
                old_endpoints=old_endpoints,
                mixed_endpoints=mixed_endpoints,
            )
        )
    return edge_orbits


def find_nz6_flow(antipode: Dict[int, int]) -> Dict[int, int]:
    reps = sorted(v for v in COORDS if v < antipode[v])
    constraints: List[Tuple[Tuple[int, int], Tuple[int, int], Tuple[int, int]]] = []
    seen = set()
    for triple in TRIPLES:
        signed = []
        for vertex in triple:
            if vertex < antipode[vertex]:
                signed.append((vertex, 1))
            else:
                signed.append((antipode[vertex], -1))
        normalized = tuple(sorted(signed))
        flipped = tuple(sorted((rep, -sign) for rep, sign in signed))
        key = min(normalized, flipped)
        if key in seen:
            continue
        seen.add(key)
        constraints.append(tuple(signed))

    supported_tuples = [
        [
            values
            for values in _all_zero_sum_triples(ALLOWED_FLOW_VALUES)
        ]
        for _ in constraints
    ]

    domains = {rep: set(ALLOWED_FLOW_VALUES) for rep in reps}
    assignments: Dict[int, int] = {}
    constraints_by_rep: Dict[int, List[int]] = defaultdict(list)
    for idx, constraint in enumerate(constraints):
        for rep, _ in constraint:
            constraints_by_rep[rep].append(idx)

    def revise_constraint(idx: int) -> bool | None:
        constraint = constraints[idx]
        feasible = []
        for values in supported_tuples[idx]:
            ok = True
            for (rep, sign), value in zip(constraint, values):
                signed_value = sign * value
                if rep in assignments:
                    if assignments[rep] != signed_value:
                        ok = False
                        break
                elif signed_value not in domains[rep]:
                    ok = False
                    break
            if ok:
                feasible.append(values)
        if not feasible:
            return None
        changed = False
        for pos, (rep, sign) in enumerate(constraint):
            if rep in assignments:
                continue
            supported = {sign * values[pos] for values in feasible}
            new_domain = domains[rep] & supported
            if not new_domain:
                return None
            if new_domain != domains[rep]:
                domains[rep] = new_domain
                changed = True
        return changed

    def propagate(queue: List[int]) -> bool:
        while queue:
            idx = queue.pop()
            changed = revise_constraint(idx)
            if changed is None:
                return False
            if changed:
                for rep, _ in constraints[idx]:
                    if rep in assignments:
                        continue
                    queue.extend(j for j in constraints_by_rep[rep] if j != idx)
        return True

    if not propagate(list(range(len(constraints)))):
        raise RuntimeError("Initial nz6 propagation failed")

    def dfs() -> bool:
        forced: List[int] = []
        progress = True
        while progress:
            progress = False
            for rep in reps:
                if rep in assignments or len(domains[rep]) != 1:
                    continue
                assignments[rep] = next(iter(domains[rep]))
                forced.append(rep)
                progress = True
                if not propagate(constraints_by_rep[rep][:]):
                    for forced_rep in reversed(forced):
                        assignments.pop(forced_rep, None)
                    return False

        if len(assignments) == len(reps):
            return True

        target = min((rep for rep in reps if rep not in assignments), key=lambda rep: len(domains[rep]))
        snapshot = {rep: set(values) for rep, values in domains.items()}
        for value in sorted(domains[target], key=lambda value: (abs(value), value)):
            assignments[target] = value
            domains[target] = {value}
            if propagate(constraints_by_rep[target][:]) and dfs():
                return True
            assignments.pop(target, None)
            domains.clear()
            domains.update({rep: set(values) for rep, values in snapshot.items()})

        for forced_rep in reversed(forced):
            assignments.pop(forced_rep, None)
        domains.clear()
        domains.update(snapshot)
        return False

    if not dfs():
        raise RuntimeError("No nz6 flow found")
    return assignments


def _all_zero_sum_triples(values: Sequence[int]) -> Iterable[Tuple[int, int, int]]:
    for a in values:
        for b in values:
            for c in values:
                if a + b + c == 0:
                    yield (a, b, c)


def pretty_expr(expr: str) -> str:
    return expr


def coord_label(coord: Tuple[str, str, str]) -> str:
    x, y, z = (pretty_expr(part) for part in coord)
    return f"{x}, {y}, {z}"


def q_label(edge: EdgeOrbit, flow: Dict[int, int]) -> str:
    return f"q {abs(flow[edge.rep])}"


def base_incidence_sign(
    tp: TriplePair, edge: EdgeOrbit, antipode: Dict[int, int]
) -> int:
    if edge.rep in tp.triple:
        return 1
    if antipode[edge.rep] in tp.triple:
        return -1
    raise ValueError(f"Triple-pair {tp.pair_id} is not incident to edge orbit {edge.rep}")


def solve_vertex_switches(
    triple_pairs: Sequence[TriplePair],
    edge_orbits: Sequence[EdgeOrbit],
    antipode: Dict[int, int],
    endpoint_attr: str,
) -> Dict[int, int]:
    by_id = {tp.pair_id: tp for tp in triple_pairs}
    adjacency: Dict[int, List[Tuple[int, int]]] = defaultdict(list)

    for edge in edge_orbits:
        endpoints = getattr(edge, endpoint_attr)
        if len(endpoints) != 2:
            continue
        u, v = endpoints
        sign_u = base_incidence_sign(by_id[u], edge, antipode)
        sign_v = base_incidence_sign(by_id[v], edge, antipode)
        relation = -sign_u * sign_v
        adjacency[u].append((v, relation))
        adjacency[v].append((u, relation))

    switches: Dict[int, int] = {}
    for start in adjacency:
        if start in switches:
            continue
        switches[start] = 1
        stack = [start]
        while stack:
            u = stack.pop()
            for v, relation in adjacency[u]:
                wanted = switches[u] * relation
                if v in switches:
                    if switches[v] != wanted:
                        raise RuntimeError(
                            f"Inconsistent switching constraints on {endpoint_attr}"
                        )
                    continue
                switches[v] = wanted
                stack.append(v)
    return switches


def render_dot(
    triple_pairs: Sequence[TriplePair],
    edge_orbits: Sequence[EdgeOrbit],
    flow: Dict[int, int],
    antipode: Dict[int, int],
) -> str:
    by_id = {tp.pair_id: tp for tp in triple_pairs}
    edge_by_old = {edge.old_endpoints: edge for edge in edge_orbits if len(edge.old_endpoints) == 2}
    edge_by_mixed = {
        edge.mixed_endpoints: edge for edge in edge_orbits if len(edge.mixed_endpoints) == 2
    }
    old_switches = solve_vertex_switches(
        triple_pairs, edge_orbits, antipode, "old_endpoints"
    )
    mixed_switches = solve_vertex_switches(
        triple_pairs, edge_orbits, antipode, "mixed_endpoints"
    )

    # This isomorphism identifies the all-old quotient with the numbering used
    # in dot/petersen_uvf.dot, so the left panel can reuse that exact layout.
    p_tp_by_label = {
        0: 2,
        1: 5,
        2: 11,
        3: 17,
        4: 13,
        5: 12,
        6: 19,
        7: 16,
        8: 8,
        9: 18,
    }
    p_labels = list(range(10))
    p_names = {label: f"p{label}" for label in p_labels}
    p_positions = {
        p_names[label]: PETERSEN_UVF_POSITIONS[label] for label in p_labels
    }
    p_uvf_edges = [
        (8, 9, 5.0, "tail"),
        (5, 4, 3.0, "tail"),
        (7, 6, 3.0, "tail"),
        (9, 1, 3.0, "tail"),
        (0, 4, 3.0, "tail"),
        (7, 3, 2.5, "tail"),
        (1, 5, 3.0, "tail"),
        (2, 7, 3.0, "tail"),
        (0, 8, 3.0, "tail"),
        (4, 2, 3.0, "tail"),
        (1, 6, 3.0, "head"),
        (3, 8, 5.0, "head"),
        (6, 0, 3.0, "tail"),
        (2, 9, 3.5, "head"),
        (5, 3, 2.0, "tail"),
    ]

    m_cycle = [0, 9, 4, 1, 15, 3, 14, 6, 7, 10]
    m_names = {tp_id: f"m{idx}" for idx, tp_id in enumerate(m_cycle)}
    m_labels = {tp_id: f"{idx}" for idx, tp_id in enumerate(m_cycle)}
    m_positions = {m_names[tp_id]: MIXED_POSITIONS[idx] for idx, tp_id in enumerate(m_cycle)}
    m_cycle_edges = list(zip(m_cycle, m_cycle[1:] + m_cycle[:1]))
    m_match = [(0, 3), (9, 14), (4, 6), (1, 7), (15, 10)]

    lines: List[str] = []
    lines.append("digraph {")
    lines.append('  graph [fontname="georgia" pad="0.22,0.06" bgcolor=white outputorder="edgesfirst" splines=spline dpi=300 margin=0]')
    lines.append('  node [fontname="georgia" style="filled" shape=circle fixedsize=shape width=0.45 height=0.45 fontsize="14"]')
    lines.append('  edge [fontname="georgia" color="gray29" minlen=3 labeldistance=3 penwidth=1.8 fontsize="9" arrowhead=normal arrowsize=0.7]')
    lines.append("")
    lines.append("  // Left side: the 10 all-old triple pairs form the Petersen graph.")
    lines.append("  // Right side: the 10 mixed triple pairs form a Mobius ladder M10.")
    lines.append("  // The 5 shared old point-pairs appear in both quotients.")
    lines.append("")
    lines.append("  // Point-pair orbit labels u_k used on the edges")
    for edge in edge_orbits:
        x, y, z = (pretty_expr(part) for part in edge.coord)
        lines.append(
            f"  // u{edge.rep}: ({x}, {y}, {z}), q = {flow[edge.rep]:+d}, kind = {edge.kind}"
        )
    lines.append("")
    lines.append("  // Petersen vertex mapping")
    for label in p_labels:
        tp_id = p_tp_by_label[label]
        tp = by_id[tp_id]
        lines.append(
            f"  // {label} = T{tp_id}: {tp.triple} / {tp.antipodal_triple}"
        )
    lines.append("  // Möbius vertex mapping")
    for tp_id in m_cycle:
        tp = by_id[tp_id]
        lines.append(
            f"  // {m_labels[tp_id]} = T{tp_id}: {tp.triple} / {tp.antipodal_triple}"
        )
    lines.append("")

    for label in p_labels:
        node_name = p_names[label]
        x, y = p_positions[node_name]
        lines.append(
            f'  {node_name} [label="{label}" fillcolor="lightskyblue" color="white" pos="{x:.5f},{y:.5f}!", pin=true];'
        )
    for tp_id in m_cycle:
        node_name = m_names[tp_id]
        x, y = m_positions[node_name]
        lines.append(
            f'  {node_name} [label="{m_labels[tp_id]}" fillcolor="hotpink" color="white" pos="{x:.5f},{y:.5f}!", pin=true];'
        )
    lines.append('  t_left [shape=plaintext, style=solid, color=white, fillcolor=white, fontsize=18, label="Petersen graph", pos="-10,3.3!", pin=true];')
    lines.append('  t_right [shape=plaintext, style=solid, color=white, fillcolor=white, fontsize=18, label="Mixed quotient", pos="-2.8,3.3!", pin=true];')
    lines.append("")

    def add_edge(
        node_map: Dict[int, str],
        u: int,
        v: int,
        edge: EdgeOrbit,
        color: str,
        switches: Dict[int, int],
        *,
        incidence_u: int | None = None,
        incidence_v: int | None = None,
        coord_side: str = "tail",
        labeldistance: float = 3.0,
    ) -> None:
        attrs = [f'color="{color}"', f'labeldistance={labeldistance}']
        if incidence_u is None:
            incidence_u = u
        if incidence_v is None:
            incidence_v = v

        sign_u = switches[incidence_u] * base_incidence_sign(by_id[incidence_u], edge, antipode)
        sign_v = switches[incidence_v] * base_incidence_sign(by_id[incidence_v], edge, antipode)
        if sign_u != -sign_v:
            raise RuntimeError(f"Inconsistent edge orientation data for orbit {edge.rep}")

        tail_tp = u
        head_tp = v
        tail_sign = sign_u
        if sign_u * flow[edge.rep] < 0:
            tail_tp = v
            head_tp = u
            tail_sign = sign_v
            coord_side = "head" if coord_side == "tail" else "tail"

        tail_coord = edge.coord if tail_sign > 0 else COORDS[edge.antipode]
        coord = coord_label(tail_coord)
        q = q_label(edge, flow)
        edge_text = f"{coord} / {q}"
        if coord_side == "tail":
            attrs.append(f'taillabel="{edge_text}"')
        else:
            attrs.append(f'headlabel="{edge_text}"')
        return f'  {node_map[tail_tp]} -> {node_map[head_tp]} [{", ".join(attrs)}];'

    lines.append("  // Petersen edges in the same order as dot/petersen_uvf.dot")
    for u_label, v_label, labeldistance, coord_side in p_uvf_edges:
        tp_u = p_tp_by_label[u_label]
        tp_v = p_tp_by_label[v_label]
        edge = edge_by_old[tuple(sorted((tp_u, tp_v)))]
        color = "purple" if edge.kind == "shared" else "steelblue3"
        lines.append(
            add_edge(
                p_names,
                u_label,
                v_label,
                edge,
                color,
                old_switches,
                incidence_u=tp_u,
                incidence_v=tp_v,
                coord_side=coord_side,
                labeldistance=labeldistance,
            )
        )

    lines.append("")
    lines.append("  // Möbius ladder: 10-cycle on the mixed triple pairs (new-only edges)")
    cycle_coord_side = ["tail", "head", "tail", "head", "tail", "head", "tail", "head", "tail", "head"]
    cycle_dist = [2.4, 2.6, 2.8, 2.8, 2.4, 2.4, 2.8, 2.8, 2.6, 2.4]
    for idx, (u, v) in enumerate(m_cycle_edges):
        lines.append(add_edge(
            m_names,
            u,
            v,
            edge_by_mixed[tuple(sorted((u, v)))],
            "orangered3",
            mixed_switches,
            coord_side=cycle_coord_side[idx],
            labeldistance=cycle_dist[idx],
        ))
    lines.append("  // Möbius ladder: matching coming from the 5 shared old point-pairs")
    match_coord_side = {
        (0, 3): "tail",
        (9, 14): "head",
        (4, 6): "tail",
        (1, 7): "head",
        (15, 10): "tail",
    }
    for u, v in m_match:
        lines.append(add_edge(
            m_names,
            u,
            v,
            edge_by_mixed[tuple(sorted((u, v)))],
            "purple",
            mixed_switches,
            coord_side=match_coord_side[(u, v)],
            labeldistance=2.4,
        ))

    lines.append("}")
    return "\n".join(lines) + "\n"


def main() -> int:
    repo_root = Path(__file__).resolve().parent.parent
    dot_dir = repo_root / "dot"

    antipode = build_antipode_map()
    triple_pairs = build_triple_pairs(antipode)
    edge_orbits = build_edge_orbits(antipode, triple_pairs)
    flow = find_nz6_flow(antipode)

    dot_text = render_dot(triple_pairs, edge_orbits, flow, antipode)
    output_path = dot_dir / "counterexample1_petersen_mobius.dot"
    output_path.write_text(dot_text, encoding="utf-8")

    print("Wrote:")
    print(f"  {output_path}")
    print("\nRepresentative nz6 values:")
    for edge in edge_orbits:
        print(f"  {edge.rep:2d}/{edge.antipode:2d}: q = {flow[edge.rep]:+d} ({edge.kind})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
