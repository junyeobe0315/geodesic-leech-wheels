#!/usr/bin/env python3
"""Verify the explicit geodesic-Leech certificates for W_7,...,W_13.

This is the computational side of Theorem 9.1 of the paper. The primary check is
deliberately graph-generic and uses no wheel-specific knowledge at all: build each
weighted wheel, run BFS from every source, enumerate every shortest path between
every unordered pair, and verify that the weights are exactly 1,...,t_gp(W_n),
once each. A bug in the paper's classification of wheel geodesics could therefore
not hide here.

Independent secondary checks cover:
  * the four geodesic classes of Proposition 2.1;
  * the Gallai graph identity Gamma(W_{m+1}) ~= complement(C_m) disjoint-union C_m;
  * the total-sum, square-sum, parity, and residual-sum identities of Appendix B;
  * the large-spoke lemma (Lemma 4.1): at most two spokes exceed N/2, and two
    only if their rim endpoints are adjacent.

No third-party package is required.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter, deque
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Set, Tuple

Node = str
PathTuple = Tuple[Node, ...]
Edge = frozenset[Node]


@dataclass(frozen=True)
class Certificate:
    n: int
    spokes: Tuple[int, ...]
    rims: Tuple[int, ...]

    @property
    def m(self) -> int:
        return self.n - 1

    @property
    def expected_geodesics(self) -> int:
        return self.m * (self.m + 3) // 2


def default_certificate_path() -> Path:
    return Path(__file__).resolve().parents[1] / "data" / "wheel_certificates_W5_W13.json"


def load_certificates(path: Path) -> Tuple[Certificate, ...]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    certificates = tuple(
        Certificate(
            n=int(item["n"]),
            spokes=tuple(int(x) for x in item["spokes"]),
            rims=tuple(int(x) for x in item["rims"]),
        )
        for item in raw["certificates"]
    )
    if not certificates:
        raise ValueError("certificate file contains no certificates")
    return certificates


def build_wheel(certificate: Certificate) -> tuple[
    Dict[Node, Set[Node]], Dict[Edge, int], List[Node]
]:
    m = certificate.m
    if len(certificate.spokes) != m or len(certificate.rims) != m:
        raise ValueError(f"W_{certificate.n} needs {m} spokes and {m} rim edges")

    all_labels = certificate.spokes + certificate.rims
    if any(label <= 0 for label in all_labels):
        raise ValueError("all edge labels must be positive")
    if len(set(all_labels)) != len(all_labels):
        raise ValueError("one-edge path weights repeat")

    hub = "h"
    rim = [f"v{i}" for i in range(m)]
    nodes = [hub, *rim]
    adjacency: Dict[Node, Set[Node]] = {node: set() for node in nodes}
    labels: Dict[Edge, int] = {}

    for i, vertex in enumerate(rim):
        adjacency[hub].add(vertex)
        adjacency[vertex].add(hub)
        labels[frozenset((hub, vertex))] = certificate.spokes[i]

    for i, vertex in enumerate(rim):
        successor = rim[(i + 1) % m]
        adjacency[vertex].add(successor)
        adjacency[successor].add(vertex)
        labels[frozenset((vertex, successor))] = certificate.rims[i]

    return adjacency, labels, nodes


def bfs_distances(adjacency: Mapping[Node, Set[Node]], source: Node) -> Dict[Node, int]:
    distances = {source: 0}
    queue: deque[Node] = deque([source])
    while queue:
        current = queue.popleft()
        for neighbor in adjacency[current]:
            if neighbor not in distances:
                distances[neighbor] = distances[current] + 1
                queue.append(neighbor)
    return distances


def all_shortest_paths(
    adjacency: Mapping[Node, Set[Node]], source: Node, target: Node
) -> List[PathTuple]:
    """Enumerate every shortest source-target path exactly once."""
    distances = bfs_distances(adjacency, source)
    target_distance = distances[target]
    paths: List[PathTuple] = []

    def dfs(current: Node, path: List[Node]) -> None:
        if current == target:
            if len(path) - 1 == target_distance:
                paths.append(tuple(path))
            return
        if distances[current] >= target_distance:
            return
        for neighbor in sorted(adjacency[current]):
            if distances.get(neighbor) == distances[current] + 1:
                dfs(neighbor, [*path, neighbor])

    dfs(source, [source])
    return paths


def path_weight(path: Sequence[Node], labels: Mapping[Edge, int]) -> int:
    return sum(
        labels[frozenset((left, right))]
        for left, right in zip(path, path[1:])
    )


def wheel_formula_weights(certificate: Certificate) -> List[int]:
    a = certificate.spokes
    b = certificate.rims
    m = certificate.m
    weights = [*a, *b]
    weights.extend(b[i] + b[(i + 1) % m] for i in range(m))

    cycle_edges = {
        tuple(sorted((i, (i + 1) % m)))
        for i in range(m)
    }
    weights.extend(
        a[i] + a[j]
        for i, j in combinations(range(m), 2)
        if (i, j) not in cycle_edges
    )
    return weights


def gallai_edges_of_wheel(certificate: Certificate) -> Set[Tuple[str, str]]:
    """Construct Gallai adjacency directly from pairs of wheel edges."""
    m = certificate.m
    wheel_edges: List[Tuple[str, Tuple[str, str]]] = []
    for i in range(m):
        wheel_edges.append((f"s{i}", ("h", f"v{i}")))
    for i in range(m):
        wheel_edges.append((f"r{i}", (f"v{i}", f"v{(i + 1) % m}")))

    original_edge_set = {
        frozenset(endpoints) for _, endpoints in wheel_edges
    }
    gallai: Set[Tuple[str, str]] = set()

    for (name1, e1), (name2, e2) in combinations(wheel_edges, 2):
        common = set(e1) & set(e2)
        if len(common) != 1:
            continue
        center = next(iter(common))
        outer1 = e1[0] if e1[1] == center else e1[1]
        outer2 = e2[0] if e2[1] == center else e2[1]
        spans_triangle = frozenset((outer1, outer2)) in original_edge_set
        if not spans_triangle:
            gallai.add(tuple(sorted((name1, name2))))
    return gallai


def expected_wheel_gallai_edges(m: int) -> Set[Tuple[str, str]]:
    expected: Set[Tuple[str, str]] = set()
    cycle_pairs = {
        tuple(sorted((i, (i + 1) % m))) for i in range(m)
    }
    for i, j in combinations(range(m), 2):
        if (i, j) not in cycle_pairs:
            expected.add(tuple(sorted((f"s{i}", f"s{j}"))))
    for i in range(m):
        expected.add(tuple(sorted((f"r{i}", f"r{(i + 1) % m}"))))
    return expected


def parity_transitions(values: Sequence[int]) -> int:
    return sum(
        (values[i] & 1) != (values[(i + 1) % len(values)] & 1)
        for i in range(len(values))
    )


def audit_identities(certificate: Certificate) -> None:
    a = certificate.spokes
    b = certificate.rims
    m = certificate.m
    N = certificate.expected_geodesics

    # Total-sum identity.
    lhs_sum = (m - 2) * sum(a) + 3 * sum(b)
    rhs_sum = N * (N + 1) // 2
    assert lhs_sum == rhs_sum, (certificate.n, "sum", lhs_sum, rhs_sum)

    # Residual set and its divisibility-by-3 identity.
    spoke_weights = set(a)
    cycle_pairs = {
        tuple(sorted((i, (i + 1) % m))) for i in range(m)
    }
    spoke_weights.update(
        a[i] + a[j]
        for i, j in combinations(range(m), 2)
        if (i, j) not in cycle_pairs
    )
    residual = set(range(1, N + 1)) - spoke_weights
    assert residual == set(b) | {b[i] + b[(i + 1) % m] for i in range(m)}
    assert len(residual) == 2 * m
    assert sum(residual) == 3 * sum(b)

    # Square-sum identity.
    lhs_sq = N * (N + 1) * (2 * N + 1) // 6
    rhs_sq = (
        sum(a) ** 2
        + (m - 3) * sum(x * x for x in a)
        - 2 * sum(a[i] * a[(i + 1) % m] for i in range(m))
        + 3 * sum(x * x for x in b)
        + 2 * sum(b[i] * b[(i + 1) % m] for i in range(m))
    )
    assert lhs_sq == rhs_sq, (certificate.n, "squares", lhs_sq, rhs_sq)

    # Parity identity.
    p = sum(x & 1 for x in a)
    r = sum(x & 1 for x in b)
    tau_a = parity_transitions(a)
    tau_b = parity_transitions(b)
    odd_count = p + r + p * (m - p) - tau_a + tau_b
    assert odd_count == (N + 1) // 2
    assert tau_a % 2 == 0 and tau_b % 2 == 0

    # Upper-half spoke structure.
    upper = [i for i, value in enumerate(a) if 2 * value > N]
    assert len(upper) <= 2
    if len(upper) == 2:
        i, j = sorted(upper)
        assert (j - i) % m in {1, m - 1}


def verify(certificate: Certificate, *, verbose: bool = False) -> List[Tuple[int, PathTuple]]:
    adjacency, labels, nodes = build_wheel(certificate)
    weighted_paths: List[Tuple[int, PathTuple]] = []

    for source, target in combinations(nodes, 2):
        for path in all_shortest_paths(adjacency, source, target):
            weighted_paths.append((path_weight(path, labels), path))

    weighted_paths.sort(key=lambda item: (item[0], item[1]))
    expected_count = certificate.expected_geodesics
    expected_weights = list(range(1, expected_count + 1))
    actual_weights = [weight for weight, _ in weighted_paths]

    if len(weighted_paths) != expected_count:
        raise AssertionError(
            f"W_{certificate.n}: enumerated {len(weighted_paths)} geodesics; "
            f"expected {expected_count}"
        )
    if actual_weights != expected_weights:
        counts = Counter(actual_weights)
        missing = [value for value in expected_weights if counts[value] == 0]
        repeated = sorted(value for value, count in counts.items() if count > 1)
        outside = sorted(value for value in counts if not 1 <= value <= expected_count)
        raise AssertionError(
            f"W_{certificate.n}: invalid certificate; "
            f"missing={missing}, repeated={repeated}, outside={outside}"
        )

    formula_weights = wheel_formula_weights(certificate)
    assert sorted(formula_weights) == expected_weights

    gallai = gallai_edges_of_wheel(certificate)
    expected_gallai = expected_wheel_gallai_edges(certificate.m)
    assert gallai == expected_gallai
    assert 2 * certificate.m + len(gallai) == expected_count

    audit_identities(certificate)

    if verbose:
        print(f"W_{certificate.n}: PASS ({expected_count} geodesics)")
        print(f"  spokes = {list(certificate.spokes)}")
        print(f"  rims   = {list(certificate.rims)}")
        for weight, path in weighted_paths:
            print(f"  {weight:>3}: {'-'.join(path)}")
    else:
        print(
            f"W_{certificate.n}: PASS — {expected_count} shortest paths have "
            f"weights 1..{expected_count} exactly once"
        )

    return weighted_paths


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--certificates",
        type=Path,
        default=default_certificate_path(),
        help="JSON certificate file",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="print every shortest path and its weight",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    certificates = load_certificates(args.certificates)
    for certificate in certificates:
        verify(certificate, verbose=args.verbose)
    print("PASS: all certificates, Gallai structures, and arithmetic identities verified.")


if __name__ == "__main__":
    main()
