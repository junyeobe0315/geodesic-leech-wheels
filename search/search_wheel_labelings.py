#!/usr/bin/env python3
"""Search for geodesic Leech labelings of the wheel W_n.

This is an independent reimplementation written for the companion repository.
It is *not* the program that originally produced the labelings in Table 5 of the
paper; that code was not kept. On one core it finds valid labelings for
W_5,...,W_12 -- the first six of them in under a minute, W_12 in a few. W_13 it
also reaches, but only just: eight seeds in parallel turned one up after four
and a half minutes on a 2024 desktop core, some 37 CPU-minutes in all (the logs
are in timing-runs/). A C++ port with the same command line, about 87x faster
per step, is in search_wheel_labelings.cpp. It is offered as a starting point for the
open cases W_14,...,W_40 (the conclusion of the paper asks whether W_14 is
geodesic Leech); see "What limits it" below for how the cost stops paying at
that size. Nothing in the paper depends on this program: the existence claims are
certified by the finite integer partitions of Appendix A, which ../code/ checks
directly.

Method
------
Proposition 3.2 (frame-completion) splits the problem into choosing an admissible
spoke frame and completing it with a rim cycle over the 2m residual integers.
The two halves are wildly uneven, so they get different treatment.

*Phase 2, the rim completion, is free.* The residual set T = [N] \\ Phi(A) has
exactly 2m elements, and W_{m+1} is geodesic Leech precisely when T splits as
{b_i} disjoint-union {b_i + b_{i+1}} for a cyclic B. Exact backtracking settles
this in well under a millisecond for every frame in the paper: min(T) must be a
rim label, and once b_{k-1} is fixed the only admissible b_k are those with
b_{k-1} + b_k again in T.

*Phase 1, the spoke frame, is the whole problem.* It is a dense Sidon-type
packing -- m(m-3)/2 pair sums must fit without repetition into the m(m+1)/2 slots
the labels leave free. Two things make it tractable here.

  Search the set, not the sequence. A set S of m labels underlies an admissible
  frame under *some* cyclic order iff there is a Hamiltonian cycle C on S with
  every pair outside C having a distinct sum in [1,N] \\ S, that is, iff

    (i)  every pair whose sum exceeds N or lands in S is an edge of C, and
    (ii) for each value v, at most one pair summing to v lies outside C.

  Writing total_v for the number of pairs summing to v and f_v for how many of
  those are forced by (i), condition (ii) says C must hold at least
  max(f_v, total_v - 1) of each group, so required = sum_v max(f_v, total_v - 1)
  cannot exceed m = |C|, and no label may carry more than two forced pairs.
  Those two necessary conditions are the annealing cost. Dropping the cyclic
  order shrinks the search space by a factor of (m-1)!/2, and it shows: measured
  against the same annealer run over sequences instead, this one produces
  feasible states two orders of magnitude faster at m = 10, and still produces
  them at m = 13 -- the frame size for W_14 -- where the sequence search produced
  none even at m = 12.

  Anneal, do not branch. Depth-first search commits to its first few labels and
  burns its whole budget without backtracking to them. Simulated annealing with
  a cycled temperature keeps emitting distinct feasible sets from one run.

A feasible set is then assembled into cyclic orders: the forced pairs form
disjoint paths, the collision groups are topped up to their quota, and the paths
are linked into a Hamiltonian cycle. Different completions give different Phi,
hence different T, so one set yields several independent Phase 2 attempts.

The sum identity of Appendix B enters the cost as well. sum(T) = 3 sum(b) and
sum(Phi) = (m-2) sum(A) depend only on the set, so sets with the wrong residue
mod 3 can never complete and are penalized rather than generated and discarded.

What limits it
--------------
The cost is a *necessary* condition, not a sufficient one, and the gap widens
with m. Whether the quota edges can actually be placed alongside the forced ones
-- degree at most two, no premature cycle -- is a small simultaneous-placement
problem that the cost only approximates, so a cost-zero set need not assemble.
The share of feasible sets that assemble falls from roughly 1% at m = 10 to 0.2%
at m = 11, and at m = 12 -- where W_13 sits -- a five-minute run turns up about
500 feasible sets of which 4 assemble, yielding 30 rim attempts, one every ten
seconds. That is the whole reason W_13 costs hours while W_12 costs minutes.

Two attempts to close the gap by charging for a closed cycle in the forced graph,
and for collision groups with no room left, were both tried and both measured
*worse*: they cost more per move than they saved in wasted assemblies. Raising
--selections does not help either -- at m = 12, 8 and 128 give the identical 4
assembled sets and 30 orders, because the sets that assemble at all have few
collision groups to choose within. The lever that would matter is a cost
reflecting the placement problem itself rather than bounding its pieces
separately.

Runs are deterministic given --seed, so any reported labeling is reproducible.
The search is incomplete by design: finding nothing is not evidence of
nonexistence, and the script says so.

No third-party package is required.
"""
from __future__ import annotations

import argparse
import json
import math
import random
import sys
import time
from itertools import combinations
from pathlib import Path
from typing import Callable, Dict, Iterator, List, Optional, Sequence, Set, Tuple

Pair = Tuple[int, int]
Labeling = Tuple[Tuple[int, ...], Tuple[int, ...]]


def geodesic_count(m: int) -> int:
    """t_gp(W_{m+1}) = m(m+3)/2."""
    return m * (m + 3) // 2


# ---------------------------------------------------------------------------
# Phase 1a: annealing over spoke sets
# ---------------------------------------------------------------------------

class SpokeSet:
    """A set of m labels, with the cost described in the module docstring.

    ``cost`` is zero exactly when the two necessary conditions hold -- the
    required-edge count fits inside a Hamiltonian cycle and no label carries more
    than two forced pairs -- and the residue condition of Appendix B is met.
    Everything is maintained incrementally: replacing one label touches the m-1
    pairs at that label, plus any pair whose sum happens to equal the label
    leaving or entering the set.
    """

    def __init__(self, m: int, N: int, values: Sequence[int]) -> None:
        self.m = m
        self.N = N
        self.s = list(values)
        self.in_set = bytearray(N + 1)
        for value in self.s:
            self.in_set[value] = 1

        self.total: Dict[int, int] = {}          # sum -> pairs with that sum
        self.forced_count: Dict[int, int] = {}   # sum -> forced pairs with that sum
        self.pairs_of: Dict[int, Set[Pair]] = {}
        self.is_forced: Dict[Pair, bool] = {}
        self.forced_degree = [0] * m
        self.required = 0

        self.label_sum = sum(self.s)
        self.residue_target = (N * (N + 1) // 2) % 3

        for i, j in combinations(range(m), 2):
            self._add_pair(i, j)
        self.cost = self._recompute_cost()

    # -- the "required" running total -------------------------------------
    def _quota(self, total: int) -> int:
        return max(0, total - 1)

    def _required_of(self, value: int) -> int:
        return max(self.forced_count.get(value, 0), self._quota(self.total.get(value, 0)))

    def _add_pair(self, i: int, j: int) -> None:
        value = self.s[i] + self.s[j]
        before = self._required_of(value)
        self.total[value] = self.total.get(value, 0) + 1
        self.pairs_of.setdefault(value, set()).add((i, j))
        forced = value > self.N or bool(self.in_set[value])
        self.is_forced[(i, j)] = forced
        if forced:
            self.forced_count[value] = self.forced_count.get(value, 0) + 1
            self._link(i, j, +1)
        self.required += self._required_of(value) - before

    def _remove_pair(self, i: int, j: int) -> None:
        value = self.s[i] + self.s[j]
        before = self._required_of(value)
        self.total[value] -= 1
        self.pairs_of[value].discard((i, j))
        if self.is_forced.pop((i, j)):
            self.forced_count[value] -= 1
            self._link(i, j, -1)
        self.required += self._required_of(value) - before

    def _refresh_forced(self, value: int) -> None:
        """Membership of ``value`` changed, so pairs summing to it may flip."""
        if value > self.N:
            return
        should_be = bool(self.in_set[value])
        for (i, j) in list(self.pairs_of.get(value, ())):
            if self.is_forced[(i, j)] == should_be:
                continue
            before = self._required_of(value)
            self.is_forced[(i, j)] = should_be
            step = 1 if should_be else -1
            self.forced_count[value] = self.forced_count.get(value, 0) + step
            self._link(i, j, step)
            self.required += self._required_of(value) - before

    def _link(self, i: int, j: int, step: int) -> None:
        self.forced_degree[i] += step
        self.forced_degree[j] += step

    def _recompute_cost(self) -> int:
        over_budget = max(0, self.required - self.m)
        over_degree = sum(max(0, d - 2) for d in self.forced_degree)
        residue = 0 if ((self.m - 2) * self.label_sum) % 3 == self.residue_target else 1
        return over_budget + over_degree + residue

    # -- the one move ------------------------------------------------------
    def replace(self, index: int, value: int) -> int:
        """Swap s[index] out for ``value``, which must not already be present."""
        previous = self.s[index]
        for j in range(self.m):
            if j != index:
                self._remove_pair(*sorted((index, j)))
        self.in_set[previous] = 0
        self.s[index] = value
        self.in_set[value] = 1
        for j in range(self.m):
            if j != index:
                self._add_pair(*sorted((index, j)))
        self._refresh_forced(previous)
        self._refresh_forced(value)
        self.label_sum += value - previous
        self.cost = self._recompute_cost()
        return previous


def anneal_spoke_sets(
    m: int,
    rng: random.Random,
    deadline: float,
    on_set: Callable[[SpokeSet], bool],
    t_hot: float = 2.0,
    t_cold: float = 0.05,
    cycle: int = 20000,
    progress: Optional[Callable[[int, int, int], None]] = None,
) -> Tuple[int, int]:
    """Stream feasible spoke sets to ``on_set`` until it returns True or time runs out.

    Returns (feasible sets produced, annealing steps taken).
    """
    N = geodesic_count(m)
    state = SpokeSet(m, N, rng.sample(range(1, N + 1), m))
    started = time.monotonic()
    produced = 0
    steps = 0
    ratio = t_cold / t_hot

    while True:
        if time.monotonic() - started >= deadline:
            return produced, steps
        if progress is not None:
            progress(produced, steps, state.cost)

        for _ in range(200):
            steps += 1
            # Cycle the temperature rather than cooling once, so a single run
            # keeps visiting new basins after emptying the current one.
            temperature = t_hot * ratio ** ((steps % cycle) / cycle)

            index = rng.randrange(m)
            value = rng.randint(1, N)
            if state.in_set[value]:
                continue
            before = state.cost
            previous = state.replace(index, value)
            rise = state.cost - before
            if rise > 0 and rng.random() > math.exp(-rise / temperature):
                state.replace(index, previous)

            if state.cost == 0:
                produced += 1
                if on_set(state):
                    return produced, steps
                # nudge off this set so the next one is different
                index = rng.randrange(m)
                value = rng.randint(1, N)
                if not state.in_set[value]:
                    state.replace(index, value)


# ---------------------------------------------------------------------------
# Phase 1b: assembling a feasible set into cyclic orders
# ---------------------------------------------------------------------------

class PathSet:
    """Edges forming vertex-disjoint simple paths on 0..m-1, with cheap undo."""

    __slots__ = ("m", "adjacency")

    def __init__(self, m: int) -> None:
        self.m = m
        self.adjacency: List[List[int]] = [[] for _ in range(m)]

    def _reaches(self, i: int, j: int) -> bool:
        previous, current = -1, i
        while True:
            forward = [x for x in self.adjacency[current] if x != previous]
            if not forward:
                return False
            previous, current = current, forward[0]
            if current == j or current == i:
                return True               # j is on i's path, or it is a cycle

    def can_add(self, i: int, j: int) -> bool:
        return (
            len(self.adjacency[i]) < 2
            and len(self.adjacency[j]) < 2
            and not self._reaches(i, j)
        )

    def add(self, i: int, j: int) -> None:
        self.adjacency[i].append(j)
        self.adjacency[j].append(i)

    def remove(self, i: int, j: int) -> None:
        self.adjacency[i].remove(j)
        self.adjacency[j].remove(i)

    def components(self) -> List[List[int]]:
        """The paths, each read end to end; isolated labels count as paths."""
        seen = [False] * self.m
        out: List[List[int]] = []
        for start in range(self.m):
            if seen[start] or len(self.adjacency[start]) == 2:
                continue                  # interior of a path; start from an end
            path = [start]
            seen[start] = True
            previous, current = -1, start
            while True:
                forward = [x for x in self.adjacency[current] if x != previous]
                if not forward:
                    break
                previous, current = current, forward[0]
                seen[current] = True
                path.append(current)
            out.append(path)
        return out


def edge_selections(state: SpokeSet, limit: int) -> Iterator[PathSet]:
    """Yield edge sets that every admissible cyclic order of ``state`` must contain.

    The forced pairs go in unconditionally; each collision group then hands over
    all but one of its members. Which ones is a choice, and choosing by
    backtracking -- with the degree bound and acyclicity maintained as we go --
    never wastes an attempt, where sampling them at random mostly does.
    """
    m = state.m
    paths = PathSet(m)
    for (i, j), forced in state.is_forced.items():
        if not forced:
            continue
        if not paths.can_add(i, j):
            return
        paths.add(i, j)

    quotas: List[Tuple[int, List[Pair]]] = []
    for value, total in state.total.items():
        deficit = (total - 1) - state.forced_count.get(value, 0)
        if deficit > 0:
            free = sorted(p for p in state.pairs_of[value] if not state.is_forced[p])
            if len(free) < deficit:
                return                    # the cost function should have caught this
            quotas.append((deficit, free))
    quotas.sort(key=lambda quota: len(quota[1]))

    produced = 0

    def descend(k: int) -> Iterator[PathSet]:
        nonlocal produced
        if produced >= limit:
            return
        if k == len(quotas):
            produced += 1
            yield paths
            return
        deficit, free = quotas[k]
        for combination in combinations(free, deficit):
            added: List[Pair] = []
            feasible = True
            for i, j in combination:
                if not paths.can_add(i, j):
                    feasible = False
                    break
                paths.add(i, j)
                added.append((i, j))
            if feasible:
                yield from descend(k + 1)
            for i, j in reversed(added):
                paths.remove(i, j)
            if produced >= limit:
                return

    yield from descend(0)


def cyclic_orders(
    state: SpokeSet,
    rng: random.Random,
    selections: int = 8,
    linkings: int = 4,
) -> Iterator[List[int]]:
    """Yield cyclic orders of ``state.s`` that make it an admissible spoke frame.

    Each selection leaves vertex-disjoint paths; linking those into a Hamiltonian
    cycle can only exempt further pairs, so it never breaks admissibility, and
    different linkings give different Phi and hence different residual sets.
    """
    for paths in edge_selections(state, selections):
        components = paths.components()
        for _ in range(linkings):
            order: List[int] = []
            shuffled = list(components)
            rng.shuffle(shuffled)
            for path in shuffled:
                order.extend(path[::-1] if rng.random() < 0.5 else path)
            if len(order) != state.m:
                break                     # should not happen; leave this selection
            yield [state.s[i] for i in order]


# ---------------------------------------------------------------------------
# Phase 2: rim completion
# ---------------------------------------------------------------------------

def residual_set(spokes: Sequence[int]) -> List[int]:
    """T(A) = [N] \\ Phi(A), of size 2m for an admissible frame."""
    m = len(spokes)
    N = geodesic_count(m)
    cycle = {tuple(sorted((i, (i + 1) % m))) for i in range(m)}
    phi = set(spokes)
    phi.update(
        spokes[i] + spokes[j]
        for i, j in combinations(range(m), 2)
        if (i, j) not in cycle
    )
    return sorted(set(range(1, N + 1)) - phi)


def complete_rim(residual: Sequence[int], m: int) -> Optional[List[int]]:
    """Split T into {b_i} and {b_i + b_{i+1}} for a cyclic B, if possible."""
    if len(residual) != 2 * m:
        return None
    total = sum(residual)
    if total % 3:
        return None                      # sum(T) = 3 sum(b), Appendix B
    rim_sum = total // 3

    available = set(residual)
    smallest = min(residual)             # min(T) cannot be a two-edge sum
    rim: List[int] = [smallest]
    available.discard(smallest)

    def descend(k: int, running: int) -> Optional[List[int]]:
        if running > rim_sum:
            return None
        if k == m:
            closing = rim[-1] + rim[0]
            if len(available) == 1 and closing in available:
                if m < 3 or rim[1] < rim[m - 1]:
                    return list(rim)
            return None
        previous = rim[-1]
        # Only values whose sum with the previous rim label is itself residual.
        for value in sorted(available):
            link = previous + value
            if link not in available:
                continue
            available.discard(value)
            available.discard(link)
            rim.append(value)
            found = descend(k + 1, running + value)
            rim.pop()
            available.add(value)
            available.add(link)
            if found is not None:
                return found
        return None

    return descend(1, smallest)


# ---------------------------------------------------------------------------
# Independent check of a finished labeling
# ---------------------------------------------------------------------------

def is_geodesic_leech(spokes: Sequence[int], rims: Sequence[int]) -> bool:
    """Recompute the four geodesic classes from scratch and compare with [N]."""
    m = len(spokes)
    if len(rims) != m:
        return False
    N = geodesic_count(m)
    cycle = {tuple(sorted((i, (i + 1) % m))) for i in range(m)}
    weights = list(spokes) + list(rims)
    weights += [rims[i] + rims[(i + 1) % m] for i in range(m)]
    weights += [
        spokes[i] + spokes[j]
        for i, j in combinations(range(m), 2)
        if (i, j) not in cycle
    ]
    return sorted(weights) == list(range(1, N + 1))


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def search(
    n: int,
    seed: int,
    time_limit: float,
    selections_per_set: int = 8,
    verbose: bool = False,
) -> Tuple[Optional[Labeling], dict]:
    m = n - 1
    rng = random.Random(seed)
    stats = {"sets": 0, "orders": 0, "steps": 0, "seconds": 0.0}
    solution: List[Optional[Labeling]] = [None]
    started = time.monotonic()
    last_report = [started]

    def on_set(state: SpokeSet) -> bool:
        stats["sets"] += 1
        seen: Set[Tuple[int, ...]] = set()
        for order in cyclic_orders(state, rng, selections_per_set):
            key = tuple(order)
            if key in seen:
                continue
            seen.add(key)
            stats["orders"] += 1
            residual = residual_set(order)
            rim = complete_rim(residual, m)
            if rim is None or not is_geodesic_leech(order, rim):
                continue
            solution[0] = (tuple(order), tuple(rim))
            return True
        return False

    def report_progress(sets: int, steps: int, cost: int) -> None:
        now = time.monotonic()
        if now - last_report[0] < 15.0:
            return
        last_report[0] = now
        print(
            f"    W_{n}: {now - started:5.0f}s  {steps:>10,} steps  "
            f"{sets:>7} feasible sets  {stats['orders']:>8} orders tried  "
            f"current cost {cost}",
            file=sys.stderr,
        )

    _, steps = anneal_spoke_sets(
        m, rng, time_limit, on_set,
        progress=report_progress if verbose else None,
    )
    stats["steps"] = steps
    stats["seconds"] = time.monotonic() - started
    return solution[0], stats


def canonical_form(spokes: Sequence[int], rims: Sequence[int]) -> Tuple[int, ...]:
    """Least representative of a labeling under rotation and reflection of the rim.

    A rotation by r sends (a_k, b_k) to (a_{k+r}, b_{k+r}); the reflection sends
    them to (a_{-k}, b_{-k-1}), since reversing the rim turns the edge
    v_k v_{k+1} into v_{-k-1} v_{-k}. Two labelings related this way are the same
    labeling of the same wheel.
    """
    m = len(spokes)
    images = []
    for r in range(m):
        images.append(
            tuple(spokes[(k + r) % m] for k in range(m))
            + tuple(rims[(k + r) % m] for k in range(m))
        )
        images.append(
            tuple(spokes[(r - k) % m] for k in range(m))
            + tuple(rims[(r - k - 1) % m] for k in range(m))
        )
    return min(images)


def published_certificates() -> dict:
    path = Path(__file__).resolve().parents[1] / "data" / "wheel_certificates_W5_W13.json"
    if not path.exists():
        return {}
    raw = json.loads(path.read_text(encoding="utf-8"))
    return {item["n"]: item for item in raw["certificates"]}


def report(n: int, result: Optional[Labeling], stats: dict, known: dict) -> bool:
    if result is None:
        print(
            f"W_{n:<3} not found   {stats['seconds']:.0f}s, "
            f"{stats['sets']} feasible sets, {stats['orders']} orders tried"
        )
        return False

    spokes, rims = result
    verified = is_geodesic_leech(spokes, rims)
    same = ""
    if n in known:
        matches = canonical_form(spokes, rims) == canonical_form(
            known[n]["spokes"], known[n]["rims"]
        )
        same = (
            "  (Table 5's labeling, up to rotation/reflection)"
            if matches
            else "  (a labeling not in Table 5)"
        )
    print(
        f"W_{n:<3} {'PASS' if verified else 'INVALID'}      "
        f"{stats['seconds']:.1f}s, {stats['sets']} feasible sets{same}"
    )
    print(f"      A = {list(spokes)}")
    print(f"      B = {list(rims)}")
    return verified


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-n",
        type=int,
        nargs="*",
        default=list(range(7, 14)),
        help="wheel orders to search (default: 7 8 9 10 11 12 13)",
    )
    parser.add_argument("--seed", type=int, default=20260831, help="RNG seed")
    parser.add_argument(
        "--time-limit",
        type=float,
        default=120.0,
        help="seconds to spend per wheel before giving up (default: 120)",
    )
    parser.add_argument(
        "--selections",
        type=int,
        default=8,
        help="edge selections to explore per feasible set (default: 8)",
    )
    parser.add_argument("--verbose", action="store_true", help="log progress")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    known = published_certificates()
    print(f"seed={args.seed}, time limit {args.time_limit:g}s per wheel\n")

    outcomes = []
    for n in args.n:
        if n < 5:
            raise SystemExit("the wheel W_n is defined here for n >= 5")
        result, stats = search(
            n, args.seed + n, args.time_limit, args.selections, args.verbose
        )
        outcomes.append(report(n, result, stats, known))

    found = sum(1 for ok in outcomes if ok)
    print(f"\n{found}/{len(outcomes)} wheels solved and independently verified.")
    if found < len(outcomes):
        print(
            "A wheel not solved here is undecided, not shown to be impossible: "
            "this search is incomplete by design."
        )


if __name__ == "__main__":
    main()
