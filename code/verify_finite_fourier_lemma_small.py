#!/usr/bin/env python3
"""Exhaustive small-range check of Lemma 5.3 (distinct pair sum upper bound).

This is not a proof.  For every interval length 2 <= w <= --max-w it enumerates
every subset S of {0,1,...,w-1} and checks the inequality of Lemma 5.3 wherever
its hypothesis L_w(r) >= 1/2 holds.  The proof itself is in Section 5 of the
paper.

Lemma 5.3 (paper, eq. (5.4)).  Let w >= 2, let S be a subset of {0,...,w-1} and
put r = |S|, M_w = 2w - 1, rho_{w-1} = 1 - 1/M_w and

    L_w(r) = 2 rho_{w-1} r / (pi + 2 rho_{w-1}).

Write s^x(S) = |{ x + y : x, y in S, x < y }| for the number of sums of
*distinct* element pairs.  If L_w(r) >= 1/2, then

    s^x(S)  <=  (M_w + C(r,2))/2  -  (L_w(r)^2 - L_w(r))/4.

Both sides are evaluated here exactly as stated: the left-hand side counts sums
of distinct pairs only (no doubles x + x) and is not truncated, and the modulus
is M_w = 2w - 1, which is odd, as the proof of the lemma requires.
"""
from __future__ import annotations

import argparse
import itertools
import math

TOL = 1e-10


def check_w(w: int) -> tuple[int, float]:
    """Check every S subseteq {0,...,w-1}; return (subsets tested, min slack)."""
    M_w = 2 * w - 1
    rho = 1.0 - 1.0 / M_w
    tested = 0
    minimum_slack = float("inf")

    for r in range(2, w + 1):
        L = 2.0 * rho * r / (math.pi + 2.0 * rho)
        if L < 0.5:
            continue
        rhs = (M_w + r * (r - 1) // 2) / 2.0 - (L * L - L) / 4.0
        for S in itertools.combinations(range(w), r):
            s_cross = len({a + b for i, a in enumerate(S) for b in S[i + 1:]})
            slack = rhs - s_cross
            minimum_slack = min(minimum_slack, slack)
            tested += 1
            if slack < -TOL:
                raise AssertionError(
                    f"violation at w={w}, S={S}, s^x={s_cross}, rhs={rhs}, slack={slack}"
                )

    return tested, minimum_slack


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--max-w",
        "--max-m",
        dest="max_w",
        type=int,
        default=18,
        help="largest interval length to exhaust (default: 18)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 2 <= args.max_w <= 24:
        raise SystemExit("--max-w must be between 2 and 24; the running time is exponential")

    total = 0
    for w in range(2, args.max_w + 1):
        tested, slack = check_w(w)
        shown = "inf" if slack == float("inf") else f"{slack:.12g}"
        print(f"w={w:2d}: M_w={2*w-1:3d}, tested={tested:8d}, minimum slack={shown}")
        total += tested
    print(f"PASS: Lemma 5.3 holds on all {total} subsets satisfying its hypothesis.")


if __name__ == "__main__":
    main()
