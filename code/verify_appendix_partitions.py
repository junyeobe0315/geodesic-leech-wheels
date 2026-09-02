#!/usr/bin/env python3
"""Check Appendix A of the paper against the certificate data.

Appendix A prints, for each of W_7,...,W_13, the three sets

    S_1 = the one-edge geodesic weights          {a_i} u {b_i}
    S_2 = the two-rim-edge geodesic weights      {b_i + b_{i+1}}
    S_3 = the hub-type two-edge geodesic weights {a_i + a_j}, {i,j} not a rim edge

as explicit lists of integers. Those lists were typed into the LaTeX source by
hand, so they are exactly the kind of thing that drifts. This script parses them
straight out of paper/main.tex, recomputes the same three sets from
data/wheel_certificates_W5_W13.json, and compares them elementwise. It also
re-checks that the three sets are pairwise disjoint with union [N], which is the
content of Theorem 9.1.

No third-party package is required.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Sequence, Set, Tuple

REPO_ROOT = Path(__file__).resolve().parents[1]

WHEEL_HEADING = re.compile(r"\\subsection\*\{\$W_\{?(\d+)\}?\$\s*\(\$N=(\d+)\$\)\}")
SET_ENTRY = re.compile(r"\\mathcal\s*S_(\d)\s*=\s*\{\}\s*&")
INTEGER = re.compile(r"\d+")


def appendix_section(source: str) -> str:
    """The text of Appendix A, from its \\section to the next \\section."""
    start = source.rindex(r"\section{", 0, source.index(r"\label{app:partitions}"))
    tail = source.index(r"\section{", start + 1)
    return source[start:tail]


def parse_appendix(path: Path) -> Dict[int, Dict[int, List[int]]]:
    """Return {n: {1: S_1, 2: S_2, 3: S_3}} exactly as printed in Appendix A."""
    body = appendix_section(path.read_text(encoding="utf-8"))
    headings = list(WHEEL_HEADING.finditer(body))
    if not headings:
        raise SystemExit(f"no wheel subsections found in Appendix A of {path}")

    parsed: Dict[int, Dict[int, List[int]]] = {}
    for index, heading in enumerate(headings):
        n = int(heading.group(1))
        stop = headings[index + 1].start() if index + 1 < len(headings) else len(body)
        block = body[heading.end():stop]

        entries = list(SET_ENTRY.finditer(block))
        sets: Dict[int, List[int]] = {}
        for position, entry in enumerate(entries):
            which = int(entry.group(1))
            end = entries[position + 1].start() if position + 1 < len(entries) else len(block)
            payload = block[entry.end():end]
            # Drop the alignment and line-break markup, then take every integer.
            payload = payload.replace(r"\qquad", " ").replace(r"\\", " ")
            payload = payload.replace("&", " ").replace(r"\end{align*}", " ")
            sets[which] = [int(v) for v in INTEGER.findall(payload)]
        parsed[n] = sets
    return parsed


def expected_sets(spokes: Sequence[int], rims: Sequence[int]) -> Dict[int, Set[int]]:
    m = len(spokes)
    cycle = {tuple(sorted((i, (i + 1) % m))) for i in range(m)}
    return {
        1: set(spokes) | set(rims),
        2: {rims[i] + rims[(i + 1) % m] for i in range(m)},
        3: {
            spokes[i] + spokes[j]
            for i, j in combinations(range(m), 2)
            if (i, j) not in cycle
        },
    }


def load_certificates(path: Path) -> Dict[int, Tuple[List[int], List[int]]]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    return {
        item["n"]: (item["spokes"], item["rims"]) for item in raw["certificates"]
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tex",
        type=Path,
        default=REPO_ROOT / "paper" / "main.tex",
        help="LaTeX source containing Appendix A",
    )
    parser.add_argument(
        "--certificates",
        type=Path,
        default=REPO_ROOT / "data" / "wheel_certificates_W5_W13.json",
        help="JSON certificate file",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    printed = parse_appendix(args.tex)
    certificates = load_certificates(args.certificates)

    missing = sorted(set(certificates) - set(printed))
    extra = sorted(set(printed) - set(certificates))
    if missing or extra:
        raise SystemExit(
            f"Appendix A and the certificate file disagree on which wheels appear: "
            f"missing from the appendix {missing}, missing from the data {extra}"
        )

    failures: List[str] = []
    for n in sorted(certificates):
        spokes, rims = certificates[n]
        m = len(spokes)
        N = m * (m + 3) // 2
        expected = expected_sets(spokes, rims)
        problems: List[str] = []

        for which in (1, 2, 3):
            listed = printed[n].get(which)
            if listed is None:
                problems.append(f"S_{which} absent from Appendix A")
                continue
            if len(listed) != len(set(listed)):
                repeated = sorted({v for v in listed if listed.count(v) > 1})
                problems.append(f"S_{which} lists {repeated} more than once")
            if set(listed) != expected[which]:
                only_paper = sorted(set(listed) - expected[which])
                only_data = sorted(expected[which] - set(listed))
                problems.append(
                    f"S_{which} differs: in the paper only {only_paper}, "
                    f"in the data only {only_data}"
                )

        union: Set[int] = set()
        overlap: Set[int] = set()
        for which in (1, 2, 3):
            block = expected[which]
            overlap |= union & block
            union |= block
        if overlap:
            problems.append(f"the three classes overlap on {sorted(overlap)}")
        if union != set(range(1, N + 1)):
            problems.append(f"the union is not [1,{N}]")

        if problems:
            failures.append(f"W_{n}: " + "; ".join(problems))
            print(f"W_{n}: FAIL")
            for problem in problems:
                print(f"  - {problem}")
        else:
            sizes = ", ".join(f"|S_{k}|={len(expected[k])}" for k in (1, 2, 3))
            print(f"W_{n}: PASS — Appendix A matches the data ({sizes}, union = [1,{N}])")

    if failures:
        print(file=sys.stderr)
        raise SystemExit(f"FAILED {len(failures)} wheel(s):\n" + "\n".join(failures))
    print("PASS: every integer printed in Appendix A is confirmed by the certificates.")


if __name__ == "__main__":
    main()
