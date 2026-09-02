#!/usr/bin/env python3
"""Run every reproducibility check in this repository.

    python3 code/run_all_checks.py

Each step is independent of the others and exits nonzero on failure, so this
script is also what CI runs.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import List, Sequence, Tuple

HERE = Path(__file__).resolve().parent

STEPS: Sequence[Tuple[str, str, Tuple[str, ...]]] = (
    (
        "certificates",
        "verify_wheel_certificates.py",
        (),
    ),
    (
        "appendix A",
        "verify_appendix_partitions.py",
        (),
    ),
    (
        "sections 8-9 numerics",
        "verify_boundary_numerics.py",
        (),
    ),
    (
        "lemma 6.3 small cases",
        "verify_finite_fourier_lemma_small.py",
        ("--max-w", "18"),
    ),
)


def run(script: str, args: Sequence[str]) -> None:
    command = [sys.executable, str(HERE / script), *args]
    print("\n$ " + " ".join(command), flush=True)
    subprocess.run(command, check=True, cwd=HERE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--quick",
        action="store_true",
        help="shrink the exponential Fourier sweep, for a fast smoke test",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    failed: List[str] = []
    for label, script, extra in STEPS:
        if args.quick and script == "verify_finite_fourier_lemma_small.py":
            extra = ("--max-w", "14")
        try:
            run(script, extra)
        except subprocess.CalledProcessError:
            failed.append(label)

    print()
    if failed:
        raise SystemExit(f"FAILED: {', '.join(failed)}")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
