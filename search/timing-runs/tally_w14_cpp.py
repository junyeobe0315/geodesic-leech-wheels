#!/usr/bin/env python3
"""Tally the C++ W_14 campaign: per-seed outcome, steps, sets, orders, and totals."""
import re, glob, os
here = os.path.dirname(os.path.abspath(__file__))
rows = []
for path in sorted(glob.glob(os.path.join(here, "w14_cpp_seed*.log")), key=lambda p: int(re.search(r"seed(\d+)", p).group(1))):
    seed = int(re.search(r"seed(\d+)", path).group(1))
    text = open(path).read()
    final = re.search(r"^W_14\s+(PASS|not found)\s+(\d+(?:\.\d+)?)s, (\d+) feasible sets(?:, (\d+) orders tried)?", text, re.M)
    prog = re.findall(r"^\s+W_14:\s+(\d+)s\s+([\d,]+) steps\s+(\d+) feasible sets\s+(\d+) orders tried", text, re.M)
    steps = int(prog[-1][1].replace(",", "")) if prog else 0
    status, secs, sets, orders = (final.group(1), float(final.group(2)), int(final.group(3)), int(final.group(4) or 0)) if final else ("running", float(prog[-1][0]) if prog else 0.0, int(prog[-1][2]) if prog else 0, int(prog[-1][3]) if prog else 0)
    rows.append((seed, status, secs, steps, sets, orders))
print(f"{'seed':>4} {'status':>9} {'seconds':>8} {'steps':>18} {'feasible sets':>14} {'orders':>7}")
for r in rows:
    print(f"{r[0]:>4} {r[1]:>9} {r[2]:>8.0f} {r[3]:>18,} {r[4]:>14,} {r[5]:>7}")
T = [sum(r[i] for r in rows) for i in (2, 3, 4, 5)]
print(f"TOTAL seconds={T[0]:,.0f} ({T[0]/3600:.1f} CPU-h)  steps={T[1]:,} ({T[1]:.3e})  sets={T[2]:,} ({T[2]:.3e})  orders={T[3]}")
