# A geodesic Leech wheel has at most 40 vertices

[![checks](https://github.com/junyeobe0315/geodesic-leech-wheels/actions/workflows/checks.yml/badge.svg)](https://github.com/junyeobe0315/geodesic-leech-wheels/actions/workflows/checks.yml)
![python](https://img.shields.io/badge/python-3.9%2B-blue)
![c++17](https://img.shields.io/badge/C%2B%2B-C%2B%2B17-blue)
![category](https://img.shields.io/badge/category-combinatorics_(math.CO)-blueviolet)
![reviewed](https://img.shields.io/badge/peer%20reviewed-no-orange)
[![arXiv](https://img.shields.io/badge/arXiv-2609.02544-b31b1b)](https://arxiv.org/abs/2609.02544)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.22254583-blue)](https://doi.org/10.5281/zenodo.22254583)

Companion repository for the paper

> **A finiteness theorem for geodesic Leech wheels**
> — Junyeop Yim, [arXiv:2609.02544](https://arxiv.org/abs/2609.02544) (2026).
> [paper/main.pdf](paper/main.pdf) · source: [paper/](paper/)

A positive-integer edge labeling of a graph is a *geodesic Leech labeling* if the
weights of its geodesics are exactly 1, 2, …, `t_gp(G)`, each occurring once.
Lakshmanan S. and Manattu ([arXiv:2502.16628](https://arxiv.org/abs/2502.16628))
exhibited such labelings of the wheels `W_5` and `W_6`, and stated as an open
problem their belief that every `W_n` with `n ≥ 7` is a non-geodesic Leech graph.
They singled out `W_7`, for which they gave an *almost* geodesic Leech labeling
and wrote that the geodesic Leech case remained to be seen. The paper settles
that problem in the negative — `W_7` included — and bounds the phenomenon from
above:

> **Theorem.** If `n ≥ 5` and `W_n` is geodesic Leech, then `n ≤ 40`.
>
> **Corollary.** With `E = { n ≥ 5 : W_n is geodesic Leech }`,
>
> ```
> {5, 6, ..., 13}  ⊆  E  ⊆  {5, 6, ..., 40}
> ```

Here `W_n` is the wheel on `n` vertices **in total** — a hub joined to a cycle of
length `n − 1`. The cases `W_14, ..., W_40` remain open, and since nonexistence at
one order is not known to imply nonexistence at larger orders, they have to be
decided one at a time.

The upper bound is proved by hand, using a finite Fourier kernel and a
six-variable Parseval argument. This repository exists to make that proof easy to
audit, not to supply any part of it.

---

## What is here

| path | what it is |
|---|---|
| [`paper/main.tex`](paper/main.tex) | the article. Self-contained: inline bibliography, no BibTeX pass. `cd paper && make` |
| [`data/wheel_certificates_W5_W13.json`](data/wheel_certificates_W5_W13.json) | the labelings of `W_5 … W_13` from Table 5, machine-readable (`W_5`, `W_6` are Lakshmanan S.–Manattu’s; the rest are new). See [`data/README.md`](data/README.md) for the schema |
| [`code/`](code/) | the reproducibility checks, described below |
| [`search/`](search/) | a labeling search for the open cases `W_14 … W_40`: the reference implementation in Python, a C++ port of the same search with the same command line, and in [`search/timing-runs/`](search/timing-runs/) the logs of the runs reported in the paper |

### The checks

| script | what it certifies |
|---|---|
| [`verify_wheel_certificates.py`](code/verify_wheel_certificates.py) | **Theorem 9.1.** Builds each weighted wheel, runs BFS from every vertex, enumerates every shortest path between every pair, and confirms the weights are `1..t_gp` exactly once. Deliberately graph-generic — the wheel-specific formula of Proposition 2.1, the Gallai-graph identity (a cross-check internal to this script; it does not appear in the paper), and the arithmetic identities of Appendix B are only cross-checks. |
| [`verify_appendix_partitions.py`](code/verify_appendix_partitions.py) | **Appendix A.** Parses the printed geodesic weight classes straight out of `paper/main.tex` and compares them elementwise against the certificate data, so a transcription slip in the LaTeX cannot go unnoticed. |
| [`verify_boundary_numerics.py`](code/verify_boundary_numerics.py) | **Sections 5–8 and Appendices C–D.** Re-derives the occupancy bound and the two position formulas in exact integer arithmetic, then every decimal comparison behind `n ≤ 40` in exact rational interval arithmetic — no floating point anywhere. See below. |
| [`verify_finite_fourier_lemma_small.py`](code/verify_finite_fourier_lemma_small.py) | **Lemma 5.3**, on small cases. Exhausts every `S ⊆ {0,…,w−1}` for `w ≤ 18` and checks the lemma's own inequality — distinct-pair sums on the left, odd modulus `M_w = 2w−1` — wherever its hypothesis `L_w(ℓ) ≥ 1/2` holds. A sanity check, not a proof; the proof is in §5. |
| [`interval.py`](code/interval.py) | the rational interval arithmetic the numeric audit is built on: exact `Fraction` endpoints, outward rounding on every operation, certified integer square roots. |

`verify_boundary_numerics.py` is the substantial one. The nonexistence argument
ends in a few dozen explicit decimals, some of them tight — the case
`(m, h) = (42, 2)` closes with a margin of `0.0072` and `(40, 2)` with `0.0069`.
The script rebuilds all of them on exact `Fraction` intervals with outward
rounding, including the constant `Ψ` of Section 8, its full six-variable gradient
from the closed forms of Appendix D, the three tangent planes, and the `S_5`
lower bounds for `Λ`. It also certifies its own `π` bounds from Machin's formula
rather than assuming them, and — where the paper reaches a bound *through* a
rational bracket for `π`, as it does for `c*_m` — it checks that the stated
bracket actually delivers it, not merely that an exact `π` would. The two are
not the same claim, and only the first is what the proof needs — even though, as
it happens, `c*_m` is insensitive to the substitution: its elasticity in `π` is
about `-0.029`, so the bracket `355/113` moves `c*_40` from `78.868916838` to
`78.868916647`, damping π's relative error by a factor of roughly 35 rather than
amplifying it.

The same script also verifies Section 6, and for the same reason. The three
boundary cases are closed by comparing an upper bound for `Ψ` against a lower
bound for `Λ` that depends on the position bound `p ≤ Π_m` of Lemma 6.2, and
`(40, 2)` has **no slack at all**: `p = 101` closes it and `p = 102` does not.
So the script takes the position inequality to its fixed point separately for
each `d` — the constraint `d < κ/2` is what rules out the large-`d` branch — and
confirms `κ ≤ 3m` and `p ≤ Π_m` over a long range of `m`, together with the
auxiliary inequalities and the two asymptotic square-root arguments. It also
confirms the paper's own remark that near `w = 110` the relaxation
`λ_w ≥ 3/8` is exactly the difference between `μ(w) ≤ 21` and `μ(w) ≤ 22`, and
would leave only `p ≤ 102` at `m = 40`.

## Reproducing

Python ≥ 3.9. No third-party package, no configuration.

```bash
python3 code/run_all_checks.py
```

Expected tail:

```
ALL CHECKS PASSED
```

To compile the paper (needs `latexmk` and a TeX distribution):

```bash
cd paper && make
```

Everything is also run on each push by [`.github/workflows/checks.yml`](.github/workflows/checks.yml).

## Searching for new labelings

[`search/search_wheel_labelings.py`](search/search_wheel_labelings.py) is built on
the frame-completion criterion (Proposition 3.2), which splits the problem exactly
into choosing an admissible spoke frame `A` and completing it with a rim cycle
over the `2m` residual integers. The two halves turn out to be wildly uneven:

- **Phase 2 (rim completion) is free.** It is exact backtracking, and it recovers
  a valid rim for every frame in the paper in well under a millisecond. `min(T)`
  must be a rim label, and once `b_{k−1}` is fixed the only admissible `b_k` are
  those with `b_{k−1} + b_k` again in `T`, which leaves almost nothing to branch on.
- **Phase 1 (the spoke frame) is the whole problem.** It is a dense Sidon-type
  packing: `m(m−3)/2` pair sums must fit, without repetition, into the `m(m+1)/2`
  slots the labels leave free.

Two things make Phase 1 tractable.

**Search the set, not the sequence.** A set `S` of `m` labels underlies an
admissible frame under *some* cyclic order exactly when there is a Hamiltonian
cycle `C` on `S` with every pair outside `C` having a distinct sum in `[1,N] \ S`
— that is, when (i) every pair whose sum exceeds `N` or lands in `S` is an edge
of `C`, and (ii) for each value `v`, at most one pair summing to `v` lies outside
`C`. With `total_v` pairs summing to `v` and `f_v` of them forced by (i),
condition (ii) says `C` must hold `max(f_v, total_v − 1)` of each group, so

```
required = Σ_v max(f_v, total_v − 1)  ≤  m = |C|,
```

and no label may carry more than two forced pairs. Those two necessary
conditions are the cost. Dropping the cyclic order shrinks the space by
`(m−1)!/2`, and it shows: measured against the same annealer run over sequences
instead, this one produces feasible states about 240× faster at `m = 10`, and it
still produces them at `m = 13` — the frame size for `W_14` — where the sequence
search produced none even at `m = 12`.

**Anneal, don't branch.** Depth-first search commits to its first few labels and
spends its whole budget without ever backtracking to them. Simulated annealing
with a cycled temperature keeps emitting distinct feasible sets from one run.

A feasible set is then assembled into cyclic orders. The forced pairs go in
unconditionally; each collision group hands over all but one of its members, and
*which* ones is a choice made by backtracking with the degree bound and
acyclicity maintained as it goes — sampling them at random instead wastes most
attempts, because the cost function bounds only the *forced* degree. The
resulting disjoint paths are then linked into a Hamiltonian cycle; linking edges
only exempt further pairs, so they can never break admissibility. Different
completions give different `Φ`, hence different `T`, so one set yields several
independent Phase 2 attempts. The sum
identity of Appendix B is folded into the cost as well: `sum(Φ) = (m−2)·sum(A)`
and `sum(T) = 3·sum(b)` depend only on the set, so sets with the wrong residue
mod 3 can never complete and are penalized rather than generated and discarded.

```bash
python3 search/search_wheel_labelings.py -n 7 8 9
python3 search/search_wheel_labelings.py -n 14 --time-limit 3600 --verbose
```

Runs are deterministic given `--seed`, and every labeling found is re-verified
from the graph before being printed. `W_14` is left open in Section 10 of the paper; a run that finds nothing leaves it undecided, and the script says so.

### What it reaches, and what stops it

Measured on one core of an Intel Core i5-14400F (WSL2, Python 3.13), default
seed; the full logs are in [`search/timing-runs/`](search/timing-runs/) and the
per-run rows in [`summary.tsv`](search/timing-runs/summary.tsv):

| wheel | `W_5`–`W_9` | `W_10` | `W_11` | `W_12` | `W_13` |
|---|---|---|---|---|---|
| time to a verified labeling | under a second | 9.4 s | 47 s | 232 s | 275 s wall, eight seeds in parallel (≈ 37 CPU-minutes) |

Single core with the default seed, except `W_13`: no single seed was run to
completion there; seeds 1–8 were raced in parallel and seed 2 landed after 275
seconds. Wall-clock figures move with load and hardware; the labelings
themselves are fixed by the seed and reproduce exactly.

The jump from minutes to hours has a specific cause. The cost is a *necessary*
condition, not a sufficient one, and the gap widens with `m`. Whether the quota
edges can actually be placed alongside the forced ones — degree at most two, no
premature cycle — is a small simultaneous-placement problem that the cost only
approximates, so a cost-zero set need not assemble. The share that does assemble
falls from roughly 1% at `m = 10` to 0.2% at `m = 11`, and at `m = 12` — where
`W_13` sits — a five-minute run turns up about 500 feasible sets of which 4
assemble, giving 30 rim attempts, one every ten seconds.

Three things were tried against that and measured, all negative: charging for a
closed cycle in the forced graph, charging for collision groups with no room
left, and raising `--selections` (at `m = 12`, 8 and 128 give the identical 4
assembled sets and 30 orders, because the sets that assemble at all have few
collision groups to choose within). The first two cost more per move than they
saved in wasted assemblies. The lever that would matter is a cost reflecting the
placement problem itself rather than bounding its pieces separately.

### The C++ engine, and the `W_14` campaign

[`search_wheel_labelings.cpp`](search/search_wheel_labelings.cpp) is a port of
the same search — same cost, same moves, same step-based temperature cycle,
same assembly and rim completion, same command line and output — for the long
open cases, where interpreter overhead is the whole cost:

```bash
g++ -O3 -march=native -std=c++17 -o search/wheel_search_cpp search/search_wheel_labelings.cpp
./search/wheel_search_cpp -n 14 --seed 1 --time-limit 64800 --verbose
```

It runs about 2.0 million annealing steps per second against CPython's 23
thousand on the same core, roughly 87×. The two programs use different random
streams, so a given seed does not trace the same trajectory in both; but on
identical inputs they must agree, and that was checked: on 3,240 spoke sets at
`m = 11, 12, 13` (800 of them random, the rest cost-zero output of the annealer)
the cost values and the number of edge selections agreed line for line, with no
mismatch, and every labeling the C++ engine found for `W_7 … W_13` was
re-verified by the Python `is_geodesic_leech`. Two flags exist for exactly this
kind of check: `--assemble-stdin` reads one spoke set per line and prints its
cost and its number of edge selections, and `--dump-sets` prints every cost-zero
set the annealer emits, so any C++ run can be replayed through the Python
implementation.

The `W_14` figures in Section 9 of the paper come from this engine: eight seeds
(1–8), 18 hours each, 144 CPU-hours in all, on the machine above.

| | per seed | all eight seeds |
|---|---|---|
| annealing steps | ≈ 1.47 × 10¹¹ | 1,180,453,763,200 |
| cost-zero spoke sets | ≈ 7.5 × 10⁵ | 6,027,008 |
| sets that assembled into an admissible frame | 0 | **0** |
| rim completions attempted | 0 | 0 |
| labelings found | 0 | 0 |

Not one of six million cost-zero sets at `m = 13` survived assembly. That is not
a property of the port: a sample of 1,246 of them, dumped with `--dump-sets` and
replayed through the Python `edge_selections`, also assembled 0 of 1,246, and
the failures split into two kinds — in 129 the forced pairs already close a
cycle, and in the other 1,117 the forced paths are fine but the quota edges of
the collision groups cannot all be placed at once. The Hamiltonian cycle has 13
edges; these sets need 2–10 of them for forced pairs and a further 3–11 for
quotas, and the cost bounds the two demands separately. So at `m = 13` the
necessary conditions in the cost are far from sufficient, and the search never
reaches the rim. Assembly is an exact finite problem, and it is the natural
place for a complete method to decide `W_14`.

**This is not the program that produced Table 5.** That code was not kept, and
this is an independent reimplementation written for the repository. The labelings
it finds are valid but different from the published ones — they are far from
unique, and `data/README.md` explains one reason why. `W_13` is a nice
illustration: the search arrived at Table 5's spoke set and Table 5's rim set,
but with the two rotated by different amounts, so the labeling is not an
automorphic image of the published one. Nothing in the paper rests on any of
this: the existence claims rest on the explicit labelings of Table 5, whose
geodesic weight classes Appendix A lists and `code/` checks directly.

## License

Code in `code/` and `search/`: [MIT](LICENSE).
The article in `paper/` and the data in `data/`: [CC BY 4.0](LICENSE-CC-BY-4.0).

## Citing

See [`CITATION.cff`](CITATION.cff), or use GitHub's "Cite this repository" button.
Tagged releases are archived on Zenodo; the concept DOI
[10.5281/zenodo.22254583](https://doi.org/10.5281/zenodo.22254583) always
resolves to the latest release, and each release also has a DOI of its own.
