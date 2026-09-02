# Data

## `wheel_certificates_W5_W13.json`

The explicit geodesic Leech labelings of `W_5, ..., W_13` printed in Table 5 of
the paper, in machine-readable form. These are the certificates behind
Theorem 9.1 and the left-hand inclusion of Corollary 1.2.

`W_7, ..., W_13` are new to the paper. `W_5` and `W_6` are the two labelings of
Lakshmanan S. and Manattu ([arXiv:2502.16628](https://arxiv.org/abs/2502.16628),
Figure 3), transcribed here into the same format so that every element of the
left-hand inclusion is certified from one place rather than by reference.

```json
{
  "certificates": [
    { "n": 7, "m": 6, "t_gp": 27,
      "spokes": [4, 5, 2, 14, 22, 7],
      "rims":   [15, 8, 17, 3, 10, 1] }
  ]
}
```

| field | meaning |
|---|---|
| `n` | the wheel `W_n` has `n` vertices **in total**: a hub joined to a cycle of length `m = n - 1`. This is the convention of Lakshmanan S. and Manattu, not the one in which `W_n` denotes a hub plus an `n`-cycle. |
| `m` | `n - 1`, the rim length |
| `t_gp` | the geodesic path number `m(m+3)/2`; the geodesic weights must be exactly `1..t_gp` |
| `spokes` | `spokes[i]` labels the edge from the hub to rim vertex `v_i` |
| `rims` | `rims[i]` labels the rim edge `v_i v_{i+1}`, indices mod `m` |

Both lists are in cyclic order along the rim and have length `m`. A labeling is
determined only up to rotation and reflection of the rim — the `2m` images of the
wheel's automorphism group, which move `spokes` and `rims` together.

There is more slack than that, though, and it is worth knowing before comparing
two labelings. In the frame-completion criterion (Proposition 3.2) the spokes and the rims
never interact through their *relative* alignment: the spoke side contributes
`{a_i}` and the sums over rim-nonadjacent pairs, the rim side contributes `{b_i}`
and `{b_i + b_{i+1}}`, and the two only have to stay disjoint as sets. So
rotating `rims` while leaving `spokes` fixed turns one geodesic Leech labeling
into `m` of them, none related to the others by an automorphism.

To check a certificate yourself, without using anything in this repository: build
the weighted wheel, enumerate every shortest path between every pair of vertices,
and confirm the multiset of path weights is exactly `{1, ..., t_gp}`. That is what
`code/verify_wheel_certificates.py` does, from the graph up, with no wheel-specific
shortcut in the primary check.
