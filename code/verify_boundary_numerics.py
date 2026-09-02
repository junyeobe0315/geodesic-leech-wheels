#!/usr/bin/env python3
"""Exact audit of every numeric claim behind the bound n <= 40.

The nonexistence proof in the paper is finite and hand-checkable, but it ends in
a few dozen explicit decimal comparisons, some of them tight (the case
(m, h) = (42, 2) has a margin of 0.0072, and (m, h) = (40, 2) a margin of
0.0069). This script re-derives all of them in exact rational interval
arithmetic -- no floating point anywhere -- so that each printed decimal is
certified rather than trusted.

What is checked, section by section:

  Section 7  (uniform Fourier inequalities, Proposition 7.2)
      * the rational bracket 333/106 < pi < 355/113, from Machin's formula;
      * U_0(40) < 39.943 < 40.9523 < L_0(40)                       (h = 0)
      * U_0(41) < 40.4409 < 40.6746 < L_1(41)                      (h = 1)
      * U_2(42) < 40.7355 < 40.7427 < L_2(42)                      (h = 2)
      * the derivative bounds u_m' < 7/32, (v_m)' < 1/5, and
        U' < 0.54 < 1 < L', which push each base case to all larger m.

  Sections 5-6  (occupancy bound and the position formulas)
      * the constant 14(w-1)/(36w-25) that Proposition 5.4 relaxes to 3/8, and the
        closed form (15/8)sqrt(w) + 12/5 that follows;
      * kappa <= 3m and p <= Pi_m for two large spokes (Lemma 6.2), by taking
        the fixed point of the position inequality separately for each d;
      * e <= 3m - 6 for one large spoke (Lemma 6.4);
      * the auxiliary inequalities and the two asymptotic square-root arguments.
    Two of these are tight enough that the relaxation to 3/8 loses them at
    m = 40 and 41, which the script reports as notes.

  Section 8 + Appendix D  (six-variable Parseval lemma)
      * the constants M, R, sigma_M, E_m, U_m^*, omega_m, c_r;
      * the coefficient square-sum bounds 267/10000 and 219/5000;
      * c*_40 < 78.87 and c*_41 < 80.90;
      * the intermediate-value table of Appendix D to the stated 1e-5;
      * the gradient of Psi at each of the three tangent points, hence the
        tangent-plane table (A_0, g_1+g_2+g_3, ||(g_4,g_5,g_6)||_2, Psi);
      * the three lower bounds for Lambda, and Lambda > Psi in each case.

This script does not prove the Fourier lemmas themselves; those are proved in
the paper. It certifies the arithmetic they are combined with.

No third-party package is required.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as Q
from typing import Dict, List, Sequence, Tuple

from interval import (
    PI,
    PI_MINUS,
    PI_PLUS,
    Iv,
    ceil_to_grid,
    sin_lower_S3,
    sin_upper_S5,
)

FAILURES: List[str] = []
NOTES: List[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    status = "ok  " if condition else "FAIL"
    print(f"  [{status}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def note(message: str) -> None:
    """Record a discrepancy that is not a failure of the theorem itself."""
    print(f"  [note] {message}")
    if message not in NOTES:
        NOTES.append(message)


def section(title: str) -> None:
    print(f"\n{title}\n" + "-" * len(title))


# ---------------------------------------------------------------------------
# Part 1: the rational bracket for pi
# ---------------------------------------------------------------------------

def part1_pi() -> None:
    section("Part 1 - rational bounds for pi (Machin's formula)")
    check(
        f"{PI_MINUS.numerator}/{PI_MINUS.denominator} < pi",
        PI.lo > PI_MINUS,
        f"pi > {float(PI.lo):.12f} > {float(PI_MINUS):.12f}",
    )
    check(
        f"pi < {PI_PLUS.numerator}/{PI_PLUS.denominator} = pi_+",
        PI.hi < PI_PLUS,
        f"pi < {float(PI.hi):.12f} < {float(PI_PLUS):.12f}",
    )


# ---------------------------------------------------------------------------
# Part 2: Section 7, the uniform inequalities
# ---------------------------------------------------------------------------

HALF_PI_PLUS = Q(355, 226)  # pi_+ / 2


def u_of(m: int) -> Iv:
    """u_m = sqrt(7m + 3)."""
    return Iv.exact(7 * m + 3).sqrt()


def v0_of(m: int) -> Iv:
    """v_m^(0) = sqrt(6m + 2 + u_m)."""
    return (Iv.exact(6 * m + 2) + u_of(m)).sqrt()


def v2_of(m: int) -> Iv:
    """v_m^(2) = sqrt(6m + 1 + u_m)."""
    return (Iv.exact(6 * m + 1) + u_of(m)).sqrt()


def U0_of(m: int) -> Iv:
    v = v0_of(m)
    return (
        HALF_PI_PLUS * (Iv.exact(6 * m + 2) - v).sqrt()
        + Q(2, 3) * v
        + Q(1, 3) * u_of(m)
    )


def U2_of(m: int) -> Iv:
    v = v2_of(m)
    return (
        HALF_PI_PLUS * (Iv.exact(Q(12 * m - 3, 2)) - v).sqrt()
        + Q(2, 3) * v
        + Q(1, 3) * u_of(m)
    )


def L0_of(m: int) -> Iv:
    return Iv.exact(Q(m + 1) - Q(2 * (m + 1), m * (m + 3)))


def L1_of(m: int) -> Iv:
    return L0_of(m) - Iv.exact(2 * PI_PLUS ** 2 * Q(3 * m - 6, (m + 1) * (m + 2)))


def L2_of(m: int) -> Iv:
    return L0_of(m) - Iv.exact(2 * PI_PLUS ** 2 * Q(5 * m + 2, (m + 1) * (m + 2)))


def shift_polynomial(coeffs: Sequence[int], shift: int) -> List[Q]:
    """Coefficients of p(m0 + t) in t, given p in m (ascending order)."""
    out = [Q(0)] * len(coeffs)
    for degree, coefficient in enumerate(coeffs):
        # expand (shift + t)^degree
        binom = 1
        for k in range(degree + 1):
            out[k] += coefficient * binom * shift ** (degree - k)
            binom = binom * (degree - k) // (k + 1)
    return out


def positive_for_all_m_at_least(coeffs: Sequence[int], m0: int) -> bool:
    """True if p(m) > 0 for every real m >= m0, certified coefficientwise."""
    shifted = shift_polynomial(coeffs, m0)
    return shifted[0] > 0 and all(c >= 0 for c in shifted[1:])


def part2_uniform_ranges(scan: int) -> None:
    section("Part 2 - Section 7, uniform Fourier inequalities (Prop. 7.2)")

    # --- base cases, exact ------------------------------------------------
    cases = [
        ("h=0", 40, U0_of(40), L0_of(40), Q(39943, 1000), Q(409523, 10000)),
        ("h=1", 41, U0_of(41), L1_of(41), Q(404409, 10000), Q(406746, 10000)),
        ("h=2", 42, U2_of(42), L2_of(42), Q(407355, 10000), Q(407427, 10000)),
    ]
    for label, m, upper, lower, printed_u, printed_l in cases:
        check(
            f"{label}: U({m}) < {float(printed_u)}",
            upper.hi < printed_u,
            f"U = {upper.show(7)}",
        )
        check(
            f"{label}: L({m}) > {float(printed_l)}",
            lower.lo > printed_l,
            f"L = {lower.show(7)}",
        )
        check(
            f"{label}: L({m}) - U({m}) > 0",
            lower.lo > upper.hi,
            f"margin >= {float(lower.lo - upper.hi):.7f}",
        )

    # --- the u_m and v_m values printed in the text -----------------------
    printed_intermediates = [
        ("u_40", u_of(40), "16.8226", 4),
        ("v_40^(0)", v0_of(40), "16.0880", 4),
        ("u_41", u_of(41), "17.0293", 4),
        ("v_41^(0)", v0_of(41), "16.2797", 4),
        ("u_42", u_of(42), "17.2336879", 7),
        ("v_42^(2)", v2_of(42), "16.4387860", 7),
    ]
    for name, value, text, decimals in printed_intermediates:
        # The paper writes these as e.g. u_40 = 16.8226..., mostly truncating but
        # occasionally rounding to nearest. Accept either, and say which applies.
        printed = Q(text)
        step = Q(1, 10 ** decimals)
        truncates = printed <= value.lo and value.hi < printed + step
        rounds = printed - step / 2 <= value.lo and value.hi <= printed + step / 2
        style = "truncated" if truncates else "rounded" if rounds else "MISMATCH"
        check(
            f"printed value {name} = {text}...",
            truncates or rounds,
            f"= {value.show(decimals + 3)} ({style})",
        )

    # --- derivative bounds: base case propagates to all larger m ----------
    # u' = 7 / (2 sqrt(7m+3)) is decreasing, so it is largest at m = m0.
    for m0, u_claim, tail_claim in ((40, Q(7, 32), Q(7, 96)), (42, Q(7, 34), Q(7, 102))):
        u_prime = Iv.exact(7) / (2 * u_of(m0))
        check(
            f"u'(m) < {u_claim} for m >= {m0}",
            u_prime.hi < u_claim,
            f"u'({m0}) = {u_prime.show(7)}",
        )
        check(
            f"(1/3) u' < {tail_claim} for m >= {m0}",
            (Q(1, 3) * u_prime).hi < tail_claim,
        )

    # v' = (6 + u') / (2v), with v increasing, so largest at m = m0.
    u40p = Iv.exact(7) / (2 * u_of(40))
    v0p = (6 + u40p) / (2 * v0_of(40))
    check("(v_m^(0))' < 1/5 for m >= 40", v0p.hi < Q(1, 5), f"= {v0p.show(7)}")

    u42p = Iv.exact(7) / (2 * u_of(42))
    v2p = (6 + u42p) / (2 * v2_of(42))
    check("(v_m^(2))' < 1/5 for m >= 42", v2p.hi < Q(1, 5), f"= {v2p.show(7)}")

    # (sqrt(6m+2-v))' = (6 - v') / (2 sqrt(6m+2-v)) <= 6 / (2 sqrt(...)) at m0.
    s0 = (Iv.exact(6 * 40 + 2) - v0_of(40)).sqrt()
    s0p = Iv.exact(6) / (2 * s0)
    check("(sqrt(6m+2-v^(0)))' < 1/5 for m >= 40", s0p.hi < Q(1, 5), f"= {s0p.show(7)}")

    s2 = (Iv.exact(Q(12 * 42 - 3, 2)) - v2_of(42)).sqrt()
    s2p = Iv.exact(6) / (2 * s2)
    check("(sqrt(6m-3/2-v^(2)))' < 1/5 for m >= 42", s2p.hi < Q(1, 5), f"= {s2p.show(7)}")

    U0p = HALF_PI_PLUS * s0p + Q(2, 3) * v0p + Q(1, 3) * u40p
    check("U_0'(m) < 0.54 for m >= 40", U0p.hi < Q(54, 100), f"= {U0p.show(7)}")
    U2p = HALF_PI_PLUS * s2p + Q(2, 3) * v2p + Q(1, 3) * u42p
    check("U_2'(m) < 1 for m >= 42", U2p.hi < 1, f"= {U2p.show(7)}")

    # L' > 1 reduces to the subtracted rational functions being decreasing.
    # d/dm [2(m+1)/(m(m+3))] has numerator -2(m^2+2m+3): need m^2+2m+3 > 0.
    check(
        "L_0'(m) > 1 for m >= 40   (m^2+2m+3 > 0)",
        positive_for_all_m_at_least([3, 2, 1], 40),
    )
    # d/dm [(3m-6)/((m+1)(m+2))] has numerator -3m^2+12m+24: need 3m^2-12m-24 > 0.
    check(
        "L_1'(m) > 1 for m >= 41   (3m^2-12m-24 > 0)",
        positive_for_all_m_at_least([-24, -12, 3], 41),
    )
    # d/dm [(5m+2)/((m+1)(m+2))] has numerator -5m^2-4m+4: need 5m^2+4m-4 > 0.
    check(
        "L_2'(m) > 1 for m >= 42   (5m^2+4m-4 > 0)",
        positive_for_all_m_at_least([-4, 4, 5], 42),
    )

    # --- belt and braces: a direct exact scan over a finite window --------
    bad: List[Tuple[str, int]] = []
    for m in range(40, 40 + scan):
        if not L0_of(m).lo > U0_of(m).hi:
            bad.append(("h=0", m))
    for m in range(41, 41 + scan):
        if not L1_of(m).lo > U0_of(m).hi:
            bad.append(("h=1", m))
    for m in range(42, 42 + scan):
        if not L2_of(m).lo > U2_of(m).hi:
            bad.append(("h=2", m))
    check(
        f"direct exact scan of {scan} consecutive m per case finds no violation",
        not bad,
        f"violations: {bad}" if bad else "",
    )

    print("  => Prop. 7.2: h=0 forces m<=39, h=1 forces m<=40, h=2 forces m<=41.")


# ---------------------------------------------------------------------------
# Part 3: Section 5 occupancy and the Section 6 position formulas
# ---------------------------------------------------------------------------

def occupancy(w: int, sharp: bool = True) -> int:
    """Largest interval occupancy not excluded by the bounds behind Proposition 5.4.

    Three constraints apply to r labels inside an interval of w consecutive
    integers whose non-path pair sums must be pairwise distinct:

      * r <= w;
      * (r-1)(r-2)/2 <= 2w-3, the elementary count of available sums;
      * Lemma 5.3 with L_w(r) >= lambda_w * r.

    ``sharp`` selects which constant is carried into the third one. The proof of
    Proposition 5.4 derives lambda_w = 14(w-1)/(36w-25) and then relaxes it to 3/8,
    which is what produces the printed 73r^2 - 344r + 384 - 256w <= 0. The
    relaxation is harmless in the asymptotic range but loses the boundary cases
    m = 40 and 41, so the sharp form is the default here; part3 checks both and
    reports the difference.
    """
    if w <= 0:
        return 0

    def admissible(r: int) -> bool:
        if r > w or (r - 1) * (r - 2) // 2 > 2 * w - 3:
            return False
        if not sharp:
            return 73 * r * r - 344 * r + 384 - 256 * w <= 0
        # (1 + lam^2) r^2 - (5 + lam) r + 6 - 4w <= 0 with lam = A/B, denominators cleared
        A, B = 14 * (w - 1), 36 * w - 25
        return (B * B + A * A) * r * r - (5 * B + A) * B * r + (6 - 4 * w) * B * B <= 0

    r = 1
    while admissible(r + 1):
        r += 1
    return r


def two_large_bounds(m: int, sharp: bool = True) -> Tuple[int, int, int, int]:
    """Bounds from eq. (7.z-fixedpoint), iterated separately for each d.

    kappa <= floor(3m/2) + rho(d) + 2 rho(kappa) is taken from the a priori bound
    kappa <= 5m-7 down to its fixed point kappa(d). A value of d survives only
    while the accompanying constraint d < kappa/2 stays satisfiable, and that
    self-consistency is what kills the large-d branch: iterating with a single
    worst-case d instead stalls at kappa <= 122 for m = 40.

    Returns (max kappa, the d attaining it, max p, the d attaining it).
    """
    base = (3 * m) // 2
    best_k = best_kd = best_p = best_pd = -1
    for d in range(1, 5 * m):
        kappa = 5 * m - 7
        for _ in range(100):
            nxt = base + occupancy(d, sharp) + 2 * occupancy(kappa, sharp)
            if nxt >= kappa:
                break
            kappa = nxt
        if 2 * d >= kappa:
            continue                     # d < kappa/2 is no longer satisfiable
        if kappa > best_k:
            best_k, best_kd = kappa, d
        if kappa - d - 1 > best_p:
            best_p, best_pd = kappa - d - 1, d
    return best_k, best_kd, best_p, best_pd


def one_large_bound(m: int, sharp: bool = True) -> int:
    """Fixed point of eq. (7.one-large-fixedpoint); returns the bound on kappa."""
    kappa = 5 * m - 7
    for _ in range(100):
        nxt = (3 * m) // 2 + 1 + occupancy(kappa, sharp) + 2 * occupancy(kappa // 2, sharp)
        if nxt >= kappa:
            break
        kappa = nxt
    return kappa


def part3_position_formulas(scan: int) -> None:
    section("Part 3 - Section 5 occupancy and the Section 6 position formulas")

    # --- the constant that the proof of Proposition 5.4 relaxes -----------------
    check(
        "14(w-1)/(36w-25) >= 3/8 for w >= 10, so the printed quadratic follows",
        all(Q(14 * (w - 1), 36 * w - 25) >= Q(3, 8) for w in range(10, 4000))
        and Q(14 * 8, 36 * 9 - 25) < Q(3, 8),
        "and it fails at w = 9, which is why Proposition 5.4 asks for w >= 10",
    )
    check(
        "the sharp constraint is at least as strong as the printed one",
        all(occupancy(w, True) <= occupancy(w, False) for w in range(1, 4000)),
    )
    worst = min(
        Q(41, 64) * w + 12 * Iv.exact(w).sqrt().lo - Q(528, 25) for w in range(10, 4000)
    )
    check(
        "closed form of Proposition 5.4: (41/64)w + 12 sqrt(w) - 528/25 > 0 for w >= 10",
        worst > 0,
        f"least value over 10 <= w < 4000 is {float(worst):.4f}, at w = 10",
    )

    # --- the two-large-spoke position formula -----------------------------
    print("\n  two large spokes: kappa <= floor(3m/2) + rho(d) + 2 rho(kappa), d < kappa/2")
    print(f"    {'m':>4} {'3m':>5} {'kappa<=':>8} {'Pi_m':>6} {'p<=':>5}   verdict")
    failures: List[int] = []
    for m in range(40, 40 + scan):
        Pi = (5 * m) // 2 + 1
        k, kd, p, pd = two_large_bounds(m)
        ok = k <= 3 * m and p <= Pi
        if not ok:
            failures.append(m)
        if m < 44:
            print(
                f"    {m:>4} {3 * m:>5} {k:>8} {Pi:>6} {p:>5}   "
                f"{'OK' if ok else 'FAILS'}   (attained at d = {kd} and d = {pd})"
            )
    check(
        f"eq. (7.z-3m) kappa <= 3m and eq. (7.p-Pm) p <= Pi_m, for 40 <= m < {40 + scan}",
        not failures,
        f"failures at m = {failures}" if failures else "",
    )

    # --- the remark after Proposition 5.4, on why the relaxation is not affordable
    check(
        "near w = 110 the relaxation is the difference between rho <= 21 and rho <= 22",
        occupancy(110, True) == 21 and occupancy(110, False) == 22,
        f"rho(110) = {occupancy(110, True)} from (9), {occupancy(110, False)} from (10)",
    )
    relaxed = {m: two_large_bounds(m, sharp=False) for m in (40, 41)}
    check(
        "the relaxed form (10) would stall at kappa <= 122 for m = 40, as the text says",
        relaxed[40][0] == 122,
        f"and would give only p <= {relaxed[40][2]} there, against Pi_40 = 101; "
        f"at m = 41, kappa <= {relaxed[41][0]} and p <= {relaxed[41][2]} against Pi_41 = 103",
    )

    # --- the one-large-spoke position formula -----------------------------
    print("\n  one large spoke: kappa <= floor(3m/2)+1+rho(kappa)+2rho(floor(kappa/2))")
    failures_one: List[int] = []
    for m in range(40, 40 + scan):
        kappa = one_large_bound(m)
        if kappa > 3 * m - 5:
            failures_one.append(m)
        if m < 44:
            print(f"    m = {m}: kappa <= {kappa} <= 3m-5 = {3 * m - 5}, so e <= {kappa - 1}")
    check(
        f"eq. (7.e-3m): e = kappa - 1 <= 3m - 6, for 40 <= m < {40 + scan}",
        not failures_one,
        f"failures at m = {failures_one}" if failures_one else "",
    )

    # --- the auxiliary inequalities quoted for m >= 42 --------------------
    aux = [
        m
        for m in range(42, 42 + scan)
        if not (
            2 * occupancy(3 * m) <= m + 4
            and 2 * occupancy((5 * m) // 2 + 11) <= m + 2
        )
    ]
    check(
        "eq. (7.rho-two-aux): 2rho(3m) <= m+4 and 2rho(Pi_m+10) <= m+2 for m >= 42",
        not aux,
        f"failures at m = {aux}" if aux else "",
    )

    # --- the asymptotic square-root arguments -----------------------------
    check(
        "1/sqrt(2) < 71/100 and (15/8)(2 + 71/100) = 813/160",
        2 * 71 ** 2 > 100 ** 2 and Q(15, 8) * (2 + Q(71, 100)) == Q(813, 160),
    )
    check(
        "two-spoke discriminant = (57600m^2 - 2459067m + 323095)/25600 > 0 for m >= 43",
        all(
            (Q(3 * m, 2) - Q(31, 5)) ** 2 - Q(813, 160) ** 2 * (3 * m + 1)
            == Q(57600 * m * m - 2459067 * m + 323095, 25600)
            for m in (43, 60, 200)
        )
        and positive_for_all_m_at_least([323095, -2459067, 57600], 43),
    )
    check(
        "1/sqrt(2) < 99/140 and (15/8)(1 + 2*99/140) = 507/112",
        2 * 99 ** 2 > 140 ** 2 and Q(15, 8) * (1 + 2 * Q(99, 140)) == Q(507, 112),
    )
    check(
        "507/112 genuinely needs the sharper 99/140, not the 71/100 used two-spoke-side",
        Q(15, 8) * (1 + 2 * Q(71, 100)) > Q(507, 112),
        f"71/100 would give {float(Q(15, 8) * (1 + 2 * Q(71, 100))):.6f} "
        f"> {float(Q(507, 112)):.6f}",
    )
    check(
        "one-spoke discriminant = (705600m^2 - 30756435m + 72381124)/313600 > 0 for m >= 42",
        all(
            (Q(3 * m, 2) - Q(61, 5)) ** 2 - Q(507, 112) ** 2 * (3 * m - 4)
            == Q(705600 * m * m - 30756435 * m + 72381124, 313600)
            for m in (42, 60, 200)
        )
        and positive_for_all_m_at_least([72381124, -30756435, 705600], 42),
    )

    print("\n  => Lemma 6.2: kappa <= 3m and p <= Pi_m.   Lemma 6.4: e <= 3m - 6.")


# ---------------------------------------------------------------------------
# Part 4: Section 8 and Appendix D, the six-variable Parseval lemma
# ---------------------------------------------------------------------------

KAPPA_X = Q(267, 10000)
KAPPA_Y = Q(219, 5000)
C2 = Q(2, 15)   # c_2
C3 = Q(2, 35)   # c_3
C4 = Q(2, 63)   # c_4
C6 = Q(2, 143)  # c_6


def c_of(r: int) -> Q:
    return Q(2, 4 * r * r - 1)


@dataclass(frozen=True)
class BoundaryCase:
    m: int
    h: int
    x0: Tuple[Q, ...]

    @property
    def label(self) -> str:
        return f"({self.m},{self.h})"

    @property
    def M(self) -> int:
        return (self.m + 1) * (self.m + 2) // 2

    @property
    def R(self) -> int:
        return (self.M - 1) // 2

    @property
    def E(self) -> int:
        return self.m * (self.M - 2 * self.m)

    @property
    def omega(self) -> Q:
        return Q(1, 7) - Q(1, self.M)

    @property
    def sigma(self) -> int:
        for s in range(1, self.R + 1):
            if (4 * s) % self.M in (1, self.M - 1):
                return s
        raise ValueError("no sigma_M found")

    @property
    def D1_D2(self) -> Tuple[Q, Q]:
        if self.h == 1:
            return Q(self.m + 1), Q(self.m + 1)
        return Q(self.m) - Q(3, 4), Q(self.m) + Q(1, 2)

    @property
    def Ustar(self) -> Iv:
        return (1 + Iv.exact(24 * self.m + 9).sqrt()) / 2


def rational_points(values: Sequence[str]) -> Tuple[Q, ...]:
    return tuple(Q(v) for v in values)


CASES = (
    BoundaryCase(40, 2, rational_points(("78.86", "78.86", "78.86", "75.64", "60.64", "27.20"))),
    BoundaryCase(40, 1, rational_points(("78.86", "78.86", "78.86", "75.61", "60.66", "27.24"))),
    BoundaryCase(41, 2, rational_points(("80.89", "80.89", "80.89", "79.75", "64.12", "28.80"))),
)

# Table of Section 8: A_0 <, g1+g2+g3 <, ||(g4,g5,g6)||_2 <, Psi <
TANGENT_TABLE: Dict[str, Tuple[Q, Q, Q, Q]] = {
    "(40,2)": (Q("27.444"), Q("0.14247"), Q("0.00000359"), Q("38.682")),
    "(40,1)": (Q("27.717"), Q("0.14159"), Q("0.00000503"), Q("38.886")),
    "(41,2)": (Q("27.782"), Q("0.14107"), Q("0.00000399"), Q("39.196")),
}

# Appendix D intermediate values, stated to five decimals within 1e-5.
APPENDIX_C: Dict[str, Dict[str, Q]] = {
    "(40,2)": {
        "eta": Q("49.44534"), "X": Q("25.58936"), "Y": Q("15.29318"),
        "w4": Q("15.99326"), "v2": Q("15.95974"), "sqrtT": Q("14.84117"),
        "W1": Q("15.95401"), "W2": Q("15.67251"), "JW": Q("11.34243"),
        "W": Q("12.02170"), "Psi": Q("38.67847"),
    },
    "(40,1)": {
        "eta": Q("49.44468"), "X": Q("25.59767"), "Y": Q("15.29297"),
        "w4": Q("15.99326"), "v2": Q("15.99104"), "sqrtT": Q("14.95757"),
        "W1": Q("15.95405"), "W2": Q("15.67070"), "JW": Q("11.34241"),
        "W": Q("12.02239"), "Psi": Q("38.88219"),
    },
    "(41,2)": {
        "eta": Q("52.25463"), "X": Q("27.04591"), "Y": Q("16.15482"),
        "w4": Q("16.18752"), "v2": Q("16.15449"), "sqrtT": Q("15.03747"),
        "W1": Q("16.15228"), "W2": Q("16.00200"), "JW": Q("11.51371"),
        "W": Q("12.23133"), "Psi": Q("39.19229"),
    },
}

# Lower bounds for Lambda printed in Section 8, and the multiplier of pi_+.
LAMBDA_ROWS: Dict[str, Tuple[int, int, int, int, Q]] = {
    #        numerator/denominator of (M-1)/M * (m+1),  multiple, argument, printed
    "(40,2)": (860, 861, 41, 101, Q("38.6889")),
    "(40,1)": (860, 861, 41, 114, Q("39.6829")),
    "(41,2)": (902, 903, 42, 103, Q("39.7498")),
}
LAMBDA_PI_MULTIPLE = {"(40,2)": 2, "(40,1)": 1, "(41,2)": 2}


def sum_of_squares_upper(indices, terms) -> Q:
    """Upper bound for sum of c_r^2 over the given indices, on the output grid."""
    total = Q(0)
    for r in indices:
        total = ceil_to_grid(total + terms(r) ** 2)
    return total


def c_star_upper(m: int, M: int, pi_lower: Q | None = None) -> Iv:
    """Upper bound for c*_m = sin(2 m pi / M) / sin(pi / M).

    Both arguments lie well inside (0, pi/2), where sin is increasing, so an
    upper bound follows from S_5 at the largest possible numerator argument
    and S_3 at the smallest possible denominator argument.

    c*_m is of size M/pi, so a relative error in pi is amplified by roughly M.
    Passing ``pi_lower`` reproduces what a given rational lower bound for pi
    yields; the default uses the Machin enclosure, which is exact for this
    purpose.
    """
    assert 2 * Q(m) * PI_PLUS / M < PI.lo / 2, "2 m pi / M must stay below pi/2"
    x_hi = Q(2 * m) * PI.hi / M
    x_lo = (PI.lo if pi_lower is None else pi_lower) / M
    return Iv.exact(sin_upper_S5(x_hi)) / Iv.exact(sin_lower_S3(x_lo))


def evaluate_case(case: BoundaryCase) -> Dict[str, object]:
    """Evaluate Psi and its gradient at the tangent point, in interval arithmetic."""
    m, M, R, E = case.m, case.M, case.R, case.E
    omega, sigma = case.omega, case.sigma
    D1, D2 = case.D1_D2
    Ustar = case.Ustar
    q1, q2, q4, q6, q8, q12 = (Iv.exact(v) for v in case.x0)

    eta = (Iv.exact(E) - sum((Iv.exact(v) ** 2 for v in case.x0), Iv.exact(0))).sqrt()

    def two_step(inner: Iv) -> Iv:
        return (2 * inner + (2 * m + 2) + Ustar).sqrt()

    def V(X: Iv, Y: Iv) -> Tuple[Iv, Iv]:
        J = two_step(Y)
        return (2 * X + (2 * m + 2) + J).sqrt(), J

    w4 = two_step(q4)
    v2 = (2 * q2 + 2 * D2 + w4).sqrt()
    T = 2 * q1 + 2 * D1 - v2
    sqrtT = T.sqrt()

    W1, J1 = V(q4, q8)
    W2, J2 = V(q6, q12)

    cR, cRm1, cSig = c_of(R), c_of(R - 1), c_of(sigma)
    aX = {"q1": cR, "q2": Q(0), "q4": Q(0), "q6": Q(0), "q8": C4, "q12": C6}
    aY = {"q1": cSig, "q2": cR, "q4": Q(0), "q6": cRm1, "q8": Q(0), "q12": Q(0)}

    X = (cR * q1 + C4 * q8 + C6 * q12 + KAPPA_X * eta) / omega
    Y = (cSig * q1 + cR * q2 + cRm1 * q6 + KAPPA_Y * eta) / omega
    W, JW = V(X, Y)

    Psi = HALF_PI_PLUS * sqrtT + Q(2, 3) * v2 + C2 * W1 + C3 * W2 + omega * W

    coordinates = {"q1": q1, "q2": q2, "q4": q4, "q6": q6, "q8": q8, "q12": q12}

    def tail(name: str) -> Iv:
        dX = (Iv.exact(aX[name]) - KAPPA_X * coordinates[name] / eta) / omega
        dY = (Iv.exact(aY[name]) - KAPPA_Y * coordinates[name] / eta) / omega
        return omega * (dX / W + dY / (2 * W * JW))

    common = Q(2, 3) - HALF_PI_PLUS / (2 * sqrtT)
    gradient = [
        HALF_PI_PLUS / sqrtT + tail("q1"),
        common / v2 + tail("q2"),
        common / (2 * v2 * w4) + C2 / W1 + tail("q4"),
        C3 / W2 + tail("q6"),
        C2 / (2 * W1 * J1) + tail("q8"),
        C3 / (2 * W2 * J2) + tail("q12"),
    ]

    dot = sum((g * Iv.exact(x) for g, x in zip(gradient, case.x0)), Iv.exact(0))
    A0 = Psi - dot

    return {
        "eta": eta, "X": X, "Y": Y, "w4": w4, "v2": v2, "sqrtT": sqrtT,
        "W1": W1, "W2": W2, "JW": JW, "W": W, "Psi": Psi,
        "gradient": gradient, "A0": A0, "sigma": sigma,
    }


def lambda_lower(case: BoundaryCase) -> Q:
    """Exact rational lower bound for Lambda from the second table of Section 8."""
    numerator, denominator, factor, argument, _ = LAMBDA_ROWS[case.label]
    multiple = LAMBDA_PI_MULTIPLE[case.label]
    u = Q(argument) * PI_PLUS / denominator
    return Q(numerator, denominator) * factor - multiple * PI_PLUS * sin_upper_S5(u)


def part4_boundary_lemma() -> None:
    section("Part 4 - Section 8 + Appendix D, six-variable Parseval lemma")

    for case in CASES:
        print(f"\n  case (m,h) = {case.label}")
        m, M, R = case.m, case.M, case.R

        # --- constants ----------------------------------------------------
        check(f"    M = {M}, R = {R}, 2R+1 = M", 2 * R + 1 == M)
        check(f"    E_{m} = m(M-2m) = {case.E}", case.E == m * (M - 2 * m))
        check(
            f"    sigma_M = {case.sigma}, 4 sigma = +-1 mod M",
            (4 * case.sigma) % M in (1, M - 1),
        )
        telescoped = sum(c_of(r) for r in range(4, R + 1))
        check(f"    omega_m = 1/7 - 1/M = {case.omega}", telescoped == case.omega)

        # --- coefficient square sums --------------------------------------
        I_X = [5] + list(range(7, R))
        I_Y = [r for r in range(4, R + 1) if r not in (case.sigma, R - 1, R)]
        sX = sum_of_squares_upper(I_X, c_of)
        sY = sum_of_squares_upper(I_Y, c_of)
        check(
            "    (sum_{I_X} c_r^2)^{1/2} < 267/10000",
            sX < KAPPA_X ** 2,
            f"sqrt = {float(sX) ** 0.5:.8f}",
        )
        check(
            "    (sum_{I_Y} c_r^2)^{1/2} < 219/5000",
            sY < KAPPA_Y ** 2,
            f"sqrt = {float(sY) ** 0.5:.8f}",
        )

        # --- c*_m ----------------------------------------------------------
        printed_cstar = Q("78.87") if m == 40 else Q("80.90")
        cstar = c_star_upper(m, M)
        check(
            f"    c*_{m} < {float(printed_cstar)}",
            cstar.hi < printed_cstar,
            f"c* = {cstar.show(6)}",
        )
        # The paper attributes this bound to "Taylor inequalities and
        # 103993/33102 < pi < 355/113", so check that this bracket -- and not
        # just an exact pi -- really delivers it. c*_m is of size M/pi, so it
        # amplifies a relative error in pi by a factor of about M; the earlier
        # convergent 333/106, off by 8.3e-5, leaves c*_40 < 78.8711 and misses.
        stated = c_star_upper(m, M, pi_lower=PI_MINUS)
        check(
            f"    ... and the paper's own bracket for pi suffices",
            stated.hi < printed_cstar,
            f"gives {stated.show(6)}",
        )
        if not stated.hi < printed_cstar:
            note(
                f"c*_{m} < {float(printed_cstar)} is true ({cstar.show(6)}), but the "
                f"paper's stated lower bound pi > {PI_MINUS.numerator}/"
                f"{PI_MINUS.denominator} only gives {stated.show(6)}."
            )
        check(
            "    tangent point respects q_1,q_2,q_4 <= c*_m",
            all(x <= printed_cstar for x in case.x0[:3]),
        )

        # --- phase-pair hypotheses in the worst case ----------------------
        D1, D2 = case.D1_D2
        beta = Q(2, 3)
        A_worst = Iv.exact(2 * D1)  # q_1 = 0
        w4_worst = (2 * Iv.exact(printed_cstar) + (2 * m + 2) + case.Ustar).sqrt()
        V_worst = (2 * Iv.exact(printed_cstar) + 2 * D2 + w4_worst).sqrt()
        tau = HALF_PI_PLUS / (2 * (A_worst - V_worst).sqrt())
        check("    tau < 0.1 < beta", tau.hi < Q(1, 10) < beta, f"tau = {tau.show(6)}")
        lhs = A_worst * (beta - tau)
        rhs = V_worst * (2 * beta - tau)
        check(
            "    A(beta-tau) > 44 > 21 > V(2beta-tau)",
            lhs.lo > 44 and 21 > rhs.hi,
            f"{lhs.show(4)} > 44 > 21 > {rhs.show(4)}",
        )

        # --- Appendix D intermediate values -------------------------------
        values = evaluate_case(case)
        tolerance = Q(1, 10 ** 5)
        stale = []
        for name, printed in APPENDIX_C[case.label].items():
            iv = values[name]
            if not (printed - tolerance <= iv.lo and iv.hi <= printed + tolerance):
                stale.append((name, printed, iv.show(6)))
        check(
            "    Appendix D intermediate values agree to 1e-5",
            not stale,
            f"off: {stale}" if stale else "",
        )

        # --- the tangent-plane table --------------------------------------
        gradient = values["gradient"]
        A0 = values["A0"]
        head = gradient[0] + gradient[1] + gradient[2]
        tail_norm = sum((g ** 2 for g in gradient[3:]), Iv.exact(0)).sqrt()
        printed_A0, printed_head, printed_tail, printed_psi = TANGENT_TABLE[case.label]

        check(
            "    g_1, g_2, g_3 > 0 (needed for the c*_m substitution)",
            all(g.lo > 0 for g in gradient[:3]),
        )
        check(f"    A_0 < {float(printed_A0)}", A0.hi < printed_A0, f"A_0 = {A0.show(6)}")
        check(
            f"    g_1+g_2+g_3 < {float(printed_head)}",
            head.hi < printed_head,
            f"= {head.show(8)}",
        )
        check(
            f"    ||(g_4,g_5,g_6)||_2 < {float(printed_tail)}",
            tail_norm.hi < printed_tail,
            f"= {float(tail_norm.hi):.3e}",
        )

        root_E = Iv.exact(case.E).sqrt()
        printed_root = Q(177) if case.E == 31240 else Q(184)
        check(
            f"    sqrt(E_{m}) = sqrt({case.E}) < {printed_root}",
            root_E.hi < printed_root,
            f"= {root_E.show(4)}",
        )

        # Psi <= A_0 + c*(g1+g2+g3) + sqrt(E) ||(g4,g5,g6)||, using the paper's
        # own rational substitutes for c* and sqrt(E).
        bound = (
            Iv.exact(printed_A0)
            + Iv.exact(printed_cstar) * Iv.exact(printed_head)
            + Iv.exact(printed_root) * Iv.exact(printed_tail)
        )
        check(
            f"    tangent-plane bound Psi < {float(printed_psi)}",
            bound.hi < printed_psi,
            f"= {bound.show(6)}",
        )

        # --- the lower bound for Lambda, and the contradiction ------------
        lower = lambda_lower(case)
        printed_lambda = LAMBDA_ROWS[case.label][4]
        check(
            f"    Lambda > {float(printed_lambda)}",
            lower > printed_lambda,
            f"= {float(lower):.7f}",
        )
        check(
            "    contradiction: lower bound for Lambda exceeds upper bound for Psi",
            lower > printed_psi,
            f"{float(lower):.6f} > {float(printed_psi)}  (margin {float(lower - printed_psi):.6f})",
        )

    print("\n  => Proposition 8.2: (40,1), (40,2), (41,2) are all impossible.")
    print("  => With Prop. 7.2: every m >= 40 is impossible, so n = m+1 <= 40.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scan",
        type=int,
        default=200,
        help="how many consecutive m to check directly per case (default: 200)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    print(__doc__.strip().splitlines()[0])
    part1_pi()
    part2_uniform_ranges(args.scan)
    part3_position_formulas(args.scan)
    part4_boundary_lemma()

    print()
    if FAILURES:
        raise SystemExit(f"FAILED {len(FAILURES)} check(s): {FAILURES}")
    print("PASS: every numeric claim behind Theorem 1.1 (n <= 40) is certified exactly.")
    if NOTES:
        print(
            f"\n{len(NOTES)} note(s) on the paper's stated justifications "
            "(the claims themselves hold):"
        )
        for message in NOTES:
            print(f"  * {message}")


if __name__ == "__main__":
    main()
