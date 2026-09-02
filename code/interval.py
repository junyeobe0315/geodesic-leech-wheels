#!/usr/bin/env python3
"""Rigorous rational interval arithmetic.

Every quantity is a closed interval with exact ``Fraction`` endpoints, and every
operation rounds outward, so a computed interval always contains the true value.
This lets the numeric claims of the paper be checked without floating point: a
statement such as "U_0(40) < 39.943" is verified by evaluating U_0(40) as an
interval and comparing its upper endpoint against the exact rational 39943/1000.

To keep denominators from growing without bound, endpoints are snapped to a
fixed grid after every operation -- the lower endpoint downward and the upper
endpoint upward, so the enclosure property is preserved.

No third-party package is required.
"""
from __future__ import annotations

from fractions import Fraction as Q
from math import isqrt

# Grid spacing for outward rounding: 10^-40 is far finer than any comparison in
# the paper, whose tightest margin is of order 10^-5.
GRID = 10 ** 40


def floor_to_grid(x: Q) -> Q:
    return Q(x.numerator * GRID // x.denominator, GRID)


def ceil_to_grid(x: Q) -> Q:
    return Q(-((-x.numerator * GRID) // x.denominator), GRID)


def sqrt_lower(x: Q) -> Q:
    """Largest grid rational r with r^2 <= x."""
    if x < 0:
        raise ValueError("sqrt of a negative lower bound")
    k = isqrt(x.numerator * GRID * GRID // x.denominator)
    r = Q(k, GRID)
    assert r * r <= x
    return r


def sqrt_upper(x: Q) -> Q:
    """Smallest grid rational r with r^2 >= x."""
    r = sqrt_lower(x) + Q(1, GRID)
    assert r * r >= x
    return r


class Iv:
    """A closed interval [lo, hi] of rationals."""

    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None) -> None:
        lo = Q(lo)
        hi = lo if hi is None else Q(hi)
        if lo > hi:
            raise ValueError(f"empty interval [{lo}, {hi}]")
        self.lo = floor_to_grid(lo)
        self.hi = ceil_to_grid(hi)

    # -- construction ------------------------------------------------------
    @staticmethod
    def exact(value) -> "Iv":
        """An interval that is a single exact rational (no widening)."""
        q = Q(value)
        out = Iv.__new__(Iv)
        out.lo = q
        out.hi = q
        return out

    @staticmethod
    def _raw(lo: Q, hi: Q) -> "Iv":
        out = Iv.__new__(Iv)
        out.lo = floor_to_grid(lo)
        out.hi = ceil_to_grid(hi)
        return out

    @staticmethod
    def _coerce(other) -> "Iv":
        return other if isinstance(other, Iv) else Iv.exact(other)

    # -- arithmetic --------------------------------------------------------
    def __add__(self, other) -> "Iv":
        o = Iv._coerce(other)
        return Iv._raw(self.lo + o.lo, self.hi + o.hi)

    __radd__ = __add__

    def __neg__(self) -> "Iv":
        return Iv._raw(-self.hi, -self.lo)

    def __sub__(self, other) -> "Iv":
        return self + (-Iv._coerce(other))

    def __rsub__(self, other) -> "Iv":
        return Iv._coerce(other) + (-self)

    def __mul__(self, other) -> "Iv":
        o = Iv._coerce(other)
        products = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return Iv._raw(min(products), max(products))

    __rmul__ = __mul__

    def __truediv__(self, other) -> "Iv":
        o = Iv._coerce(other)
        if o.lo <= 0 <= o.hi:
            raise ZeroDivisionError(f"divisor interval [{o.lo}, {o.hi}] straddles 0")
        return self * Iv._raw(1 / o.hi, 1 / o.lo)

    def __rtruediv__(self, other) -> "Iv":
        return Iv._coerce(other) / self

    def __pow__(self, n: int) -> "Iv":
        if n < 0:
            raise ValueError("negative powers are not needed here")
        result = Iv.exact(1)
        for _ in range(n):
            result = result * self
        return result

    def sqrt(self) -> "Iv":
        if self.lo < 0:
            raise ValueError(f"sqrt of interval [{self.lo}, {self.hi}] with negative lower end")
        return Iv._raw(sqrt_lower(self.lo), sqrt_upper(self.hi))

    # -- comparison --------------------------------------------------------
    def __lt__(self, other) -> bool:
        """True only if the whole interval is strictly below the whole of other."""
        return self.hi < Iv._coerce(other).lo

    def __gt__(self, other) -> bool:
        return self.lo > Iv._coerce(other).hi

    def contains_zero(self) -> bool:
        return self.lo <= 0 <= self.hi

    @property
    def width(self) -> Q:
        return self.hi - self.lo

    # -- display -----------------------------------------------------------
    def __repr__(self) -> str:
        return f"[{float(self.lo):.10f}, {float(self.hi):.10f}]"

    def show(self, digits: int = 6) -> str:
        return f"{float(self.lo):.{digits}f}"


def arctan_bounds(x: Q, terms: int) -> Iv:
    """Enclose arctan(x) for 0 < x <= 1 by an alternating partial sum.

    The series sum (-1)^k x^(2k+1)/(2k+1) is alternating with strictly
    decreasing terms, so the truncation error is bounded by the first omitted
    term.
    """
    if not 0 < x <= 1:
        raise ValueError("arctan_bounds expects 0 < x <= 1")
    total = Q(0)
    power = x
    x2 = x * x
    for k in range(terms):
        total += (-1) ** k * power / (2 * k + 1)
        power *= x2
    tail = power / (2 * terms + 1)
    return Iv(total - tail, total + tail)


def pi_bounds() -> Iv:
    """Enclose pi via Machin's formula pi/4 = 4 arctan(1/5) - arctan(1/239)."""
    a5 = arctan_bounds(Q(1, 5), 40)
    a239 = arctan_bounds(Q(1, 239), 12)
    return 4 * (4 * a5 - a239)


PI = pi_bounds()
PI_PLUS = Q(355, 113)         # the paper's rational upper bound for pi
PI_MINUS = Q(103993, 33102)   # the paper's rational lower bound for pi


def sin_upper_S5(u: Q) -> Q:
    """S_5(u) = u - u^3/6 + u^5/120 >= sin(u) for 0 <= u <= 1 (exact rational)."""
    if not 0 <= u <= 1:
        raise ValueError("S_5 is used by the paper only on [0, 1]")
    return u - u ** 3 / 6 + u ** 5 / 120


def sin_lower_S3(u: Q) -> Q:
    """S_3(u) = u - u^3/6 <= sin(u) for u >= 0 (exact rational)."""
    if u < 0:
        raise ValueError("S_3 is used here only for u >= 0")
    return u - u ** 3 / 6
