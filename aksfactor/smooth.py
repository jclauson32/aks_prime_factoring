"""Round 31: the smoothness wall, measured against Dickman's function.

Every sub-exponential factoring method is priced by the density of smooth
numbers.  This module computes Dickman's ``rho`` numerically and measures the
three smoothness probabilities the mechanisms actually depend on:

* random integers of a given size -- the baseline, ``rho(u)``;
* ``p - 1`` for random primes ``p`` -- what Pollard's method needs;
* ``#E(F_p)`` over random curves -- what ECM needs, one fresh draw per curve.

The last two are *not* the baseline: shifted primes are even and small-factor
rich, and elliptic orders fill the Hasse interval with a bias of their own.
"""

from __future__ import annotations

from math import exp, log


def dickman_rho(u: float, steps_per_unit: int = 400) -> float:
    """``rho(u)``: the density of ``x^(1/u)``-smooth numbers near ``x``.

    ``rho = 1`` on ``[0, 1]`` and ``rho(u) = 1 - int_1^u rho(t-1)/t dt``,
    integrated with the trapezoid rule on a fixed grid.
    """
    if u <= 0:
        return 0.0
    if u <= 1:
        return 1.0
    n = int(u * steps_per_unit) + 1
    h = u / n
    values = [1.0]                      # values[i] = rho(i * h)
    for i in range(1, n + 1):
        t = i * h
        if t <= 1:
            values.append(1.0)
            continue
        total = 0.0
        j = int(1 / h)
        prev_t = j * h
        prev = values[max(0, j - int(1 / h))] / prev_t if prev_t else 0.0
        while j < i:
            cur_t = (j + 1) * h
            idx = j + 1 - int(1 / h)
            cur = values[max(0, idx)] / cur_t
            total += (prev + cur) * h / 2
            prev, prev_t = cur, cur_t
            j += 1
        values.append(max(0.0, 1.0 - total))
    return values[n]


def is_smooth(m: int, bound: int, primes=None) -> bool:
    from .arith import sieve

    for p in (primes if primes is not None else sieve(bound)):
        while m % p == 0:
            m //= p
        if m == 1:
            return True
    return m == 1


def smooth_rate(sample, bound: int) -> float:
    from .arith import sieve

    primes = sieve(bound)
    return sum(is_smooth(m, bound, primes) for m in sample) / len(sample)


def curve_order(a: int, b: int, p: int) -> int:
    """``#E(F_p)`` for ``y^2 = x^3 + a x + b`` by counting (small ``p``)."""
    total = 1
    for x in range(p):
        r = (x ** 3 + a * x + b) % p
        total += 1 if r == 0 else (2 if pow(r, (p - 1) // 2, p) == 1 else 0)
    return total


def theory(size: int, bound: int) -> float:
    """``rho(log size / log bound)``."""
    return dickman_rho(log(size) / log(bound))
