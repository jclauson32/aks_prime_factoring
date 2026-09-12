"""Round 33: continued fraction factorisation -- relations without a sieve.

The convergents of ``sqrt(N)`` give
``A_i^2 = (-1)^(i+1) Q_(i+1) (mod N)`` with ``|Q| < 2 sqrt(N)`` -- residues as
small as the quadratic sieve's, for free, one per continued-fraction step.
Collect the smooth ones and finish by linear algebra mod 2, exactly as the
sieve does.

CFRAC and the quadratic sieve share a mechanism, a relation size and an
asymptotic (`L[1/2]`).  What separates them is that the sieve *sieves*: it finds
its smooth values in bulk, while CFRAC must trial-divide every residue it
generates.
"""

from __future__ import annotations

from math import gcd, isqrt

from .arith import sieve
from .qs import _combine


def cf_residues(n: int, count: int):
    """Yield ``(A_i mod n, Q_(i+1), sign)`` from the continued fraction of ``sqrt n``."""
    a0 = isqrt(n)
    m, d, a = 0, 1, a0
    num_prev, num = 1, a0 % n
    for i in range(count):
        m = d * a - m
        d = (n - m * m) // d
        if d == 0:
            return
        a = (a0 + m) // d
        yield num, d, (-1) ** (i + 1)
        num_prev, num = num, (a * num + num_prev) % n


def factor_base(n: int, bound: int):
    """Primes ``p <= bound`` with ``N`` a square mod ``p`` (plus 2)."""
    return [2] + [p for p in sieve(bound)[1:] if pow(n, (p - 1) // 2, p) == 1]


def cfrac(n: int, bound: int | None = None, max_steps: int | None = None,
          stats: dict | None = None):
    """A non-trivial factor of the odd composite ``n``, or ``None``."""
    if n % 2 == 0:
        return 2
    r = isqrt(n)
    if r * r == n:
        return r
    from math import exp, log, sqrt

    if bound is None:
        ln = log(n)
        bound = int(exp(0.5 * sqrt(ln * log(ln)))) + 30
    base = factor_base(n, bound)
    spare = 20                   # each dependency splits with probability about 1/2
    steps = max_steps or 4000 * (len(base) + spare)
    relations, tried = [], 0
    stream = cf_residues(n, steps)
    for attempt in range(4):
        need = len(base) + spare * (attempt + 1)
        for a_i, q, sign in stream:
            tried += 1
            v, e = q, [1 if sign < 0 else 0]
            for p in base:
                k = 0
                while v % p == 0:
                    v //= p
                    k += 1
                e.append(k)
            if v == 1:
                relations.append((a_i, e))
                if len(relations) >= need:
                    break
        if stats is not None:
            stats.update(base=len(base), relations=len(relations), residues=tried,
                         attempts=attempt + 1)
        if len(relations) <= len(base):
            return None
        got = _combine(n, base, relations)
        if got:
            return got
    return None
