"""Elementary number theory helpers.

Pure standard library.  Nothing here is specific to the Pascal/AKS results; it
is the arithmetic substrate the rest of the package is built on.
"""

from __future__ import annotations

import random
from math import gcd, isqrt

__all__ = [
    "is_prime",
    "sieve",
    "factorize_small",
    "spf_trial",
    "pollard_rho",
    "factorize",
    "valuation",
    "carries",
]

# Deterministic Miller-Rabin bases: this set is a proven witness set for all
# n < 3_317_044_064_679_887_385_961_981.  Above that bound the test is a
# (very strong) probable-prime test, which we flag in the docstring rather
# than pretend otherwise.
_MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
_MR_DETERMINISTIC_BELOW = 3_317_044_064_679_887_385_961_981


def is_prime(n: int) -> bool:
    """Miller-Rabin primality test.

    Deterministic for ``n < 3.3e24``; a probable-prime test above that bound.
    """
    if n < 2:
        return False
    for p in _MR_BASES:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in _MR_BASES:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def sieve(limit: int) -> list[int]:
    """All primes <= ``limit`` by a plain sieve of Eratosthenes."""
    if limit < 2:
        return []
    flags = bytearray([1]) * (limit + 1)
    flags[0] = flags[1] = 0
    for i in range(2, isqrt(limit) + 1):
        if flags[i]:
            flags[i * i :: i] = bytearray(len(flags[i * i :: i]))
    return [i for i, f in enumerate(flags) if f]


def valuation(n: int, p: int) -> int:
    """The exponent of ``p`` in ``n`` (with ``valuation(0, p) == 0``)."""
    if n == 0:
        return 0
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def carries(a: int, b: int, p: int) -> list[int]:
    """Carry-out indicators, least significant digit first, for ``a + b`` base ``p``.

    By Kummer's theorem the number of carries when adding ``M`` and ``N - M``
    in base ``p`` equals the exponent of ``p`` in ``C(N, M)``.
    """
    out: list[int] = []
    carry = 0
    while a or b or carry:
        s = a % p + b % p + carry
        carry = 1 if s >= p else 0
        out.append(carry)
        a //= p
        b //= p
    return out


def factorize_small(m: int) -> dict[int, int]:
    """Prime factorization of a small integer by trial division."""
    if m <= 0:
        raise ValueError("factorize_small expects a positive integer")
    out: dict[int, int] = {}
    d = 2
    while d * d <= m:
        while m % d == 0:
            out[d] = out.get(d, 0) + 1
            m //= d
        d += 1 if d == 2 else 2
    if m > 1:
        out[m] = out.get(m, 0) + 1
    return out


_WHEEL_GAPS = (4, 2, 4, 2, 4, 6, 2, 6)


def spf_trial(n: int, bound: int | None = None) -> int | None:
    """Smallest prime factor of ``n`` by 2-3-5 wheel trial division.

    Returns ``None`` if no factor <= ``bound`` exists (when ``bound`` is given
    and smaller than ``isqrt(n)``), otherwise returns ``n`` itself for primes.
    """
    if n < 2:
        raise ValueError("spf_trial expects n >= 2")
    for p in (2, 3, 5):
        if n % p == 0:
            return p
    limit = isqrt(n)
    if bound is not None:
        limit = min(limit, bound)
    d, i = 7, 0
    while d <= limit:
        if n % d == 0:
            return d
        d += _WHEEL_GAPS[i]
        i = (i + 1) % 8
    if bound is not None and bound < isqrt(n):
        return None
    return n


def pollard_rho(n: int, rng: random.Random | None = None) -> int:
    """Return a non-trivial divisor of composite ``n`` (Brent's variant)."""
    if n % 2 == 0:
        return 2
    rng = rng or random.Random(0xC0FFEE)
    while True:
        c = rng.randrange(1, n)
        f = lambda x: (x * x + c) % n  # noqa: E731
        x = y = rng.randrange(0, n)
        d = 1
        while d == 1:
            x = f(x)
            y = f(f(y))
            d = gcd(abs(x - y), n)
        if d != n:
            return d


def factorize(n: int) -> dict[int, int]:
    """Full factorization: trial division for small factors, then Pollard rho.

    Used as the *reference* factorizer for tests and benchmarks; the Pascal-row
    factorizer lives in :mod:`aksfactor.factor`.
    """
    if n < 1:
        raise ValueError("factorize expects a positive integer")
    out: dict[int, int] = {}
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
    if n == 1:
        return out
    stack = [n]
    while stack:
        m = stack.pop()
        if m == 1:
            continue
        if is_prime(m):
            out[m] = out.get(m, 0) + 1
            continue
        d = pollard_rho(m)
        stack.append(d)
        stack.append(m // d)
    return out
