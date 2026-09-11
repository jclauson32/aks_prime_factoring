"""A small quadratic sieve: the sign mechanism with smoothness behind it.

Every Pascal-native method in this repository lands at ``N^(1/4)``.  The
quadratic sieve is here as the reference point for what the *sign* mechanism
buys once smooth numbers do the work: relations

    (x + m)^2  =  Q(x)  (mod N),        Q(x) = (x + m)^2 - N,   m = ceil(sqrt N),

with ``Q(x)`` smooth over the primes for which ``N`` is a square, are combined
by linear algebra mod 2 into ``X^2 = Y^2 (mod N)``.  ``|Q(x)| ~ 2 x sqrt N`` is
small, which is the whole point, and sieving finds the smooth ones in bulk.

Single polynomial, no large primes, pure Python: a reference, not a record.
"""

from __future__ import annotations

from math import exp, gcd, isqrt, log, sqrt

from .arith import sieve


def sqrt_mod_prime(a: int, p: int) -> int:
    """A square root of ``a`` modulo the odd prime ``p`` (Tonelli-Shanks)."""
    a %= p
    if a == 0:
        return 0
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r


def default_parameters(n: int) -> tuple[int, int]:
    """Factor-base bound ``B ~ L^(1/2)`` and block length, tuned for Python."""
    ln = log(n)
    b = int(exp(0.55 * sqrt(ln * log(ln)))) + 60
    return b, max(20000, 40 * b)


def quadratic_sieve(n: int, bound: int | None = None, half_width: int | None = None,
                    stats: dict | None = None):
    """A non-trivial factor of the odd composite ``n`` (not a prime power), or None."""
    if n % 2 == 0:
        return 2
    r = isqrt(n)
    if r * r == n:
        return r
    b, width = default_parameters(n)
    bound = bound or b
    half_width = half_width or width
    base = [2] + [p for p in sieve(bound)[1:] if pow(n, (p - 1) // 2, p) == 1]
    for p in sieve(bound):
        if n % p == 0 and p < n:
            return p
    m = r + 1
    roots = {}
    for p in base[1:]:
        t = sqrt_mod_prime(n, p)
        roots[p] = {(t - m) % p, (-t - m) % p}
    scale = 64
    logs = {p: int(round(scale * log(p))) for p in base}
    slack = scale * 2.2 * log(bound)
    need = len(base) + 10
    relations = []
    pos, neg, rounds, candidates = 0, 0, 0, 0
    while len(relations) < need and rounds < 80:
        if rounds % 2 == 0:
            lo, pos = pos, pos + half_width
        else:
            lo, neg = neg - half_width, neg - half_width
        rounds += 1
        acc = [0] * half_width
        for p in base[1:]:
            lg = logs[p]
            for root in roots[p]:
                for i in range((root - lo) % p, half_width, p):
                    acc[i] += lg
        for i in range(half_width):
            x = lo + i
            q = (x + m) * (x + m) - n
            if q == 0 or acc[i] < scale * log(abs(q)) - slack:
                continue
            candidates += 1
            e = [1 if q < 0 else 0]
            v = abs(q)
            for p in base:
                k = 0
                while v % p == 0:
                    v //= p
                    k += 1
                e.append(k)
            if v == 1:
                relations.append((x + m, e))
    if stats is not None:
        stats.update(base=len(base), relations=len(relations), blocks=rounds,
                     sieved=rounds * half_width, candidates=candidates)
    if len(relations) <= len(base):
        return None
    return _combine(n, base, relations)


def _combine(n, base, relations):
    width = len(base) + 1
    pivots = {}
    for i, (_, e) in enumerate(relations):
        mask = sum(1 << j for j, x in enumerate(e) if x & 1)
        comb = 1 << i
        for bit in range(width):
            if not (mask >> bit) & 1:
                continue
            if bit in pivots:
                pm, pc = pivots[bit]
                mask ^= pm
                comb ^= pc
            else:
                pivots[bit] = (mask, comb)
                break
        if mask:
            continue
        x_prod, total = 1, [0] * width
        for j, (a, e) in enumerate(relations):
            if (comb >> j) & 1:
                x_prod = x_prod * a % n
                total = [s + t for s, t in zip(total, e)]
        y = 1
        for j, p in enumerate(base):
            y = y * pow(p, total[j + 1] // 2, n) % n
        g = gcd(x_prod - y, n)
        if 1 < g < n:
            return g
    return None
