"""Round 20: Schnorr's prime-number lattice, tested at small sizes.

Schnorr (2021, "Fast factoring integers by SVP algorithms") proposed finding
the smooth relations of a sieve with lattice reduction.  For primes
``p_1 .. p_n`` take the lattice spanned by

    b_i = (0, .., f_i, .., 0, round(C ln p_i)),      target t = (0, .., 0, round(C ln N)),

with ``f`` a random permutation of ``ceil(i/2)``.  A lattice vector
``sum e_i b_i`` near ``t`` has ``sum e_i ln p_i ~ ln N``: with
``u = prod_{e_i > 0} p_i^e_i`` and ``v = prod_{e_i < 0} p_i^-e_i``, ``u/v ~ N``
and ``r = u - vN`` is smallish.  If ``r`` is smooth too, ``u = r (mod N)`` is a
relation, and enough of them give ``X^2 = 1 (mod N)`` by linear algebra mod 2.

The question the round measures: are lattice-found ``r`` smooth more often
than integers of the same size?  One structural fact decides a lot of it:
``u`` and ``v`` use disjoint primes, and ``r`` is divisible by none of them
(``p | u`` forces ``r = -vN != 0 mod p``), so ``r`` can only be smooth over
the primes the lattice vector did *not* use (``support_coprime``).
"""

from __future__ import annotations

import math
from math import gcd


def lll_integral(basis, a: int = 99, b: int = 100):
    """Exact LLL with ``delta = a/b`` in integer arithmetic (Cohen, Algorithm
    2.6.7): ``d_i`` are Gram determinants and ``lam[i][j] = d_j mu_ij``."""
    B = [list(r) for r in basis]
    n = len(B)

    def dot(u, v):
        return sum(x * y for x, y in zip(u, v))

    d = [0] * (n + 1)
    d[0] = 1
    lam = [[0] * n for _ in range(n)]

    def red(k, l):
        if 2 * abs(lam[k][l]) > d[l + 1]:
            q = (2 * lam[k][l] + d[l + 1]) // (2 * d[l + 1])
            B[k] = [x - q * y for x, y in zip(B[k], B[l])]
            lam[k][l] -= q * d[l + 1]
            for i in range(l):
                lam[k][i] -= q * lam[l][i]

    def swap(k, kmax):
        B[k], B[k - 1] = B[k - 1], B[k]
        for j in range(k - 1):
            lam[k][j], lam[k - 1][j] = lam[k - 1][j], lam[k][j]
        lm = lam[k][k - 1]
        bn = (d[k - 1] * d[k + 1] + lm * lm) // d[k]
        for i in range(k + 1, kmax + 1):
            t = lam[i][k]
            lam[i][k] = (d[k + 1] * lam[i][k - 1] - lm * t) // d[k]
            lam[i][k - 1] = (bn * t + lm * lam[i][k]) // d[k + 1]
        d[k] = bn

    d[1] = dot(B[0], B[0])
    k, kmax = 1, 0
    while k < n:
        if k > kmax:
            kmax = k
            for j in range(k + 1):
                u = dot(B[k], B[j])
                for i in range(j):
                    u = (d[i + 1] * u - lam[k][i] * lam[j][i]) // d[i]
                if j < k:
                    lam[k][j] = u
                else:
                    d[k + 1] = u
                    if u == 0:
                        raise ValueError("dependent vectors")
        red(k, k - 1)
        if b * d[k + 1] * d[k - 1] < a * d[k] * d[k] - b * lam[k][k - 1] ** 2:
            swap(k, kmax)
            k = max(1, k - 1)
        else:
            for l in range(k - 2, -1, -1):
                red(k, l)
            k += 1
    return B


def prime_lattice_relations(n: int, primes, c: float, rng, embed: int = 1):
    """Candidate ``(u, v, r)`` from one randomised prime-number lattice.

    Kannan's embedding: append ``(t, embed)`` to the basis, LLL-reduce, and read
    close vectors off the reduced rows whose last coordinate is ``+-embed``
    (and their sums and differences with rows ending in ``0``).
    """
    k = len(primes)
    scale = int(round(n ** c))
    f = [(i + 2) // 2 for i in range(k)]
    rng.shuffle(f)
    rows = []
    for i in range(k):
        row = [0] * (k + 2)
        row[i] = f[i]
        row[k] = round(scale * math.log(primes[i]))
        rows.append(row)
    target = [0] * k + [round(scale * math.log(n)), embed]
    rows.append(target)
    reduced = lll_integral(rows)
    pool = [w for w in reduced if abs(w[-1]) == embed]
    vectors = list(pool)
    for w in reduced:
        if w[-1] == 0:
            for z in pool[:3]:
                vectors.append([x + y for x, y in zip(z, w)])
                vectors.append([x - y for x, y in zip(z, w)])
    out = []
    for w in vectors:
        s = 1 if w[-1] == embed else -1
        v = [t - s * x for t, x in zip(target[:k + 1], w[:k + 1])]
        if any(v[i] % f[i] for i in range(k)):
            continue
        e = [v[i] // f[i] for i in range(k)]
        u = math.prod(primes[i] ** e[i] for i in range(k) if e[i] > 0)
        w_ = math.prod(primes[i] ** -e[i] for i in range(k) if e[i] < 0)
        r = u - w_ * n
        if r:
            out.append((u, w_, r))
    return out


def support_coprime(u: int, v: int, r: int, primes) -> bool:
    """``r = u - vN`` shares no prime with ``u v`` (for ``gcd(uv, N) = 1``)."""
    return all(r % p for p in primes if u % p == 0 or v % p == 0)


def smooth_exponents(m: int, primes):
    """``[sign, e_1, .., e_k]`` of ``m`` over ``primes``, or ``None``."""
    e = [1 if m < 0 else 0]
    m = abs(m)
    for p in primes:
        k = 0
        while m % p == 0:
            m //= p
            k += 1
        e.append(k)
    return e if m == 1 else None


def split_from_relations(n: int, primes, relations):
    """From relations ``u = r (mod n)`` with ``u, r`` smooth, find a factor.

    A GF(2) dependency among ``exps(u) - exps(r)`` makes ``prod u / prod r`` a
    rational square ``Z^2``, and ``prod u = prod r (mod n)`` gives
    ``Z^2 = 1 (mod n)``; ``gcd(Z - 1, n)`` splits ``n`` about half the time.
    """
    vecs = []
    for u, r in relations:
        a, b = smooth_exponents(u, primes), smooth_exponents(r, primes)
        vecs.append([x - y for x, y in zip(a, b)])
    width = len(primes) + 1
    pivots = {}
    for i, v in enumerate(vecs):
        mask = sum(1 << j for j, x in enumerate(v) if x & 1)
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
        total = [0] * width
        for j in range(len(vecs)):
            if (comb >> j) & 1:
                total = [x + y for x, y in zip(total, vecs[j])]
        z = 1
        for j, p in enumerate(primes):
            e = total[j + 1] // 2
            z = z * pow(p if e >= 0 else pow(p, -1, n), abs(e), n) % n
        g = gcd(z - 1, n)
        if 1 < g < n:
            return g
    return None


def _local_weight(x: int, used: set, small) -> float:
    """Density of the residues' local law relative to uniform integers.

    At a prime ``p`` the vector used, ``r`` is a unit (weight ``p/(p-1)`` off
    zero, none at zero).  At an unused odd ``p``, ``u`` and ``v`` are units, so
    ``p^k | r`` with probability ``1/phi(p^k)`` instead of ``p^-k``: weight
    ``p/(p-1)`` on multiples of ``p`` and ``p(p-2)/(p-1)^2`` elsewhere.  At an
    unused ``2``, ``r`` is even and its 2-adic valuation is exactly that of a
    random even number: weight ``2`` on evens.
    """
    w = 1.0
    for p in small:
        div = x % p == 0
        if p in used:
            if div:
                return 0.0
            w *= p / (p - 1)
        elif p == 2:
            if not div:
                return 0.0
            w *= 2.0
        else:
            w *= p / (p - 1) if div else p * (p - 2) / (p - 1) ** 2
    return w


def residue_baselines(r: int, u: int, v: int, primes, rng, k: int = 60,
                      small_bound: int = 60) -> tuple[float, float]:
    """Smoothness rates of integers within 10% of ``|r|``: plain, and with the
    local law of lattice residues at every prime below ``small_bound`` (and
    coprimality to every prime the vector used)."""
    small = [p for p in primes if p < small_bound]
    used = {p for p in primes if u % p == 0 or v % p == 0}
    a = abs(r)
    lo, hi = max(2, int(a * 0.9)), int(a * 1.1) + 2
    plain = sum(smooth_exponents(rng.randrange(lo, hi), primes) is not None
                for _ in range(k)) / k
    top = 1.0
    for p in small:
        top *= p / (p - 1) if (p in used or p != 2) else 2.0
    hits = got = 0
    while got < k:
        x = rng.randrange(lo, hi)
        if any(x % p == 0 for p in used if p >= small_bound):
            continue
        if rng.random() * top > _local_weight(x, used, small):
            continue
        got += 1
        hits += smooth_exponents(x, primes) is not None
    return plain, hits / k
