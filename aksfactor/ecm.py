"""Round 22: Lenstra's elliptic curve method -- the order mechanism, redrawn.

The Pascal row, the AKS fold, Pollard's ``p - 1`` and Williams' ``p + 1`` all
work in a group whose order mod ``p`` is fixed by ``p`` (``p - 1``, ``p + 1``,
``p^d - 1``): one ticket per prime, and it wins only if that one number is
smooth (Proposition 14).  Lenstra's idea is to draw the group at random: the
curve ``E(F_p)`` has order ``p + 1 - t`` with ``|t| <= 2 sqrt p`` spread across
the Hasse interval as the curve varies, so every new curve is a new ticket.

Montgomery curves ``B y^2 = x^3 + A x^2 + x`` in ``(X : Z)`` coordinates with
Suyama's parametrisation (group order divisible by 12), a Montgomery-ladder
stage 1 to ``B1`` and a baby-step giant-step stage 2 to ``B2``.
"""

from __future__ import annotations

from math import gcd, log

from .arith import sieve


class FoundFactor(Exception):
    def __init__(self, g):
        super().__init__(g)
        self.g = g


def _inv(a: int, n: int) -> int:
    g = gcd(a, n)
    if g != 1:
        raise FoundFactor(g)
    return pow(a, -1, n)


def suyama_curve(sigma: int, n: int):
    """``(a24, X, Z)`` for Suyama's curve with parameter ``sigma``."""
    u = (sigma * sigma - 5) % n
    v = 4 * sigma % n
    x = pow(u, 3, n)
    z = pow(v, 3, n)
    num = pow(v - u, 3, n) * (3 * u + v) % n
    den = 16 * x * v % n
    a24 = num * _inv(den, n) % n          # (A + 2) / 4
    return a24, x, z


def _double(x, z, a24, n):
    s = (x + z) % n
    d = (x - z) % n
    s2, d2 = s * s % n, d * d % n
    t = (s2 - d2) % n
    return s2 * d2 % n, t * (d2 + a24 * t) % n


def _add(x1, z1, x2, z2, xd, zd, n):
    """``P1 + P2`` given ``P1 - P2 = (xd : zd)``."""
    a = (x1 - z1) * (x2 + z2) % n
    b = (x1 + z1) * (x2 - z2) % n
    s, d = (a + b) % n, (a - b) % n
    return zd * s * s % n, xd * d * d % n


def ladder(k: int, x, z, a24, n):
    """``k (x : z)`` by the Montgomery ladder."""
    if k == 0:
        return 0, 0
    x0, z0 = x, z
    x1, z1 = _double(x, z, a24, n)
    for bit in bin(k)[3:]:
        if bit == "1":
            x0, z0 = _add(x1, z1, x0, z0, x, z, n)
            x1, z1 = _double(x1, z1, a24, n)
        else:
            x1, z1 = _add(x0, z0, x1, z1, x, z, n)
            x0, z0 = _double(x0, z0, a24, n)
    return x0, z0


def _stage1(n, x, z, a24, b1):
    for p in sieve(b1):
        e = int(log(b1) / log(p))
        x, z = ladder(p ** e, x, z, a24, n)
    return x, z


def _stage2(n, x, z, a24, b1, b2, w=210):
    """Baby steps ``j Q`` (odd ``j < w/2``, ``gcd(j, w) = 1``), giant steps
    ``G_k = k w Q``; a prime ``l = k w +- j`` kills ``Q`` mod ``p`` exactly when
    ``x(G_k) = x(j Q)``, so accumulate ``X_G Z_j - X_j Z_G`` over those pairs."""
    if b1 < 2 * w:
        return 1
    x2, z2 = _double(x, z, a24, n)
    baby = {1: (x, z)}
    prev, cur, j = (x, z), _add(x2, z2, x, z, x, z, n), 3      # Q, 3Q
    while j < w // 2:
        if gcd(j, w) == 1:
            baby[j] = cur
        prev, cur = cur, _add(cur[0], cur[1], x2, z2, prev[0], prev[1], n)
        j += 2
    is_prime = bytearray(b2 + w + 1)
    for q in sieve(b2 + w):
        if q > b1:
            is_prime[q] = 1
    wx, wz = ladder(w, x, z, a24, n)
    k = b1 // w
    gpx, gpz = ladder((k - 1) * w, x, z, a24, n)
    gx, gz = ladder(k * w, x, z, a24, n)
    acc = 1
    while k * w - w // 2 <= b2:
        for j, (bx, bz) in baby.items():
            if is_prime[k * w + j] or is_prime[k * w - j]:
                acc = acc * (gx * bz - bx * gz) % n
        gpx, gpz, (gx, gz) = gx, gz, _add(gx, gz, wx, wz, gpx, gpz, n)
        k += 1
    return gcd(acc, n)


def ecm_curve(n: int, sigma: int, b1: int, b2: int | None = None):
    """Run one curve; return a factor of ``n`` or ``None``."""
    try:
        a24, x, z = suyama_curve(sigma, n)
        x, z = _stage1(n, x, z, a24, b1)
    except FoundFactor as f:
        return f.g if f.g < n else None
    g = gcd(z, n)
    if 1 < g < n:
        return g
    if g == n or b2 is None:
        return None
    g = _stage2(n, x, z, a24, b1, b2)
    return g if 1 < g < n else None


def ecm(n: int, b1: int = 2000, b2: int | None = None, curves: int = 500,
        rng=None, stats: dict | None = None):
    """Lenstra's ECM: try ``curves`` random curves; a factor of ``n`` or ``None``."""
    import random

    rng = rng or random.Random(1)
    b2 = b2 if b2 is not None else 50 * b1
    for c in range(1, curves + 1):
        g = ecm_curve(n, rng.randrange(6, n - 1), b1, b2)
        if g:
            if stats is not None:
                stats["curves"] = c
            return g
    if stats is not None:
        stats["curves"] = curves
    return None


def pminus1(n: int, b1: int, b2: int | None = None, base: int = 2):
    """Pollard's ``p - 1`` with the same ``B1``/``B2`` budget: one fixed group."""
    a = base
    for p in sieve(b1):
        a = pow(a, p ** int(log(b1) / log(p)), n)
    g = gcd(a - 1, n)
    if 1 < g < n:
        return g
    if b2 is None or g == n:
        return None
    acc, ap, last = 1, 1, 0
    for q in sieve(b2):
        if q <= b1:
            continue
        ap = ap * pow(a, q - last, n) % n          # a^q, stepping between primes
        last = q
        acc = acc * (ap - 1) % n
    g = gcd(acc, n)
    return g if 1 < g < n else None


def is_b1b2_smooth(m: int, b1: int, b2: int) -> bool:
    """Whether ``m`` has every prime power ``<= b1`` except possibly one prime
    in ``(b1, b2]`` to the first power -- exactly what stage 1 + stage 2 kill."""
    from .arith import factorize

    big = 0
    for ell, k in factorize(m).items():
        if ell ** k <= b1:
            continue
        if k == 1 and ell <= b2 and not big:
            big = ell
            continue
        return False
    return True
