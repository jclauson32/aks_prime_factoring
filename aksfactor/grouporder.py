"""Counting the group orders each construction can reach at a *fixed* prime.

Proposition 14 claims that every method reachable from the AKS ring has a
**rigid** group order -- fixed once ``p`` is -- while elliptic curves supply a
*family* of orders at the same ``p``, and that this is the structural difference
that separates the two.

This module makes that countable.  For a fixed prime ``p``:

``ring_unit_order``
    the order of ``(F_p[x]/f)*`` for a squarefree ``f``, which is
    ``prod_i (p**deg(f_i) - 1)`` over the irreducible factors.  Varying ``f``
    over degree ``d`` can only permute the degree multiset, so the number of
    distinct orders is bounded by the number of partitions of ``d``.
``elliptic_order``
    ``#E(F_p)`` for ``y**2 = x**3 + a x + b``, which ranges over an interval of
    width ``4*sqrt(p)`` as the curve varies.
"""

from __future__ import annotations

from math import gcd as _igcd

__all__ = [
    "poly_gcd",
    "degree_pattern",
    "ring_unit_order",
    "elliptic_order",
    "partitions",
]


# --------------------------------------------------------------------------
# polynomials over F_p  (p prime, so coefficients are invertible)
# --------------------------------------------------------------------------

def _trim(a: list[int]) -> list[int]:
    while a and a[-1] == 0:
        a.pop()
    return a


def _mul(a: list[int], b: list[int], p: int) -> list[int]:
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    out[i + j] = (out[i + j] + ai * bj) % p
    return _trim(out)


def _divmod(a: list[int], b: list[int], p: int):
    a = a[:]
    db = len(b) - 1
    if db < 0:
        raise ZeroDivisionError
    inv = pow(b[db], -1, p)
    q = [0] * max(0, len(a) - db)
    for i in range(len(a) - 1, db - 1, -1):
        c = a[i] * inv % p
        q[i - db] = c
        if c:
            for j in range(db + 1):
                a[i - db + j] = (a[i - db + j] - c * b[j]) % p
    return _trim(q), _trim(a[:db])


def poly_gcd(a: list[int], b: list[int], p: int) -> list[int]:
    a, b = _trim(a[:]), _trim(b[:])
    while b:
        a, b = b, _divmod(a, b, p)[1]
    if a:
        inv = pow(a[-1], -1, p)
        a = [c * inv % p for c in a]
    return a


def _powmod(base: list[int], e: int, mod: list[int], p: int) -> list[int]:
    result = [1 % p]
    base = _divmod(base, mod, p)[1]
    while e:
        if e & 1:
            result = _divmod(_mul(result, base, p), mod, p)[1]
        e >>= 1
        if e:
            base = _divmod(_mul(base, base, p), mod, p)[1]
    return result


def degree_pattern(f: list[int], p: int) -> list[int] | None:
    """Degrees of the irreducible factors of a **squarefree** monic ``f``.

    Returns ``None`` when ``f`` is not squarefree (distinct-degree factorization
    assumes it, and the unit-group formula below changes when it fails).
    """
    f = _trim(f[:])
    d = len(f) - 1
    if d < 1:
        return []
    deriv = _trim([(i * f[i]) % p for i in range(1, len(f))])
    if not deriv or len(poly_gcd(f, deriv, p)) - 1 > 0:
        return None
    degrees: list[int] = []
    star = f[:]
    h = [0, 1 % p]  # x
    i = 0
    while len(star) - 1 > 0:
        i += 1
        if i > (len(star) - 1):
            degrees.append(len(star) - 1)  # remaining factor is irreducible
            break
        h = _powmod(h, p, star, p)
        hx = _trim([(h[k] if k < len(h) else 0) - (1 if k == 1 else 0)
                    for k in range(max(len(h), 2))])
        hx = [c % p for c in hx]
        g = poly_gcd(hx, star, p)
        dg = len(g) - 1
        if dg > 0:
            degrees.extend([i] * (dg // i))
            star = _divmod(star, g, p)[0]
    return sorted(degrees)


def ring_unit_order(f: list[int], p: int) -> int | None:
    """``|(F_p[x]/f)*|`` for squarefree monic ``f``; ``None`` if not squarefree."""
    pattern = degree_pattern(f, p)
    if pattern is None:
        return None
    order = 1
    for d in pattern:
        order *= p**d - 1
    return order


def partitions(d: int) -> int:
    """Number of integer partitions of ``d`` -- the ceiling on distinct ring orders."""
    table = [1] + [0] * d
    for k in range(1, d + 1):
        for v in range(k, d + 1):
            table[v] += table[v - k]
    return table[d]


# --------------------------------------------------------------------------
# elliptic curves
# --------------------------------------------------------------------------

def elliptic_order(a: int, b: int, p: int) -> int | None:
    """``#E(F_p)`` for ``y^2 = x^3 + a x + b``, by direct counting.

    ``None`` when the curve is singular.  Cost is ``O(p)``, which is fine for the
    small primes this is used on.
    """
    if (4 * a**3 + 27 * b**2) % p == 0:
        return None
    chi = [0] * p
    for y in range(1, (p + 1) // 2 + 1):
        chi[y * y % p] = 1
    total = 1  # point at infinity
    for x in range(p):
        v = (x * x % p * x + a * x + b) % p
        if v == 0:
            total += 1
        elif chi[v]:
            total += 2
    return total
