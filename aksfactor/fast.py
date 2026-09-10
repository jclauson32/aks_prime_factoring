"""Finding the first non-zero Pascal position in ``O~(n**0.25)`` instead of ``O~(n**0.5)``.

Theorem 2 says the first non-zero interior entry of row ``n`` mod ``n`` is the
first ``k`` with ``gcd(k,n) > 1``.  Read naively that is trial division, one
position at a time, costing ``Theta(spf(n))``.

But a *block* of positions can be tested with a single gcd: if

    f(X) = (X+1)(X+2)...(X+c)

then ``gcd(f(i*c) mod n, n) > 1`` exactly when some position in
``(i*c, i*c + c]`` shares a factor with ``n``.  Evaluating one degree-``c``
polynomial at ``c`` points is a *fast multipoint evaluation* -- ``O~(c)`` ring
operations via a product tree and a remainder tree -- so ``c**2`` positions cost
``O~(c)``.  Choosing ``c = sqrt(B)`` searches all positions up to ``B`` in
``O~(sqrt(B))``, i.e. ``O~(n**0.25)`` for ``B = sqrt(n)``.

This is Strassen's deterministic factoring bound, reached here through the
Pascal-row characterisation, and it is the best deterministic bound known.  It
does not make the method polynomial; it halves the exponent.
"""

from __future__ import annotations

from math import gcd, isqrt

from .arith import spf_trial
from .pascal import _pack, _unpack

__all__ = ["fast_spf", "fast_split", "poly_mul", "multipoint_eval", "product_of",
           "COUNTERS", "reset_counters"]

# Below this search bound the wheel beats the polynomial machinery outright.
_CROSSOVER = 1 << 14


# --------------------------------------------------------------------------
# Polynomial arithmetic over Z/n  (coefficient i = coefficient of x**i)
# --------------------------------------------------------------------------

def _trim(a: list[int]) -> list[int]:
    while a and a[-1] == 0:
        a.pop()
    return a


COUNTERS = {"poly_mul": 0, "coeff_mul": 0}


def reset_counters() -> None:
    COUNTERS["poly_mul"] = 0
    COUNTERS["coeff_mul"] = 0


def poly_mul(a: list[int], b: list[int], n: int) -> list[int]:
    """Full product of two polynomials over ``Z/n``."""
    if not a or not b:
        return []
    COUNTERS["poly_mul"] += 1
    COUNTERS["coeff_mul"] += len(a) * len(b)
    if min(len(a), len(b)) <= 24:
        out = [0] * (len(a) + len(b) - 1)
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b):
                    if bj:
                        out[i + j] = (out[i + j] + ai * bj) % n
        return _trim(out)
    slack = max(len(a), len(b)).bit_length()
    wbytes = (2 * n.bit_length() + slack + 8) // 8 + 1
    prod = _pack(a, wbytes) * _pack(b, wbytes)
    return _trim(_unpack(prod, wbytes, len(a) + len(b) - 1, n))


def _pad(a: list[int], k: int) -> list[int]:
    """Exactly ``k`` coefficients.  Length carries meaning wherever a coefficient
    list gets reversed, so trimmed trailing zeros must be restored."""
    a = a[:k]
    return a + [0] * (k - len(a))


def _poly_mul_trunc(a: list[int], b: list[int], k: int, n: int) -> list[int]:
    return _pad(poly_mul(a, b, n), k)


def _series_inv(a: list[int], k: int, n: int) -> list[int]:
    """Inverse of ``a`` as a power series to precision ``k``.  Requires ``a[0] == 1``.

    Every call site here inverts the reversal of a *monic* polynomial, whose
    constant term is the leading coefficient ``1`` -- so no modular inverse mod
    the composite ``n`` is ever needed.
    """
    if a[0] != 1 % n:
        raise ValueError("_series_inv needs a[0] == 1")
    g = [1 % n]
    prec = 1
    while prec < k:
        prec = min(2 * prec, k)
        # g <- g * (2 - a*g)  mod x**prec
        ag = _poly_mul_trunc(_pad(a, prec), g, prec, n)
        two_minus = [(-v) % n for v in ag]
        two_minus[0] = (two_minus[0] + 2) % n
        g = _poly_mul_trunc(g, two_minus, prec, n)
    return _pad(g, k)


def _poly_rem_monic(a: list[int], b: list[int], n: int) -> list[int]:
    """Remainder of ``a`` divided by the **monic** polynomial ``b``."""
    m, t = len(a) - 1, len(b) - 1
    if t < 0:
        raise ZeroDivisionError("division by the zero polynomial")
    if m < t:
        return a[:]
    if t == 0:
        return []
    ra = _pad(a, m + 1)[::-1]
    rb = _pad(b, t + 1)[::-1]
    inv = _series_inv(rb, m - t + 1, n)
    rq = _poly_mul_trunc(ra, inv, m - t + 1, n)
    q = rq[::-1]
    bq = poly_mul(b, q, n)
    r = [(a[i] - (bq[i] if i < len(bq) else 0)) % n for i in range(t)]
    return _trim(r)


def product_of(polys: list[list[int]], n: int) -> list[int]:
    """Balanced product of a list of polynomials (a product tree)."""
    if not polys:
        return [1 % n]
    while len(polys) > 1:
        polys = [
            poly_mul(polys[i], polys[i + 1], n) if i + 1 < len(polys) else polys[i]
            for i in range(0, len(polys), 2)
        ]
    return polys[0]


def _subproduct_tree(points: list[int], n: int) -> list[list[list[int]]]:
    """Levels of the tree of ``prod (X - x_i)``; level 0 is the leaves."""
    nodes = [[(-x) % n, 1 % n] for x in points]
    tree = [nodes]
    while len(nodes) > 1:
        nodes = [poly_mul(nodes[i], nodes[i + 1], n) for i in range(0, len(nodes), 2)]
        tree.append(nodes)
    return tree


def multipoint_eval(f: list[int], points: list[int], n: int) -> list[int]:
    """Evaluate ``f`` at every point, in ``O~(deg f + #points)`` ring operations."""
    if not points:
        return []
    m = len(points)
    size = 1 << (m - 1).bit_length()
    padded = points + [points[-1]] * (size - m)
    tree = _subproduct_tree(padded, n)
    cur = [_poly_rem_monic(f, tree[-1][0], n)]
    for level in range(len(tree) - 2, -1, -1):
        nxt = []
        for i, rem in enumerate(cur):
            nxt.append(_poly_rem_monic(rem, tree[level][2 * i], n))
            nxt.append(_poly_rem_monic(rem, tree[level][2 * i + 1], n))
        cur = nxt
    return [(c[0] if c else 0) for c in cur[:m]]


# --------------------------------------------------------------------------
# The search
# --------------------------------------------------------------------------

def fast_spf(n: int, bound: int | None = None, stats: dict | None = None):
    """Smallest prime factor of ``n`` in ``O~(sqrt(bound))`` operations.

    Equivalent to `aksfactor.factor.pascal_spf` in what it returns -- the
    position of the first non-zero entry of row ``n`` mod ``n`` -- but it tests
    ``c`` positions per gcd instead of one.

    Returns ``n`` for primes, or ``None`` if ``bound`` is set and no factor
    ``<= bound`` exists.
    """
    if n < 2:
        raise ValueError("fast_spf expects n >= 2")
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return p
    limit = isqrt(n) if bound is None else min(bound, isqrt(n))
    c = isqrt(limit) + 1
    if stats is not None:
        stats.update({"c": c, "positions_covered": c * c, "gcds": 0,
                      "path": "multipoint" if limit >= _CROSSOVER else "wheel"})
    if limit < _CROSSOVER:
        # Below the crossover the wheel wins outright; the bound is asymptotic.
        return spf_trial(n, bound=bound)


    f = product_of([[j % n, 1 % n] for j in range(1, c + 1)], n)
    nblocks = limit // c + 1
    points = [i * c % n for i in range(nblocks)]
    values = multipoint_eval(f, points, n)

    for i, v in enumerate(values):
        if stats is not None:
            stats["gcds"] += 1
        if gcd(v, n) > 1:
            for j in range(1, c + 1):
                k = i * c + j
                if k > 1 and gcd(k, n) > 1:
                    return k  # smallest k > 1 sharing a factor with n, hence prime
    if bound is not None and bound < isqrt(n):
        return None
    return n


def fast_split(n: int, bound: int | None = None):
    """Split ``n`` as ``(p, n // p)`` using :func:`fast_spf`."""
    p = fast_spf(n, bound=bound)
    if p is None or p == n:
        return None
    return p, n // p
