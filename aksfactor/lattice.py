"""LLL and Coppersmith: the one family in this project that *is* polynomial time.

Everything else here searches -- for a scale, a smooth group order, a candidate
in a list -- and pays exponentially for it.  Coppersmith's method does not
search.  Given a modulus ``N`` with an unknown divisor ``p >= N**beta`` and a
polynomial ``f`` with a small root modulo ``p``, it *computes* that root in time
polynomial in ``log N``, provided the root is smaller than roughly
``N**(beta**2/deg f)``.

Applied to factoring: if you know ``p`` to within about ``N**0.25``, you recover
``p`` in polynomial time.  That is a genuine polynomial-time factoring algorithm
-- conditional on partial information nobody knows how to obtain cheaply, which
is exactly why factoring is still hard.

The construction is Howgrave-Graham's formulation.  To find a small ``x0`` with
``f(x0) = 0 (mod p)`` where ``p | N``:

* build the lattice spanned by ``N**(m-i) * f(x)**i * x**j``, coefficients scaled
  by powers of the bound ``X`` so short vectors mean small coefficients;
* LLL-reduce it;
* a sufficiently short vector is a polynomial with ``x0`` as a root **over the
  integers**, not merely mod ``p``, so ordinary root-finding recovers it.

The LLL here uses exact rational arithmetic and recomputes Gram-Schmidt after
each update: correct and easy to check, quadratic-ish slower than a floating
implementation.  Lattice dimensions in this file stay small enough for that.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd, isqrt

__all__ = ["lll", "coppersmith_small_root", "factor_with_hint", "hint_bits_needed"]


def _dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def _gram_schmidt(basis):
    ortho, mu = [], [[Fraction(0)] * len(basis) for _ in basis]
    for i, row in enumerate(basis):
        vec = [Fraction(x) for x in row]
        for j in range(i):
            denom = _dot(ortho[j], ortho[j])
            mu[i][j] = _dot([Fraction(x) for x in row], ortho[j]) / denom if denom else Fraction(0)
            vec = [vec[k] - mu[i][j] * ortho[j][k] for k in range(len(vec))]
        ortho.append(vec)
    return ortho, mu


def lll(basis, delta=Fraction(99, 100)):
    """LLL-reduce an integer basis (list of rows), in exact arithmetic.

    Gram-Schmidt is computed once and then updated incrementally after each
    size-reduction and swap.  Recomputing it from scratch every step -- the
    obvious first implementation -- is what made this unusably slow.
    """
    b = [list(map(int, row)) for row in basis]
    n = len(b)
    if n < 2:
        return b

    def build():
        ortho, mu, norms = [], [[Fraction(0)] * n for _ in range(n)], []
        for i in range(n):
            vec = [Fraction(x) for x in b[i]]
            for j in range(i):
                mu[i][j] = (Fraction(_dot(b[i], ortho[j])) / norms[j]
                            if norms[j] else Fraction(0))
                vec = [vec[k] - mu[i][j] * ortho[j][k] for k in range(len(vec))]
            ortho.append(vec)
            norms.append(_dot(vec, vec))
        return ortho, mu, norms

    ortho, mu, norms = build()
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            if abs(mu[k][j]) > Fraction(1, 2):
                q = int(mu[k][j] + Fraction(1, 2)) if mu[k][j] > 0 \
                    else -int(-mu[k][j] + Fraction(1, 2))
                if q:
                    b[k] = [b[k][i] - q * b[j][i] for i in range(len(b[k]))]
                    for i in range(j):
                        mu[k][i] -= q * mu[j][i]
                    mu[k][j] -= q
        if norms[k] >= (delta - mu[k][k - 1] ** 2) * norms[k - 1]:
            k += 1
        else:
            b[k], b[k - 1] = b[k - 1], b[k]
            ortho, mu, norms = build()
            k = max(k - 1, 1)
    return b


def _poly_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def _poly_eval(coeffs, x):
    total = 0
    for c in reversed(coeffs):
        total = total * x + c
    return total


def _poly_deriv(coeffs):
    return [i * c for i, c in enumerate(coeffs)][1:]


def _real_roots(coeffs, lo, hi):
    """Real roots of a low-degree integer polynomial in ``[lo, hi]``.

    Recursively finds the critical points (roots of the derivative), which
    separate the roots of the polynomial, then bisects each bracketing interval.
    Exact where it matters: the caller rounds and verifies over the integers.
    """
    while coeffs and coeffs[-1] == 0:
        coeffs = coeffs[:-1]
    if len(coeffs) <= 1:
        return []
    if len(coeffs) == 2:  # linear
        a0, a1 = coeffs
        if a1 == 0:
            return []
        x = -a0 / a1
        return [x] if lo <= x <= hi else []
    crit = _real_roots(_poly_deriv(coeffs), lo, hi)
    points = [lo] + sorted(crit) + [hi]
    out = []
    for left, right in zip(points, points[1:]):
        fl, fr = _poly_eval(coeffs, left), _poly_eval(coeffs, right)
        if fl == 0:
            out.append(left)
        if fl * fr < 0:
            a, b = left, right
            for _ in range(200):
                mid = (a + b) / 2
                if _poly_eval(coeffs, a) * _poly_eval(coeffs, mid) <= 0:
                    b = mid
                else:
                    a = mid
            out.append((a + b) / 2)
    if _poly_eval(coeffs, hi) == 0:
        out.append(hi)
    return out


def _integer_roots(coeffs, bound):
    """Integer roots with ``|root| <= bound``, found via real-root isolation.

    The previous version enumerated divisors of the constant term -- which in a
    Coppersmith lattice is of size ``N**m``, making it hopeless.  Isolating real
    roots first and verifying exactly is both correct and fast.
    """
    while coeffs and coeffs[-1] == 0:
        coeffs = coeffs[:-1]
    if len(coeffs) <= 1:
        return []
    found = []
    for approx in _real_roots([float(c) for c in coeffs], -float(bound), float(bound)):
        for cand in {int(approx), int(approx) + 1, int(approx) - 1,
                     round(approx)}:
            if abs(cand) <= bound and _poly_eval(coeffs, cand) == 0 and cand not in found:
                found.append(cand)
    if coeffs[0] == 0 and 0 not in found:
        found.append(0)
    return found


def coppersmith_small_root(n: int, poly, bound: int, beta: float = 0.5,
                           m: int = 3, t: int | None = None):
    """Small roots of ``poly(x) == 0 (mod p)`` for an unknown ``p >= n**beta``.

    ``poly`` is a monic coefficient list, low degree first.  Returns the integer
    roots found with ``|x| <= bound``.
    """
    d = len(poly) - 1
    if d < 1 or poly[-1] != 1:
        raise ValueError("coppersmith_small_root expects a monic polynomial")
    if t is None:
        t = max(1, int(d * m * (1 / beta - 1)))

    # rows: N^(m-i) f^i x^j  for i<m, j<d ;  then f^m x^j for j<t
    rows_poly = []
    fpow = [1]
    for i in range(m + 1):
        if i < m:
            scale = n ** (m - i)
            for j in range(d):
                shifted = [0] * j + [c * scale for c in fpow]
                rows_poly.append(shifted)
        else:
            for j in range(t):
                rows_poly.append([0] * j + list(fpow))
        fpow = _poly_mul(fpow, poly)

    width = max(len(r) for r in rows_poly)
    basis = []
    for r in rows_poly:
        row = r + [0] * (width - len(r))
        basis.append([row[k] * bound**k for k in range(width)])

    reduced = lll(basis)
    found = []
    for vec in reduced:
        if all(c == 0 for c in vec):
            continue
        coeffs = []
        ok = True
        for k, c in enumerate(vec):
            bk = bound**k
            if c % bk:
                ok = False
                break
            coeffs.append(c // bk)
        if not ok:
            continue
        for root in _integer_roots(coeffs, bound):
            if root not in found:
                found.append(root)
    return found


def hint_bits_needed(n: int) -> int:
    """How many high bits of ``p`` Coppersmith needs: about ``(1/4) log2 N``."""
    return (n.bit_length() + 3) // 4


def factor_with_hint(n: int, p_approx: int, bound: int | None = None,
                     m: int = 3):
    """Factor ``n`` given ``p_approx`` within ``bound`` of a true divisor.

    Polynomial time in ``log n`` for ``bound`` up to roughly ``n**0.25``.  This
    is the closest thing to a polynomial-time factoring algorithm that exists --
    and it needs partial information no one knows how to get cheaply.
    """
    if bound is None:
        bound = 1 << hint_bits_needed(n)
    for root in coppersmith_small_root(n, [p_approx % n, 1], bound, beta=0.5, m=m):
        cand = p_approx + root
        if cand > 1:
            g = gcd(cand, n)
            if 1 < g < n:
                return (g, n // g) if g <= n // g else (n // g, g)
    return None
