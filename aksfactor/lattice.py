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

__all__ = ["lll", "coppersmith_small_root", "factor_with_hint", "hint_bits_needed",
           "cost_exponents", "factor_with_congruence", "congruence_modulus_needed",
           "residue_candidates", "crt_assembly_cost"]


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
                     m: int = 3, beta: float = 0.5):
    """Factor ``n`` given ``p_approx`` within ``bound`` of a true divisor.

    Polynomial time in ``log n``.  The reachable ``bound`` is ``n**(beta**2)``
    where ``beta`` is defined by ``p >= n**beta``, so a balanced semiprime
    (``beta = 1/2``) gives the familiar ``n**0.25`` and nothing gives more:
    ``beta <= 1/2`` always, since ``p`` is the smaller factor.

    This is the closest thing to a polynomial-time factoring algorithm that
    exists -- and it needs partial information no one knows how to get cheaply.
    """
    if bound is None:
        bound = 1 << hint_bits_needed(n)
    for root in coppersmith_small_root(n, [p_approx % n, 1], bound, beta=beta, m=m):
        cand = p_approx + root
        if cand > 1:
            g = gcd(cand, n)
            if 1 < g < n:
                return (g, n // g) if g <= n // g else (n // g, g)
    return None


def cost_exponents(beta: float) -> dict:
    """Exponents of the two ways to turn Coppersmith into a factoring algorithm.

    Write ``p ~ N**beta`` for the smaller factor, so ``beta <= 1/2`` always.

    * Coppersmith's window is ``N**(beta**2)`` (measured in
      ``experiments/exp20_equivalence.py``).  Covering the range of possible
      ``p``, of length ``~N**beta``, therefore needs ``N**(beta - beta**2)``
      guesses -- exponent ``beta(1 - beta)``.
    * Strassen finds ``p`` outright in ``O~(sqrt(p)) = O~(N**(beta/2))``.

    Comparing, ``beta(1-beta) <= beta/2`` iff ``beta >= 1/2``.  Since
    ``beta <= 1/2`` with equality only for a balanced semiprime, **guessing plus
    Coppersmith is never strictly better than Strassen**, and is strictly worse
    for every unbalanced semiprime.  The two coincide at exactly ``N**(1/4)``.
    """
    return {
        "beta": beta,
        "coppersmith_window": beta * beta,
        "guess_and_coppersmith": beta * (1 - beta),
        "strassen": beta / 2,
        "strassen_at_least_as_good": beta * (1 - beta) >= beta / 2 - 1e-12,
    }


def congruence_modulus_needed(n: int) -> int:
    """Smallest modulus ``M`` for which knowing ``p mod M`` suffices: ``~n**0.25``."""
    root = isqrt(isqrt(n))
    return root + 1


def factor_with_congruence(n: int, r: int, m_mod: int, m: int = 3,
                           beta: float = 0.5):
    """Factor ``n`` given ``p ≡ r (mod m_mod)`` for a known modulus.

    The *same* `N**(1/4)` budget in a different shape.  Instead of the high bits
    of ``p``, this consumes a congruence: writing ``p = r + m_mod * x`` the
    unknown ``x`` is bounded by ``sqrt(n)/m_mod``, so Coppersmith succeeds once
    ``m_mod >= n**0.25``.

    Making the polynomial monic is the only wrinkle: ``r + m_mod*x`` is not, so
    it is multiplied by ``m_mod**(-1) mod n`` first, which does not change the
    roots modulo ``p``.

    Returns ``(p, q)`` or ``None``.
    """
    if m_mod <= 1:
        return None
    g = gcd(m_mod, n)
    if 1 < g < n:
        return (g, n // g) if g <= n // g else (n // g, g)
    monic_const = (r % n) * pow(m_mod, -1, n) % n
    bound = isqrt(n) // m_mod + 1
    for root in coppersmith_small_root(n, [monic_const, 1], bound, beta=beta, m=m):
        cand = r + m_mod * root
        if cand > 1:
            g = gcd(cand, n)
            if 1 < g < n:
                return (g, n // g) if g <= n // g else (n // g, g)
    return None


def residue_candidates(n: int, ell: int) -> dict:
    """What ``N`` reveals about ``p mod ell`` for a small prime ``ell``.

    Over ``F_ell`` the residues of ``p`` and ``q`` are the roots of

        z**2 - s*z + N,       s = (p + q) mod ell,

    so ``p mod ell`` is determined exactly by ``s`` -- equivalently, since
    ``p + q = N + 1 - phi(N)``, by ``phi(N) mod ell``.  ``N mod ell`` is free;
    ``s`` is the entire unknown.

    Knowing only ``N`` leaves the possible sums ``{a + N/a : a in F_ell*}``.  The
    map ``a -> a + N/a`` is two-to-one (``a`` and ``N/a`` collide), so this set
    has about ``(ell-1)/2`` elements: **``N`` gives away exactly the ``p <-> q``
    symmetry and nothing else.**

    Returns the possible sums, the implied candidate count, and whether the
    residue happens to be free (a unique sum).
    """
    if ell < 2 or n % ell == 0:
        return {"ell": ell, "sums": [], "candidates": 0, "free": False}
    sums = sorted({(a + n * pow(a, -1, ell)) % ell for a in range(1, ell)})
    return {
        "ell": ell,
        "sums": sums,
        "candidates": len(sums),
        "symmetry_bound": (ell - 1) // 2,
        "free": len(sums) == 1,
    }


def crt_assembly_cost(n_bits: int, prime_bound: int = 4000) -> dict:
    """Cost of assembling ``p mod M`` with ``M >= N**0.25`` from small primes.

    Choosing ``p mod ell`` independently for each prime ``ell`` up to ``y``, with
    ``(ell-1)/2`` candidates apiece, the search is

        prod_{ell <= y} (ell - 1)/2  =  N**0.25 / 2**pi(y),

    so the ``p <-> q`` symmetry is worth a factor ``2**pi(y)``.  Since
    ``pi(y) = O(log N / log log N)`` that is ``N**o(1)`` -- a subexponential
    saving, never a polynomial one.
    """
    from math import log

    from .arith import sieve

    target = n_bits / 4
    primes, acc = [], 0.0
    for ell in sieve(prime_bound):
        if acc >= target:
            break
        primes.append(ell)
        acc += log(ell, 2)
    work = sum(log((ell - 1) / 2, 2) for ell in primes if ell > 2)
    return {
        "n_bits": n_bits,
        "primes": len(primes),
        "largest_prime": primes[-1] if primes else 0,
        "log2_search": work,
        "log2_target": target,
        "saving_bits": target - work,
    }
