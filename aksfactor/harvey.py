"""Harvey's exponent-one-fifth deterministic factoring algorithm.

Implements the algorithm of

    David Harvey, *An exponent one-fifth algorithm for deterministic integer
    factorisation*, arXiv:2010.05450 (Math. Comp. 90 (2021)),

which pushes Hittmeir's ``N**(2/9)`` bound down to ``N**(1/5+o(1))``, improving
on the ``N**(1/4)`` of Pollard-Strassen that had stood since the 1970s.

The three ingredients:

**Lehman's strategy.**  For a semiprime ``N = pq`` there are small ``a, b`` with
``aq + bp`` lying in a short, explicitly known interval (Lemma 3.3 of the paper).
Knowing ``u = aq + bp`` recovers ``p`` and ``q``: they come from the roots of
``y**2 - u*y + abN`` (Lemma 3.1).

**Hittmeir's congruence.**  Fermat gives, for any ``alpha`` coprime to ``N``,

    alpha**(a*q + b*p)  ==  alpha**(a*N + b)   (mod p),

so a candidate ``u`` can be tested *modulo p* without knowing ``p`` -- by taking
a gcd.  The right-hand side is computable; the left-hand side is what we search.

**A single baby-step/giant-step sweep.**  Writing the unknown offset as
``i + j*m`` turns the search into matching ``alpha**(-j*m) * t_{a,b}`` against
the table ``{alpha**i}``.  Hittmeir applies BSGS to chunks of the search space;
Harvey's improvement is to sweep the *whole* space with one search, which is
where the exponential speedup comes from.

Unmatched candidates are then swept for collisions modulo ``p`` or ``q`` alone
(Algorithm 4.1), using the product-tree and multipoint-evaluation machinery
already in :mod:`aksfactor.fast`.

Caveats, stated plainly.  The paper uses Bluestein's algorithm for the
multipoint evaluation at powers of ``alpha``; this implementation uses the
generic evaluator, costing an extra log factor.  And the crossover against
Strassen is far beyond anything reachable in Python -- see
``experiments/results/exp16_harvey.md``.  What is offered here is a correct,
tested implementation, not a fast one.
"""

from __future__ import annotations

from math import gcd, isqrt

from .arith import is_prime
from .fast import fast_spf, multipoint_eval, product_of

__all__ = [
    "lehman_recover",
    "collision_search",
    "harvey_search",
    "harvey_factor",
    "divisor_summatory",
    "candidate_counts",
    "exponent_runs",
    "run_length_bound",
    "lehman_interval",
    "lehman_intervals",
    "minimum_subcover",
]


def _isqrt_ceil(x: int) -> int:
    r = isqrt(x)
    return r if r * r == x else r + 1


def lehman_recover(n: int, a: int, b: int, u: int):
    """Lemma 3.1: test whether ``u == a*q + b*p`` and if so return ``(p, q)``.

    If ``N = pq`` and ``u = aq + bp`` then ``aq`` and ``bp`` are the roots of
    ``y**2 - u*y + abN``, so the discriminant ``u**2 - 4abN`` must be a perfect
    square.
    """
    disc = u * u - 4 * a * b * n
    if disc < 0:
        return None
    s = isqrt(disc)
    if s * s != disc:
        return None
    if (u + s) % 2:
        return None
    hi, lo = (u + s) // 2, (u - s) // 2
    for y_a, y_b in ((hi, lo), (lo, hi)):
        if y_a <= 0 or y_b <= 0 or y_a % a or y_b % b:
            continue
        q, p = y_a // a, y_b // b
        if p > 1 and q > 1 and p * q == n:
            return (p, q) if p <= q else (q, p)
    return None


def _split(n: int, g: int):
    """Normalise a divisor into an ordered pair ``(p, q)`` with ``p <= q``."""
    other = n // g
    return (g, other) if g <= other else (other, g)


def collision_search(n: int, alpha_powers: list[int], values: list[int]):
    """Algorithm 4.1: find a collision ``v_h == alpha**i`` modulo ``p`` or ``q``.

    Builds ``f(x) = prod (x - v_h)`` and evaluates it at every ``alpha**i``.  A
    non-trivial ``gcd(N, f(alpha**i))`` means some ``v_h - alpha**i`` shares a
    factor with ``N`` -- a collision mod one prime but not the other.
    """
    if not values:
        return None
    f = product_of([[(-v) % n, 1 % n] for v in values], n)
    for i, val in enumerate(multipoint_eval(f, alpha_powers, n)):
        g = gcd(n, val)
        if 1 < g < n:
            return _split(n, g)
        if g == n:
            ai = alpha_powers[i]
            for v in values:
                gg = gcd(n, (v - ai) % n)
                if 1 < gg < n:
                    return _split(n, gg)
    return None


def harvey_search(n: int, r: int, m: int, alpha: int = 2, stats: dict | None = None):
    """Algorithm 4.2: the main search.

    Assumes ``n`` is prime or a semiprime ``pq`` with ``(n/r)**0.5 <= p < n**0.5``.
    Returns ``(p, q)``, or ``None`` meaning "no factors found" (which, under the
    stated precondition, certifies primality).

    ``None`` is also returned if ``alpha`` turns out to have order below ``m``
    modulo ``n``; :func:`harvey_factor` retries with another base.
    """
    if gcd(alpha, n) != 1:
        g = gcd(alpha, n)
        return _split(n, g) if 1 < g < n else None

    # Step 1: the table of powers, which simultaneously certifies ord >= m.
    powers = [0] * m
    x = 1 % n
    for i in range(m):
        powers[i] = x
        if i >= 1:
            g = gcd(n, x - 1)
            if 1 < g < n:
                return _split(n, g)
            if g == n:
                return None  # ord_n(alpha) < m; caller should pick another alpha
        x = x * alpha % n
    index = {v: i for i, v in enumerate(powers)}
    if len(index) < m:
        return None  # powers not distinct: precondition violated

    inv_step = pow(alpha, -m, n)
    root_n = isqrt(n)
    unmatched: list[int] = []
    triples = 0

    # Step 2 and 3: sweep (a, b, j), matching against the table as we go.
    for a in range(1, r + 1):
        for b in range(1, r // a + 1):
            ab = a * b
            root = _isqrt_ceil(4 * ab * n)
            t = pow(alpha, a * n + b - root, n)
            j_bound = root_n // (4 * r * m * isqrt(ab) if isqrt(ab) else 1)
            v = t
            for j in range(0, max(j_bound, 0) + 1):
                triples += 1
                i = index.get(v)
                if i is not None:
                    got = lehman_recover(n, a, b, i + j * m + root)
                    if got:
                        if stats is not None:
                            stats["triples"] = triples
                        return got
                else:
                    unmatched.append(v)
                v = v * inv_step % n

    if stats is not None:
        stats.update({"triples": triples, "unmatched": len(unmatched), "m": m, "r": r})

    # Step 4: collisions modulo a single prime.
    return collision_search(n, powers, unmatched)


def harvey_factor(n: int, r: int | None = None, m: int | None = None,
                  stats: dict | None = None):
    """Algorithm 4.3: factor an ``n`` that is prime or semiprime.

    Returns ``(p, q)`` or ``None`` if ``n`` is prime.  ``r`` and ``m`` default to
    the paper's asymptotic choices, ``r ~ N**0.2 / lg(N)**0.8`` and
    ``m ~ N**0.2 * lg(N)**1.2``, clamped to stay sane on small inputs.
    """
    if n < 2:
        raise ValueError("harvey_factor expects n >= 2")
    if n % 2 == 0:
        return _split(n, 2)
    if is_prime(n):
        return None
    root = isqrt(n)
    if root * root == n:
        return (root, root)

    lg = max(n.bit_length(), 2)
    if r is None:
        r = max(1, round(n ** 0.2 / lg**0.8))
    if m is None:
        m = max(2, round(n ** 0.2 * lg**1.2))
    m = min(m, root)  # never worth tabulating more than sqrt(n) powers

    # Step 2: clear any prime factor below (n/r)^(1/2) with Strassen.
    bound = _isqrt_ceil(n // r) if r > 1 else root
    small = fast_spf(n, bound=bound)
    if small is not None and 1 < small < n:
        return _split(n, small)

    # Steps 3 and 4: certify a base of large order, then sweep.
    for alpha in (2, 3, 5, 7, 11, 13, 17, 19, 23):
        if gcd(alpha, n) != 1:
            return _split(n, gcd(alpha, n))
        got = harvey_search(n, r, m, alpha, stats=stats)
        if got:
            return got
    return None


# --------------------------------------------------------------------------
# Anatomy of Remark 3.4 -- where the N**(1/6) question binds
# --------------------------------------------------------------------------

def divisor_summatory(r: int) -> int:
    """``#{(a, b) : a*b <= r}``, by the hyperbola method in ``O(sqrt(r))``.

    Asymptotically ``r (ln r + 2*gamma - 1)``; this is the exact count, and it
    is the source of the ``+r`` term in Harvey's candidate bound.
    """
    if r < 1:
        return 0
    total, k = 0, isqrt(r)
    for a in range(1, k + 1):
        total += r // a
    return 2 * total - k * k


def candidate_counts(n: int, r: int, m: int) -> dict:
    """Split Harvey's candidate count into its two terms.

    ``s ~ (sqrt(N) ln r)/(2 sqrt(r) m)  +  r ln r``: the first term counts the
    ``j`` sweep for each ``(a, b)``, the second counts one candidate per pair.
    At the optimum ``r = m = N**(1/5)`` the ratio of the second to the first is
    exactly ``2``, independent of ``N`` -- so two thirds of the work is the pair
    enumeration, which is what pins the exponent at ``1/5``.
    """
    from math import log, sqrt

    nf = float(n)
    ln_r = log(r) if r > 1 else 1.0
    j_terms = (sqrt(nf) * ln_r) / (2 * sqrt(r) * m) if r and m else 0.0
    pair_terms = r * ln_r
    return {
        "r": r,
        "m": m,
        "j_candidates": j_terms,
        "pair_candidates": pair_terms,
        "pair_share": pair_terms / (j_terms + pair_terms) if j_terms + pair_terms else 0.0,
        "exact_pairs": divisor_summatory(r) if r <= 10**7 else None,
    }


def run_length_bound(b: int, n: int) -> float:
    """Longest run of constant first difference in ``e(1, b)`` near ``b``.

    ``e(1,b) = N + b - ceil(2 sqrt(bN))``.  Its first difference stays constant
    only while ``f'(b) = sqrt(N/b)`` moves by less than one, and
    ``f''(b) = -sqrt(N)/(2 b**1.5)``, giving a run length of at most about
    ``2 b**1.5 / sqrt(N)``.

    With Lehman's ``r = N**(1/3)`` and ``b <= r`` this is at most about ``2``:
    the sequence is nowhere locally an arithmetic progression, so there is no
    geometric structure for a baby-step/giant-step sweep to exploit.
    """
    from math import sqrt

    return 2 * b**1.5 / sqrt(float(n))


def exponent_runs(n: int, r: int) -> dict:
    """Measure runs of constant first difference of ``e(1, b)`` for ``b <= r``."""
    exps = [n + b - _isqrt_ceil(4 * b * n) for b in range(1, r + 1)]
    diffs = [exps[i + 1] - exps[i] for i in range(len(exps) - 1)]
    if not diffs:
        return {"r": r, "longest": 0, "mean": 0.0, "distinct": 0}
    runs, cur = [], 1
    for i in range(1, len(diffs)):
        if diffs[i] == diffs[i - 1]:
            cur += 1
        else:
            runs.append(cur)
            cur = 1
    runs.append(cur)
    return {
        "r": r,
        "longest": max(runs),
        "mean": sum(runs) / len(runs),
        "distinct": len(set(diffs)),
        "predicted_bound": run_length_bound(r, n),
    }


def lehman_interval(n: float, r: int, a: int, b: int):
    """The exact set of ``p`` that the pair ``(a, b)`` certifies, as an interval.

    Lehman's condition is ``0 <= a*N/p + b*p - 2 sqrt(abN) < W`` with
    ``W = sqrt(N)/(4 r sqrt(ab))``.  Multiplying through by ``p > 0`` turns it
    into the quadratic ``b p^2 - (2 sqrt(abN) + W) p + aN < 0``, so the certified
    ``p`` lie strictly between its roots.

    The interval is centred at ``p* = sqrt(aN/b)`` -- where the AM-GM bound
    ``aq + bp >= 2 sqrt(abN)`` is tight -- and its half-width is approximately
    ``sqrt(N)/(2 b sqrt(r))``, which depends only on ``b``.

    Returns ``(lo, hi)`` or ``None`` if the pair certifies nothing.
    """
    from math import sqrt

    ab = a * b
    w = sqrt(n) / (4 * r * sqrt(ab))
    top = 2 * sqrt(ab * n) + w
    disc = top * top - 4 * b * (a * n)
    if disc <= 0:
        return None
    root = sqrt(disc)
    return ((top - root) / (2 * b), (top + root) / (2 * b))


def lehman_intervals(n: float, r: int):
    """Every pair's certified interval, clipped to ``[sqrt(N/r), sqrt(N))``.

    Returns ``(intervals, lo, hi)`` where each interval is ``(lo, hi, a, b)``.
    By Lemma 3.3 the union must cover ``[lo, hi)`` -- which
    ``tests/test_harvey.py`` checks directly, a computational verification of
    Lehman's theorem.
    """
    from math import sqrt

    lo, hi = sqrt(n / r), sqrt(n)
    out = []
    for a in range(1, r + 1):
        for b in range(1, r // a + 1):
            iv = lehman_interval(n, r, a, b)
            if iv is None:
                continue
            left, right = max(iv[0], lo), min(iv[1], hi)
            if right > left:
                out.append((left, right, a, b))
    return out, lo, hi


def minimum_subcover(intervals, lo: float, hi: float, eps: float = 1e-9):
    """Fewest intervals needed to cover ``[lo, hi]`` (greedy, which is optimal).

    Returns ``None`` if the intervals do not cover the range at all.

    The point of measuring this: an ``N**(1/6)`` algorithm would need the
    ``Theta(r log r)`` Lehman pairs replaced by ``O(sqrt(r))`` of them.  Measured,
    the minimum subcover is a **constant fraction** (about 0.70) of the pairs, and
    that fraction is flat across a 64-fold range of ``r`` -- so the covering is
    essentially non-redundant and no such reduction exists.
    """
    ordered = sorted(intervals)
    total, i, cur, used = len(ordered), 0, lo, 0
    while cur < hi - eps:
        best = None
        while i < total and ordered[i][0] <= cur + eps:
            if best is None or ordered[i][1] > best:
                best = ordered[i][1]
            i += 1
        if best is None or best <= cur + eps:
            return None
        cur = best
        used += 1
    return used
