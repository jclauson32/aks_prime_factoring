"""Round 16: hunting for anything that sees ``p`` rather than ``N``.

Round 13 showed that ``N`` pins ``s = (p + q) mod ell`` to the set
``{a + N/a}`` -- about ``(ell-1)/2`` values -- and that one exact residue per
small prime is all Coppersmith needs.  This module measures three candidate
sources of that residue:

* **cheap functions of N** (``feature_accuracy``): a held-out lookup predictor
  for ``s`` given ``(N mod ell, f(N))``, scored against the *Bayes* baseline
  ``2/(ell-1)`` -- the best any predictor can do from ``N mod ell`` alone;
* **Ramanujan's tau** (``tau_table``, ``tau_sum_candidates``): ``tau(N) mod 691``
  would pin ``s mod 691`` exactly, but every known route to ``tau(N)`` passes
  through the divisors of ``N``;
* **four-square counts** (``r4``): ``r_4(N) = 8 sigma(N)`` for odd ``N``, so an
  exact lattice-point count on the 3-sphere of radius ``sqrt N`` gives ``p + q``.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from math import isqrt


def bayes_sum_accuracy(ell: int) -> float:
    """Best possible accuracy guessing ``(p+q) mod ell`` from ``N mod ell`` alone.

    For ``r = N mod ell`` and ``p`` uniform in ``F_ell*``, the sum ``s = p + r/p``
    has ``1 + ((s^2 - 4r)/ell)`` preimages, so the most likely value has two
    and the Bayes predictor is right with probability ``2/(ell-1)`` -- not the
    ``1/#{allowed sums}`` a uniform guess over the allowed set achieves.
    """
    if ell < 5:
        raise ValueError("ell must be a prime >= 5")
    return 2 / (ell - 1)


def feature_accuracy(data, ell: int, feature) -> float:
    """Held-out accuracy of a lookup table predicting ``(p+q) mod ell``.

    ``data`` is a list of ``(N, p, q)``; the first half trains a table keyed by
    ``(N mod ell, feature(N, ell))`` and the second half is scored.  Keys unseen
    in training count as misses, so the estimate is conservative.
    """
    half = len(data) // 2
    table = defaultdict(Counter)
    for n, p, q in data[:half]:
        table[(n % ell, feature(n, ell))][(p + q) % ell] += 1
    hits = 0
    for n, p, q in data[half:]:
        seen = table.get((n % ell, feature(n, ell)))
        if seen:
            hits += seen.most_common(1)[0][0] == (p + q) % ell
    return hits / (len(data) - half)


def _series_mul(a, b, size):
    out = [0] * size
    for i, x in enumerate(a):
        if x:
            for j in range(0, size - i):
                y = b[j]
                if y:
                    out[i + j] += x * y
    return out


def tau_table(bound: int) -> list[int]:
    """``tau(0..bound-1)`` from ``Delta = q * (eta^3)^8`` (``tau(0) = 0``).

    Jacobi's identity ``prod (1-q^n)^3 = sum (-1)^k (2k+1) q^(k(k+1)/2)`` makes
    ``eta^3`` sparse; three squarings give ``prod (1-q^n)^24``.
    """
    size = max(bound - 1, 1)
    e3 = [0] * size
    k = 0
    while k * (k + 1) // 2 < size:
        e3[k * (k + 1) // 2] = (-1) ** k * (2 * k + 1)
        k += 1
    s = e3
    for _ in range(3):
        s = _series_mul(s, s, size)
    return [0] + s[: bound - 1]


def tau_prime_mod_691(p: int) -> int:
    """``tau(p) mod 691`` for a prime ``p``: Ramanujan's ``1 + p^11``."""
    return (1 + pow(p, 11, 691)) % 691


def tau_sum_candidates(n: int, tau_mod: int, ell: int = 691) -> set[int]:
    """Values of ``(p+q) mod 691`` consistent with ``N`` and ``tau(N) mod 691``.

    ``tau(N) = tau(p) tau(q) = (1 + p^11)(1 + q^11) (mod 691)``.  With
    ``p^11 q^11 = N^11`` known, this fixes the pair ``{p^11, q^11}``, and since
    ``gcd(11, 690) = 1`` the eleventh-power map is a bijection on ``F_691`` --
    so the pair ``{p, q} mod 691`` and hence ``p + q mod 691`` are pinned.
    """
    if ell != 691:
        raise ValueError("Ramanujan's congruence is mod 691")
    e = n % ell
    out = set()
    for a in range(1, ell):
        b = e * pow(a, -1, ell) % ell
        if (1 + pow(a, 11, ell)) * (1 + pow(b, 11, ell)) % ell == tau_mod % ell:
            out.add((a + b) % ell)
    return out


def trace_divisor_term(n: int, weight: int = 12) -> int:
    """The divisor sum ``sum_{d | N} min(d, N/d)^(k-1)`` of the Eichler-Selberg
    trace formula.  For ``N = pq`` with ``p < q`` it is ``2 + 2 p^11``: the
    formula for ``tau(N)`` names the smaller factor explicitly."""
    return sum(min(d, n // d) ** (weight - 1) for d in range(1, n + 1) if n % d == 0)


def r4(n: int) -> int:
    """Number of ``(a, b, c, d) in Z^4`` with ``a^2 + b^2 + c^2 + d^2 = n``.

    Brute force, ``O(n^(3/2))``: a check on Jacobi's theorem, not a method.
    """
    count = 0
    root = isqrt(n)
    for a in range(-root, root + 1):
        ra = n - a * a
        for b in range(-isqrt(ra), isqrt(ra) + 1):
            rb = ra - b * b
            for c in range(-isqrt(rb), isqrt(rb) + 1):
                rc = rb - c * c
                d = isqrt(rc)
                if d * d == rc:
                    count += 1 if d == 0 else 2
    return count


def walsh_scan(train, test, target, bits: int, triple_top: int = 0,
               skip: tuple = (0,), keep_sigma: float = 4.5, max_degree: int = 2) -> dict:
    """Low-degree Fourier learner for one target bit of the factorisation.

    ``train``/``test`` are lists of ``(N, p, q)``; ``target(N, p, q)`` is a bool.
    Every Walsh character ``chi_S(N) = (-1)^(sum_{i in S} bit_i N)`` with
    ``|S| <= max_degree`` (and ``|S| = 3`` among the top ``triple_top`` bits) is
    correlated with the target on the training half, using bitsets over the
    samples (one AND and one popcount per character).  Characters beyond
    ``keep_sigma`` standard errors vote on the test half.

    ``skip`` drops constant bits -- bit 0 of an odd ``N`` would otherwise
    report the target's bias as a "correlation".
    Returns test accuracy, the majority baseline, the largest ``|corr|``, and
    the largest ``|corr|`` expected under the null, ``sigma sqrt(2 ln 2M)``.
    """
    from itertools import combinations
    from math import log, sqrt

    idx = [i for i in range(bits) if i not in skip]
    subsets = []
    for deg in range(1, max_degree + 1):
        subsets += list(combinations(idx, deg))
    if triple_top and max_degree < 3:
        subsets += list(combinations([i for i in idx if i >= bits - triple_top], 3))
    cols = []
    for i in range(bits):
        v = 0
        for s, (n, _, _) in enumerate(train):
            if (n >> i) & 1:
                v |= 1 << s
        cols.append(v)
    tv = 0
    for s, row in enumerate(train):
        if target(*row):
            tv |= 1 << s
    m = len(train)
    t_cnt = tv.bit_count()
    corr = {}
    for sub in subsets:
        par = 0
        for i in sub:
            par ^= cols[i]
        corr[sub] = (2 * t_cnt - m - 4 * (par & tv).bit_count() + 2 * par.bit_count()) / m
    sigma = 1 / sqrt(m)
    keep = {sub: c for sub, c in corr.items() if abs(c) > keep_sigma * sigma}
    ybar = (2 * t_cnt - m) / m
    hits = ones = 0
    for n, p, q in test:
        y = 1 if target(n, p, q) else -1
        ones += y == 1
        score = ybar
        for sub, c in keep.items():
            par = 0
            for i in sub:
                par ^= (n >> i) & 1
            score += c * (1 - 2 * par)
        hits += (1 if score >= 0 else -1) == y
    best = max(corr.items(), key=lambda kv: abs(kv[1]))
    return {
        "accuracy": hits / len(test),
        "majority": max(ones, len(test) - ones) / len(test),
        "max_corr": abs(best[1]),
        "argmax": best[0],
        "null_max": sigma * sqrt(2 * log(2 * len(subsets))),
        "kept": len(keep),
        "characters": len(subsets),
    }
