"""Machine-checkable statements of the results in ``docs/THEORY.md``.

Every function takes an integer (and sometimes a scan bound) and returns
``(ok, detail)``.  ``detail`` carries a counterexample when ``ok`` is False, so
the experiment scripts and the test-suite can share one source of truth.
"""

from __future__ import annotations

from math import comb, gcd

from .arith import factorize_small, is_prime, spf_trial
from .pascal import row_entry, row_prefix

__all__ = [
    "check_t1_shape",
    "check_t2_coprime_vanishing",
    "check_t3_exact_value",
    "check_t4_first_nonzero",
    "check_t5_closed_form",
    "check_gist",
    "check_all",
    "check_t7_aliasing",
    "check_t9_second_digit",
    "check_t15_no_row_leaks",
    "check_t16_q_kummer",
    "check_t20_factorial_threshold",
]


def check_t1_shape(n: int, kmax: int | None = None) -> tuple[bool, dict]:
    """T1: ``n / gcd(n,k)`` divides ``C(n,k)``; hence every residue is ``x * y``
    with ``x`` a proper divisor of ``n``."""
    kmax = n - 1 if kmax is None else min(kmax, n - 1)
    row = row_prefix(n, kmax)
    for k in range(1, kmax + 1):
        x = n // gcd(n, k)
        if row[k] % x:
            return False, {"n": n, "k": k, "residue": row[k], "x": x}
    return True, {"n": n, "checked": kmax}


def check_t2_coprime_vanishing(n: int, kmax: int | None = None) -> tuple[bool, dict]:
    """T2: ``gcd(k, n) == 1`` with ``0 < k < n`` forces ``C(n,k) == 0 mod n``."""
    kmax = n - 1 if kmax is None else min(kmax, n - 1)
    row = row_prefix(n, kmax)
    for k in range(1, kmax + 1):
        if gcd(k, n) == 1 and row[k]:
            return False, {"n": n, "k": k, "residue": row[k]}
    return True, {"n": n, "checked": kmax}


def check_t3_exact_value(n: int) -> tuple[bool, dict]:
    """T3: for every prime ``p | n``, ``C(n, p) mod n == n // p`` exactly."""
    if is_prime(n):
        return True, {"n": n, "prime": True}
    for p in factorize_small(n):
        got = row_entry(n, p)
        if got != n // p:
            return False, {"n": n, "p": p, "got": got, "want": n // p}
    return True, {"n": n, "primes": sorted(factorize_small(n))}


def check_t4_first_nonzero(n: int, kmax: int | None = None) -> tuple[bool, dict]:
    """T4: the first non-zero interior residue sits at ``k = spf(n)`` and equals
    ``n / spf(n)``; for prime ``n`` the whole interior vanishes."""
    kmax = n - 1 if kmax is None else min(kmax, n - 1)
    row = row_prefix(n, kmax)
    hit = next(((k, row[k]) for k in range(1, kmax + 1) if row[k]), None)
    if is_prime(n):
        return (hit is None), {"n": n, "prime": True, "hit": hit}
    p = spf_trial(n)
    if hit is None:
        return False, {"n": n, "spf": p, "hit": None}
    k, r = hit
    ok = (k == p) and (r == n // p)
    return ok, {"n": n, "spf": p, "k": k, "residue": r, "want": n // p}


def check_t5_closed_form(n: int, kmax: int | None = None) -> tuple[bool, dict]:
    """T5: ``C(n,k) mod n == (n/d) * (C(n-1,k-1) / (k/d) mod d)``, ``d = gcd(n,k)``."""
    kmax = n - 1 if kmax is None else min(kmax, n - 1)
    for k in range(1, kmax + 1):
        d = gcd(n, k)
        kp = k // d
        big = comb(n - 1, k - 1)
        if big % kp:
            return False, {"n": n, "k": k, "reason": "k/d does not divide C(n-1,k-1)"}
        want = (n // d) * ((big // kp) % d)
        if comb(n, k) % n != want:
            return False, {"n": n, "k": k, "got": comb(n, k) % n, "want": want}
    return True, {"n": n, "checked": kmax}


def check_gist(n: int, kmax: int | None = None) -> tuple[bool, dict]:
    """The original conjecture, stated precisely.

    "Every non-zero remainder splits as ``x * y`` where ``x`` is a factor of
    ``n``" -- and, stronger, ``gcd(residue, n)`` is always a proper non-trivial
    divisor of ``n``.  Also records how often ``x`` happens to be *prime*.
    """
    kmax = n - 1 if kmax is None else min(kmax, n - 1)
    row = row_prefix(n, kmax)
    nonzero = x_prime = 0
    for k in range(1, kmax + 1):
        r = row[k]
        if not r:
            continue
        nonzero += 1
        x = n // gcd(n, k)
        if r % x or not (1 < x < n):
            return False, {"n": n, "k": k, "residue": r, "x": x}
        g = gcd(r, n)
        if not (1 < g < n):
            return False, {"n": n, "k": k, "residue": r, "gcd": g}
        if is_prime(x):
            x_prime += 1
    return True, {"n": n, "nonzero": nonzero, "x_prime": x_prime}


def check_all(n: int, kmax: int | None = None) -> dict[str, tuple[bool, dict]]:
    """Run every checker on ``n``."""
    return {
        "T1_shape": check_t1_shape(n, kmax),
        "T2_coprime_vanishing": check_t2_coprime_vanishing(n, kmax),
        "T3_exact_value": check_t3_exact_value(n),
        "T4_first_nonzero": check_t4_first_nonzero(n, kmax),
        "T5_closed_form": check_t5_closed_form(n, kmax),
        "gist": check_gist(n, kmax),
    }


def check_t7_aliasing(n: int, r: int) -> tuple[bool, dict]:
    """T7: with ``gcd(r, n) == 1`` the row's support is equidistributed mod ``r``.

    The support of the interior of ``R(n)`` is ``{ k : gcd(k,n) > 1 }`` -- a union
    of the multiples of the prime factors.  If ``gcd(r, p) = 1`` then
    ``i |-> i*p mod r`` is a bijection of ``Z/r``, so the multiples of ``p`` land
    evenly in every class.  Folding at level ``r < p`` therefore destroys all
    *positional* information about ``p``: the classes are indistinguishable by
    occupancy, and the only remaining signal is the accidental vanishing of a
    class sum.  Positional information needs ``gcd(r, p) > 1``, i.e. ``r >= p``.
    """
    if gcd(r, n) != 1:
        return True, {"n": n, "r": r, "skipped": "gcd(r,n) > 1"}
    counts = [0] * r
    for k in range(1, n):
        if gcd(k, n) > 1:
            counts[k % r] += 1
    spread = max(counts) - min(counts)
    # each class holds the same count up to the partial final period
    return spread <= r, {"n": n, "r": r, "counts": counts, "spread": spread}


def check_t9_second_digit(n: int, kmax: int | None = None) -> tuple[bool, dict]:
    """T9: for ``gcd(k,n) == 1`` and ``C(n,k) = n*m``, ``m`` vanishes mod ``p``
    exactly when Lucas' digit condition holds.

    Compares the Lucas criterion against an exact evaluation of
    ``C(n,k) mod n**2``, and reports how often ``gcd(m, n)`` splits ``n``.
    """
    from .binom import binom_mod_prime

    kmax = n - 1 if kmax is None else min(kmax, n - 1)
    primes = sorted(factorize_small(n))
    if len(primes) < 2:
        return True, {"n": n, "skipped": "not composite with >= 2 distinct primes"}
    n2 = n * n
    splits = trials = 0
    for k in range(2, kmax + 1):
        if gcd(k, n) != 1:
            continue
        c = comb(n, k) % n2
        if c % n:
            return False, {"n": n, "k": k, "reason": "Theorem 2 violated"}
        m = c // n
        zero = [binom_mod_prime(n - 1, k - 1, p) == 0 for p in primes]
        for p, z in zip(primes, zero):
            if (m % p == 0) != z:
                return False, {"n": n, "k": k, "p": p, "lucas": z,
                               "actual": m % p == 0}
        trials += 1
        if 1 < gcd(m, n) < n:
            splits += 1
    return True, {"n": n, "trials": trials, "splits": splits,
                  "rate": (splits / trials) if trials else None}


def check_t15_no_row_leaks(n: int, samples: int = 6) -> tuple[bool, dict]:
    """T15: below ``spf(n)``, *no* row leaks anything about the factorization.

    For ``k < spf(n)`` every ``i <= k`` is invertible mod ``n``, so

        C(N, k) = prod_{i<k} (N - i) / k!   ==   depends only on (N mod n, k).

    Two rows congruent mod ``n`` therefore agree at every position below
    ``spf(n)``.  Shifting to a different row buys nothing there -- which closes
    off the whole "try another row" family in one line.
    """
    import random as _random

    s = spf_trial(n)
    if s is None or s == n:
        return True, {"n": n, "prime": True}
    rng = _random.Random(n)
    for k in range(0, min(s, 14)):
        for _ in range(samples):
            big = rng.randrange(n, 40 * n)
            same = big % n + n * rng.randrange(1, 9)
            if comb(big, k) % n != comb(same, k) % n:
                return False, {"n": n, "k": k, "N": big, "N2": same}
    return True, {"n": n, "spf": s, "checked_below": min(s, 14)}


def check_t16_q_kummer(n: int, bases=(2, 3, 5, 6, 7), kmax: int = 40) -> tuple[bool, dict]:
    """T16: the q-analogue of Kummer's theorem, checked against the explicit row.

    ``p | [n,k]_q``  iff  ``(k mod d) > (n mod d)``  or  ``p | C(n//d, k//d)``,
    where ``d = ord_p(q)``.
    """
    from .qpascal import q_lucas_divides, q_pascal_row, q_period

    checked = 0
    for p in factorize_small(n):
        for base in bases:
            if base % p == 0:
                continue
            d = q_period(base, p)
            if d is None:
                continue
            row = q_pascal_row(n, min(kmax, n), base, p)
            for k in range(1, len(row)):
                checked += 1
                if (row[k] % p == 0) != q_lucas_divides(n, k, d, p):
                    return False, {"n": n, "p": p, "q": base, "k": k, "d": d}
    return True, {"n": n, "checked": checked}


def check_t20_factorial_threshold(n: int, cap: int = 120) -> tuple[bool, dict]:
    """T20: ``gcd(i! mod n, n) == prod_p p**min(v_p(n), v_p(i!))``.

    In particular the gcd exceeds 1 exactly when ``i >= spf(n)`` -- a *monotone*
    threshold, which is what makes binary search on it legal.

    Read on the fractal picture: row ``p`` is where the mod-``p`` gasket first
    makes holes, so "which gaskets have started making holes by row ``i``" is
    "which primes divide ``i!``", and the least such row is ``spf(n)``.
    """
    def _legendre(i: int, p: int) -> int:
        v, q = 0, p
        while q <= i:
            v += i // q
            q *= p
        return v

    facs = factorize_small(n)
    spf = min(facs)
    value = 1
    for i in range(1, min(n, cap)):
        value = value * i % n
        want = 1
        for p, e in facs.items():
            want *= p ** min(e, _legendre(i, p))
        got = gcd(value, n)
        if got != want:
            return False, {"n": n, "i": i, "got": got, "want": want}
        if (got > 1) != (i >= spf):
            return False, {"n": n, "i": i, "monotone": False, "spf": spf}
    return True, {"n": n, "spf": spf, "checked": min(n, cap) - 1}
