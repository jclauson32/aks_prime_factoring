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
