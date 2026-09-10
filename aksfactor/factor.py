"""Factoring by locating the first non-zero entry of Pascal's row ``n`` mod ``n``.

Theorem 4 says that entry sits at ``k = spf(n)`` and *equals* ``n / spf(n)``, so
a single hit hands back both the smallest prime factor and its cofactor.

Two modes are offered, and the difference between them is the whole story of
this repository:

``scan``
    Walk ``k = 2, 3, 4, ...`` evaluating residues until one is non-zero.
    Faithful to the theorem, and ``Theta(spf(n))`` residue evaluations.
``certified``
    Theorem 2 proves a residue can only be non-zero at a position sharing a
    factor with ``n``, so the scan is *provably* the same search as trial
    division.  This mode does the cheap search (a 2-3-5 wheel) and then pays
    for one residue evaluation to emit the Pascal certificate.
``fast``
    The same search, batched: one gcd per block of ``c`` positions via fast
    multipoint evaluation, ``O~(n**0.25)`` ring operations instead of
    ``Theta(spf(n))``.  See :mod:`aksfactor.fast`.

``scan`` and ``certified`` are ``Theta(spf(n))`` and differ only in constant.
``fast`` is the one mode that moves the exponent.
"""

from __future__ import annotations

from math import gcd, isqrt

from .arith import is_prime, spf_trial
from .pascal import row_entry

__all__ = ["certificate", "pascal_spf", "pascal_split", "factor"]


def certificate(n: int, p: int) -> dict:
    """The Theorem 3/4 witness that ``p`` is the smallest prime factor of ``n``."""
    r = row_entry(n, p)
    return {
        "n": n,
        "p": p,
        "position": p,
        "residue": r,
        "expected_residue": n // p,
        "valid": r == n // p and n % p == 0 and is_prime(p),
        "cofactor": n // p,
        "identity": f"C({n}, {p}) mod {n} = {r} = {n}/{p}",
    }


def pascal_spf(n: int, bound: int | None = None, mode: str = "certified"):
    """Smallest prime factor of ``n`` via the first non-zero Pascal residue.

    Returns ``n`` itself when ``n`` is prime, or ``None`` when ``bound`` is set
    and no factor ``<= bound`` exists.
    """
    if n < 2:
        raise ValueError("pascal_spf expects n >= 2")
    limit = isqrt(n) if bound is None else min(bound, isqrt(n))
    if mode == "certified":
        p = spf_trial(n, bound=limit)
        if p is None:
            return None
        if p == n or p > limit:
            return n
        cert = certificate(n, p)
        if not cert["valid"]:  # pragma: no cover - would falsify Theorem 3
            raise AssertionError(f"Pascal certificate failed: {cert}")
        return p
    if mode == "fast":
        from .fast import fast_spf

        p = fast_spf(n, bound=limit)
        if p is None or p == n:
            return p
        cert = certificate(n, p)
        if not cert["valid"]:  # pragma: no cover - would falsify Theorem 3
            raise AssertionError(f"Pascal certificate failed: {cert}")
        return p
    if mode == "scan":
        for k in range(2, limit + 1):
            if row_entry(n, k):
                # Theorem 4: the first hit is the smallest prime factor.
                return k
        return None if bound is not None and bound < isqrt(n) else n
    raise ValueError(f"unknown mode {mode!r}")


def pascal_split(n: int, bound: int | None = None, mode: str = "certified"):
    """Split ``n`` as ``(p, n // p)`` using the Pascal residue.  ``None`` if unsplit."""
    if n < 2 or is_prime(n):
        return None
    p = pascal_spf(n, bound=bound, mode=mode)
    if p is None or p == n:
        return None
    return p, n // p


def factor(n: int, bound: int | None = None, mode: str = "certified") -> dict[int, int]:
    """Full factorization of ``n`` by repeated Pascal splitting.

    Cost is ``Theta(P)`` where ``P`` is the largest smallest-prime-factor met
    during the recursion, i.e. trial-division cost.  Raises ``ValueError`` if a
    cofactor resists splitting within ``bound``.
    """
    if n < 1:
        raise ValueError("factor expects a positive integer")
    out: dict[int, int] = {}
    stack = [n] if n > 1 else []
    while stack:
        m = stack.pop()
        if m == 1:
            continue
        if is_prime(m):
            out[m] = out.get(m, 0) + 1
            continue
        split = pascal_split(m, bound=bound, mode=mode)
        if split is None:
            raise ValueError(f"could not split {m} within bound {bound}")
        p, q = split
        out[p] = out.get(p, 0) + 1
        stack.append(q)
    return out
