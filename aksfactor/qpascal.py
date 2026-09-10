"""The q-deformed Pascal triangle: a period you can actually tune.

Round 3 ended on a requirement: *find a construction whose relevant parameter
varies at fixed ``p``*.  The classical row has none -- Theorem 7 pins its period
to ``p`` itself.  The Gaussian binomial coefficient supplies one.

For ``[n, k]_q``, divisibility by ``p`` is governed not by ``p`` but by

    d = ord_p(q),   the multiplicative order of the base,

and ``d`` moves as ``q`` moves.  The q-analogue of Kummer's theorem (verified
exhaustively in ``tests/test_qpascal.py``) is

    p | [n,k]_q   <=>   adding k and n-k in base d carries,
                        or  p | C(floor(n/d), floor(k/d)).

The two clauses behave very differently.  The first fires at

    k = (n mod d) + 1        (available only when  n mod d <= d - 2)

which can be *tiny* -- far below ``spf(n)``, where Theorem 2 says the classical
row is still identically zero.  That is a real deformation of the barrier.  The
second clause is a classical base-``p`` carry one level up, and fires only at
``Theta(p)`` scale, so it never helps.

It is also not an escape, and the reason is the point of the module: ``d``
divides ``p - 1``, so tuning the period is exactly the problem Pollard's ``p-1``
method already solves.  See ``docs/FINDINGS.md``, round 4.
"""

from __future__ import annotations

from math import gcd

from .binom import binom_mod_prime

__all__ = [
    "q_period",
    "q_clause_one_position",
    "q_pascal_row",
    "q_lucas_divides",
    "q_first_hit",
    "shift_split",
]


def q_period(q: int, p: int) -> int | None:
    """``ord_p(q)``, or ``None`` if ``p`` divides ``q``."""
    if q % p == 0:
        return None
    d, x = 1, q % p
    while x != 1:
        x = x * q % p
        d += 1
        if d > p:
            return None
    return d


def q_pascal_row(n: int, kmax: int, q: int, mod: int) -> list[int]:
    """``[n, j]_q mod mod`` for ``j = 0..kmax``, by the q-Pascal recurrence.

    ``[n,k]_q = [n-1,k-1]_q + q**k * [n-1,k]_q``.  Costs ``O(n * kmax)``, so this
    is for demonstration and testing on small ``n``.
    """
    prev = [1 % mod]
    for i in range(1, n + 1):
        cur = [1 % mod]
        for j in range(1, min(i, kmax) + 1):
            a = prev[j - 1] if j - 1 < len(prev) else 0
            b = prev[j] if j < len(prev) else 0
            cur.append((a + pow(q, j, mod) * b) % mod)
        prev = cur
    return prev


def q_lucas_divides(n: int, k: int, d: int, p: int) -> bool:
    """Does ``p`` divide ``[n,k]_q``?  ``d`` must be ``ord_p(q)``.

    The q-analogue of Kummer plus q-Lucas: a carry when adding ``k`` and ``n-k``
    in base ``d``, or a classical vanishing one level up.
    """
    if n // d - (n - k) // d > k // d:
        return True
    return binom_mod_prime(n // d, k // d, p) == 0


def q_clause_one_position(n: int, d: int) -> int | None:
    """Smallest ``k`` where clause one of :func:`q_lucas_divides` fires, if any.

    ``None`` when ``n mod d == d - 1``, which blocks the low-digit carry
    entirely.
    """
    low = n % d
    return low + 1 if low <= d - 2 else None


def q_first_hit(n: int, kmax: int, q: int) -> tuple[int, int] | None:
    """First ``k <= kmax`` where ``gcd([n,k]_q mod n, n)`` is a proper divisor.

    Uses the explicit row, so ``n`` must be small.  :func:`shift_split` is the
    same search done in a way that scales.
    """
    row = q_pascal_row(n, kmax, q, n)
    for k in range(1, len(row)):
        g = gcd(row[k], n)
        if 1 < g < n:
            return k, g
    return None


def shift_split(n: int, q: int, kmax: int) -> tuple[int, int] | None:
    """The scalable form of the *first clause*: ``gcd(q**(n-j) - 1, n)``.

    Clause one fires at ``k = (n mod d) + 1``, and ``d | n - (k-1)`` exactly
    there.  So hunting first-clause hits is hunting a ``j <= kmax`` with
    ``ord_p(q) | n - j``, which this does with one modular exponentiation and
    ``kmax`` multiplications -- no row required, any size ``n``.

    It does *not* reproduce every q-row hit: when ``n mod d = d - 1`` clause one
    is unavailable and the row instead fires via clause two, a classical base-p
    carry at ``Theta(p)`` scale.  Those hits are not worth chasing.

    Returns ``(j, factor)`` or ``None``.
    """
    if gcd(q, n) != 1:
        g = gcd(q, n)
        return (0, g) if 1 < g < n else None
    x = pow(q, n, n)
    qinv = pow(q, -1, n)
    for j in range(kmax + 1):
        g = gcd(x - 1, n)
        if 1 < g < n:
            return j, g
        x = x * qinv % n
    return None
