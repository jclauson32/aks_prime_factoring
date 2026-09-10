"""Binomial coefficients modulo small moduli, for astronomically large inputs.

The functions here compute ``C(N, M) mod m`` where ``N`` and ``M`` may have
thousands of digits but ``m`` is small enough to factor.  This is what makes
*random access* into a Pascal row possible: see :func:`aksfactor.pascal.row_entry`.

Three independent routes are provided, and the test-suite cross-checks them
against :func:`math.comb` exhaustively on small inputs:

``binom_mod_prime``
    Lucas' theorem.  Handles huge ``N`` and ``M`` for a prime modulus.
``binom_mod_prime_power``
    Granville's generalisation of Lucas to prime powers.  Also handles huge
    ``N`` and ``M``.
``binom_mod_falling``
    Direct falling-factorial evaluation with the ``p``-part split off.  Needs
    small ``M``, and exists as an independent reference implementation.
"""

from __future__ import annotations

from math import comb

from .arith import carries, factorize_small, valuation

__all__ = [
    "binom_mod_prime",
    "binom_mod_prime_power",
    "binom_mod_small",
    "binom_mod_falling",
]


def binom_mod_prime(N: int, M: int, p: int) -> int:
    """``C(N, M) mod p`` for prime ``p``, by Lucas' theorem."""
    if M < 0 or M > N:
        return 0
    res = 1
    while N or M:
        n_d, m_d = N % p, M % p
        if m_d > n_d:
            return 0
        res = res * comb(n_d, m_d) % p
        N //= p
        M //= p
    return res


def _pfree_prefix(p: int, e: int) -> list[int]:
    """Prefix products, mod ``p**e``, of the integers not divisible by ``p``."""
    q = p**e
    tab = [1] * (q + 1)
    cur = 1
    for i in range(1, q + 1):
        if i % p:
            cur = cur * i % q
        tab[i] = cur
    return tab


def _fact_pfree(n: int, p: int, e: int, tab: list[int]) -> int:
    """``prod{ i <= n : p does not divide i } mod p**e``.

    The p-free residues mod ``p**e`` form a group whose full-period product is
    ``-1``, except for ``p == 2, e >= 3`` where it is ``+1`` (the unit group
    stops being cyclic there).
    """
    if n <= 0:
        return 1
    q = p**e
    full, rem = divmod(n, q)
    sign = 1 if (p == 2 and e >= 3) else q - 1
    return pow(sign, full, q) * tab[rem] % q


def binom_mod_prime_power(N: int, M: int, p: int, e: int) -> int:
    """``C(N, M) mod p**e``, Granville-style.  ``N``, ``M`` may be huge.

    Rests on the exact factorial identity

        n! = p**v_p(n!) * prod_{j>=0} (floor(n/p**j)!)_p

    where ``(m!)_p`` is the product of the integers ``<= m`` coprime to ``p``.
    Dividing the identity for ``N!`` by those for ``M!`` and ``R! = (N-M)!``
    gives ``C(N, M) / p**c0`` as an exact ratio of ``p``-free products, where
    ``c0`` is the carry count of Kummer's theorem.  Both sides are then reduced
    mod ``p**e``; the denominator is invertible because it is ``p``-free.

    No sign correction is needed: the ``(+-1)`` appearing in textbook statements
    of Granville's theorem is absorbed by evaluating ``(m!)_p`` exactly, which
    :func:`_fact_pfree` does via the full-period product of the unit group.
    """
    if M < 0 or M > N:
        return 0
    q = p**e
    if q == 1:
        return 0
    R = N - M
    c0 = sum(carries(M, R, p))            # = v_p(C(N, M))   [Kummer]
    if c0 >= e:
        return 0
    tab = _pfree_prefix(p, e)

    num = den = 1
    n_i, m_i, r_i = N, M, R
    while n_i:
        num = num * _fact_pfree(n_i, p, e, tab) % q
        den = den * _fact_pfree(m_i, p, e, tab) % q
        den = den * _fact_pfree(r_i, p, e, tab) % q
        n_i //= p
        m_i //= p
        r_i //= p

    return num * pow(den, -1, q) % q * pow(p, c0, q) % q


def binom_mod_small(N: int, M: int, m: int) -> int:
    """``C(N, M) mod m`` for a small (factorable) modulus ``m``.

    ``N`` and ``M`` are unrestricted in size.  Works by CRT over the prime
    power components of ``m``.
    """
    if m == 1:
        return 0
    if M < 0 or M > N:
        return 0
    residues: list[tuple[int, int]] = []
    for p, e in factorize_small(m).items():
        q = p**e
        r = binom_mod_prime(N, M, p) if e == 1 else binom_mod_prime_power(N, M, p, e)
        residues.append((r % q, q))
    # Chinese remainder over pairwise coprime prime powers.
    x, mod = 0, 1
    for r, q in residues:
        # solve x' = x (mod mod), x' = r (mod q)
        t = (r - x) * pow(mod, -1, q) % q
        x += mod * t
        mod *= q
    return x % m


def binom_mod_falling(N: int, M: int, m: int) -> int:
    """``C(N, M) mod m`` by direct falling-factorial evaluation.

    Reference implementation, independent of Lucas/Granville.  Costs ``O(M)``
    big-integer operations, so ``M`` must be small; ``N`` may be huge.
    """
    if M < 0 or M > N:
        return 0
    if m == 1:
        return 0
    out = 0
    parts: list[tuple[int, int]] = []
    for p, e in factorize_small(m).items():
        q = p**e
        # C(N, M) = prod_{i<M} (N - i) / M!  -- split off the power of p.
        v = 0
        unit = 1
        for i in range(M):
            t = N - i
            k = valuation(t, p)
            v += k
            unit = unit * ((t // p**k) % q) % q
        den = 1
        for i in range(1, M + 1):
            k = valuation(i, p)
            v -= k
            den = den * ((i // p**k) % q) % q
        if v >= e:
            parts.append((0, q))
            continue
        val = unit * pow(den, -1, q) % q * pow(p, v, q) % q
        parts.append((val, q))
    x, mod = 0, 1
    for r, q in parts:
        t = (r - x) * pow(mod, -1, q) % q
        x += mod * t
        mod *= q
    out = x % m
    return out
