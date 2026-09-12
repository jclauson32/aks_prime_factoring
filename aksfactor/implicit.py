"""Round 37: an alignment *between* instances (May-Ritzenhofen).

Rounds 12 and 13 priced information about one modulus: Coppersmith finishes
from a quarter of `p`'s bits, and the bits themselves cost `N^(1/4)` to guess.
This module prices a different currency -- not bits handed to you, but bits two
moduli happen to *share*.

Let ``N_i = p_i q_i`` for ``i = 1..k``, with every ``p_i`` congruent mod
``M = 2^t`` (they share their ``t`` low bits; nobody knows what those bits are)
and every ``q_i`` about ``alpha`` bits.  Then

    q_1 N_i - q_i N_1  =  q_1 q_i (p_i - p_1)  ==  0   (mod M),

so the vector ``(q_1, ..., q_k)`` -- the thing we want -- lies in the lattice

    L = { x in Z^k : x_1 N_i == x_i N_1 (mod M) },   det L = M^(k-1).

``(q_1, ..., q_k)`` has norm about ``2^alpha sqrt k`` while the Gaussian
heuristic puts a generic shortest vector at ``det^(1/k) = 2^(t (k-1)/k)``, so
the wanted vector is the shortest one, and LLL finds it, as soon as

    t  >  alpha * k / (k - 1).

Two moduli need twice the size of the hidden cofactor in shared bits; ten need
``10/9`` of it; no number of moduli gets below ``alpha`` itself.  The secret is
never cheaper than the secret.
"""

from __future__ import annotations

from math import ceil, log


def threshold_bits(k: int, alpha: int) -> int:
    """Smallest ``t`` the heuristic expects to work: ``alpha k / (k - 1)``."""
    if k < 2:
        raise ValueError("needs at least two moduli")
    return ceil(alpha * k / (k - 1))


def shared_lsb_lattice(moduli: list[int], t: int) -> list[list[int]]:
    """A basis of ``L``; its rows span ``{x : x_1 N_i == x_i N_1 mod 2^t}``."""
    if len(moduli) < 2:
        raise ValueError("needs at least two moduli")
    m = 1 << t
    inv = pow(moduli[0] % m, -1, m)                 # the moduli are odd
    k = len(moduli)
    basis = [[1] + [moduli[i] * inv % m for i in range(1, k)]]
    for i in range(1, k):
        row = [0] * k
        row[i] = m
        basis.append(row)
    return basis


def _candidate(vector: list[int], moduli: list[int]):
    """The ``p_i`` if ``vector`` is ``+-(q_1, ..., q_k)``, else ``None``."""
    sign = next((x for x in vector if x), 0)
    if not sign:
        return None
    v = [-x for x in vector] if sign < 0 else list(vector)
    if any(x <= 1 or x >= n or n % x for x, n in zip(v, moduli)):
        return None
    return [n // x for x, n in zip(v, moduli)]


def implicit_factor(moduli: list[int], t: int):
    """Factor every ``N_i`` at once from the assumption that the ``p_i`` agree
    mod ``2^t``.  Returns the list of ``p_i``, or ``None``.

    Nothing here uses the shared bits' value -- only that they are shared.
    """
    from .schnorr import lll_integral

    for row in lll_integral(shared_lsb_lattice(moduli, t)):
        found = _candidate(row, moduli)
        if found:
            return found
    return None


def shared_low_primes(k: int, pbits: int, t: int, qbits: int, rng):
    """``k`` pairs ``(p, q)`` whose ``p`` share their ``t`` low bits.

    The shared part is drawn once and kept secret; only the moduli are needed
    to run the attack.
    """
    from .arith import is_prime

    space = 1 << max(pbits - t - 1, 0)
    expected = space / (pbits * log(2))          # primes among the candidates
    if expected < 3 * k:
        # Asking for more distinct primes than the residue class contains loops
        # forever rather than failing, which is the worse of the two outcomes.
        raise ValueError(f"{space} candidates share {t} low bits at {pbits} bits, "
                         f"about {expected:.1f} of them prime -- too few to draw "
                         f"{k} distinct ones")
    low = rng.randrange(1 << (t - 1), 1 << t) | 1
    pairs = []
    while len(pairs) < k:
        p = low + (rng.randrange(1 << (pbits - t - 1), 1 << (pbits - t)) << t)
        if not is_prime(p) or any(p == seen for seen, _ in pairs):
            continue                      # equal primes would split by gcd alone
        while True:
            q = rng.randrange(1 << (qbits - 1), 1 << qbits) | 1
            if is_prime(q) and q != p:
                break
        pairs.append((p, q))
    return pairs
