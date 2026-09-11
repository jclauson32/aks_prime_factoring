import random
from fractions import Fraction

from aksfactor.arith import is_prime, sieve
from aksfactor.schnorr import (
    lll_integral,
    prime_lattice_relations,
    smooth_exponents,
    split_from_relations,
    support_coprime,
)


def _check_reduced(b):
    n = len(b)
    ortho, norms, mu = [], [], [[Fraction(0)] * n for _ in range(n)]
    for i in range(n):
        v = [Fraction(x) for x in b[i]]
        for j in range(i):
            mu[i][j] = sum(Fraction(x) * y for x, y in zip(b[i], ortho[j])) / norms[j]
            v = [a - mu[i][j] * c for a, c in zip(v, ortho[j])]
        ortho.append(v)
        norms.append(sum(x * x for x in v))
    size = all(abs(mu[i][j]) <= Fraction(1, 2) for i in range(n) for j in range(i))
    lovasz = all(norms[k] >= (Fraction(99, 100) - mu[k][k - 1] ** 2) * norms[k - 1]
                 for k in range(1, n))
    prod = 1
    for x in norms:
        prod *= x
    return size, lovasz, prod


def test_lll_integral_reduces_and_keeps_the_lattice():
    rng = random.Random(1)
    for n, mag in [(5, 10 ** 4), (12, 10 ** 9), (20, 2 ** 50)]:
        basis = [[rng.randrange(-mag, mag) for _ in range(n + 1)] for _ in range(n)]
        red = lll_integral(basis)
        size, lovasz, gram = _check_reduced(red)
        assert size and lovasz
        assert gram == _check_reduced(basis)[2]


def test_residue_avoids_the_primes_the_vector_used():
    rng = random.Random(2)
    primes = sieve(100)[:20]
    n = 1000003 * 999983
    count = 0
    for _ in range(10):
        for u, v, r in prime_lattice_relations(n, primes, 1.0, rng):
            assert u - v * n == r
            assert support_coprime(u, v, r, primes)
            count += 1
    assert count > 50


def test_split_from_relations_on_linear_sieve_relations():
    """Relations u = u mod n with u a random product of factor-base primes."""
    rng = random.Random(3)
    n = 1009 * 1013
    primes = sieve(60)
    rels = []
    while len(rels) < len(primes) + 8:
        u = 1
        while u < n:
            u *= rng.choice(primes)
        r = u % n
        if r > n // 2:
            r -= n
        if r and smooth_exponents(r, primes):
            rels.append((u, r))
    g = split_from_relations(n, primes, rels)
    assert g in (1009, 1013)
