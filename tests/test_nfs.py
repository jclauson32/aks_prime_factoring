import random

from aksfactor.arith import is_prime
from aksfactor.nfs import (
    algebraic_sqrt,
    norm,
    number_field_sieve,
    poly_eval,
    poly_mulmod,
    select_polynomial,
)


def test_polynomial_selection():
    n = 1000003 * 999983
    for m, f in select_polynomial(n, 3, tries=5):
        assert poly_eval(f, m) == n and f[-1] == 1


def test_norm_is_the_field_norm_on_rationals():
    f = [5, 0, 0, 1]                      # alpha = -5^(1/3)
    for a, b in [(3, 1), (7, 2), (-4, 3)]:
        assert norm(f, a, b) == a ** 3 + 5 * b ** 3


def test_algebraic_square_root_recovers_planted_squares():
    rng = random.Random(1)
    for f in ([7558, 3942, 2, 1], [39, 0, 0, 1]):
        for bits in (30, 300, 1500):
            beta = [rng.randrange(-(1 << bits), 1 << bits) for _ in range(3)]
            gamma = poly_mulmod(beta, beta, f)
            got = algebraic_sqrt(gamma, f, max(abs(c) for c in gamma).bit_length() + 16)
            assert got in (beta, [-x for x in beta])


def test_number_field_sieve_factors():
    rng = random.Random(5)
    for bits in (40, 48):
        while True:
            p = rng.randrange(1 << (bits // 2 - 1), 1 << (bits // 2)) | 1
            q = rng.randrange(1 << (bits // 2 - 1), 1 << (bits // 2)) | 1
            if p != q and is_prime(p) and is_prime(q):
                break
        assert number_field_sieve(p * q, 3, 500, 500, 20, 4000, 1500) in (p, q)
