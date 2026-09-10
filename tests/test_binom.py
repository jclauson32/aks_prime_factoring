"""Cross-validate the three binomial engines against math.comb and each other."""

import sys
from math import comb
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from aksfactor.arith import sieve
from aksfactor.binom import (
    binom_mod_falling,
    binom_mod_prime,
    binom_mod_prime_power,
    binom_mod_small,
)


def test_lucas_matches_comb():
    for p in sieve(30):
        for N in range(0, 90):
            for M in range(0, N + 1):
                assert binom_mod_prime(N, M, p) == comb(N, M) % p


def test_granville_matches_comb():
    prime_powers = []
    for p in sieve(20):
        e, q = 1, p
        while q <= 128:
            prime_powers.append((p, e))
            e += 1
            q = p**e
    for N in range(0, 100):
        for M in range(0, N + 1):
            c = comb(N, M)
            for p, e in prime_powers:
                assert binom_mod_prime_power(N, M, p, e) == c % p**e, (N, M, p, e)


def test_crt_composite_modulus():
    for N in range(0, 80):
        for M in range(0, N + 1):
            c = comb(N, M)
            for m in range(1, 50):
                assert binom_mod_small(N, M, m) == c % m, (N, M, m)


def test_huge_n_agrees_with_reference():
    for N in (10**40 + 7, 2**127 - 1, 987654321987654321987654321):
        for M in range(0, 40):
            for m in range(2, 40):
                assert binom_mod_small(N, M, m) == binom_mod_falling(N, M, m), (N, M, m)


def test_out_of_range():
    assert binom_mod_small(10, 11, 7) == 0
    assert binom_mod_small(10, -1, 7) == 0
    assert binom_mod_prime(5, 9, 3) == 0
