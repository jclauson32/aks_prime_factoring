"""The central column of Pascal's triangle computes Legendre symbols."""

import random
import sys
from math import comb, gcd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from aksfactor.arith import is_prime
from aksfactor.central import (
    central_split,
    central_sum,
    central_sum_naive,
    jacobi,
    power_series_naive,
)


def _legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def test_jacobi_matches_legendre_on_primes():
    for p in (3, 5, 7, 11, 13, 101, 1009):
        for a in range(0, 60):
            assert jacobi(a, p) == _legendre(a, p), (a, p)


def test_jacobi_is_multiplicative_in_the_modulus():
    for p, q in [(11, 13), (101, 103), (7, 1009)]:
        for a in range(1, 50):
            assert jacobi(a, p * q) == jacobi(a, p) * jacobi(a, q)


def test_naive_sum_matches_definition():
    n = 10007 * 10009
    for t in range(0, 40):
        want = sum(comb(2 * k, k) * 3**k for k in range(t + 1)) % n
        assert central_sum_naive(t, 3, n) == want


def test_fast_sum_matches_naive():
    for n in (10007 * 10009, 1009 * 1013):
        for a in (3, 5, 7):
            for t in (5, 17, 64, 100, 257, 1000):
                s, g = central_sum(t, a, n)
                if g is None:
                    assert s == central_sum_naive(t, a, n), (n, a, t)


def test_the_legendre_identity():
    """P_T(a) = Legendre((1-4a)/p) mod p for every T in [(p-1)/2, p-1]."""
    for p in (11, 13, 17, 19, 23, 29, 31, 37, 41, 43):
        for a in range(1, 20):
            leg = _legendre(1 - 4 * a, p)
            if leg == 0:
                continue
            for t in range((p - 1) // 2, p):
                assert central_sum_naive(t, a, p) == leg % p, (p, a, t)


def test_splits_inside_the_window_never_below():
    """The guard against round 1's trap: with p large, below-window T must fail."""
    rng = random.Random(8)

    def randprime(lo, hi):
        while True:
            x = rng.randrange(lo, hi) | 1
            if is_prime(x):
                return x

    for _ in range(4):
        p = randprime(1500, 4000)
        q = randprime(p + 2, int(1.6 * p))
        n = p * q
        a = next(a for a in range(2, 500) if jacobi(1 - 4 * a, n) == -1)
        lo = (p - 1) // 2
        inside = central_sum_naive(lo + (p - 1 - lo) // 2, a, n)
        assert any(1 < gcd((inside - e) % n, n) < n for e in (1, n - 1))
        for t in (10, lo // 4, lo // 2, lo - 5):
            below = central_sum_naive(t, a, n)
            assert not any(1 < gcd((below - e) % n, n) < n for e in (1, n - 1)), t


def test_central_split_uses_the_legendre_mechanism_below_p():
    """Splits are always correct; the Legendre mechanism accounts for almost all.

    A split *below* the window is possible -- the partial sum can hit +-1 mod p
    by chance, a ~2/p fluke per trial -- and the classifier must say so rather
    than credit the mechanism.  An earlier version labeled every split
    "Legendre symbol"; this test caught it on p = 16249, T = 512.
    """
    rng = random.Random(21)

    def randprime(lo, hi):
        while True:
            x = rng.randrange(lo, hi) | 1
            if is_prime(x):
                return x

    legendre = 0
    trials = 8
    for _ in range(trials):
        p = randprime(10**3, 4 * 10**4)
        q = randprime(p + 2, int(1.8 * p))
        stats = {}
        assert central_split(p * q, stats=stats) == (p, q)
        if stats["mechanism"] == "Legendre symbol":
            legendre += 1
            # it fires below p, where gcd(T!, N) = 1 and Strassen sees nothing
            assert stats["T"] <= p - 1
        else:
            assert stats["mechanism"].startswith("coincidence")
    assert legendre >= trials - 2, legendre


def test_power_series_coefficients_are_integers():
    """(1 - d^d x)^(-1/d) has integer coefficients."""
    from fractions import Fraction

    for d in (2, 3, 4, 5):
        c = Fraction(1)
        for k in range(12):
            assert c.denominator == 1, (d, k, c)
            c = c * Fraction(d ** (d - 1) * (1 + d * k), k + 1)


def test_power_series_d2_is_the_central_column():
    n = 10007 * 10009
    for t in range(0, 30):
        assert power_series_naive(t, 5, n, 2) == central_sum_naive(t, 5, n)


def test_cube_character_when_three_divides_p_minus_one():
    """For p = 1 mod 3, the d = 3 series gives the cubic character at T ~ p/3."""
    for p in (7, 13, 19, 31, 37, 43, 61, 67, 73, 79):
        assert p % 3 == 1
        for a in range(1, 12):
            base = (1 - 27 * a) % p
            if base == 0:
                continue
            want = pow(base, (p - 1) // 3, p)
            for t in range((p - 1) // 3, p):
                assert power_series_naive(t, a, p, 3) == want, (p, a, t)
