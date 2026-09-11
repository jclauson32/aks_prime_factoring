import random
from math import isqrt

from aksfactor.arith import is_prime, sieve
from aksfactor.hyperbola import (
    class_number,
    divisors_up_to,
    hyperbola_factor,
    quadratic_floor_sum,
    row_sum,
    row_sum_naive,
)


def test_row_sum_matches_reference():
    rng = random.Random(1)
    for trial in range(1500):
        n = trial + 2 if trial < 60 else rng.randrange(2, 10 ** 9)
        top = rng.randrange(1, isqrt(n) + 1)
        assert row_sum(n, top)[0] == row_sum_naive(n, top), (n, top)


def test_row_sum_large():
    rng = random.Random(2)
    for _ in range(5):
        n = rng.randrange(10 ** 11, 10 ** 12)
        assert row_sum(n, isqrt(n))[0] == row_sum_naive(n, isqrt(n))


def test_divisor_count_is_monotone_predicate():
    p, q = 1009, 2003
    n = p * q
    assert divisors_up_to(n, p - 1) == 1
    assert divisors_up_to(n, p) == 2
    assert divisors_up_to(n, isqrt(n - 1)) == 2


def test_hyperbola_factor():
    rng = random.Random(3)
    for bits in (24, 30, 36):
        for _ in range(4):
            while True:
                p = rng.randrange(1 << (bits // 2 - 1), 1 << (bits // 2)) | 1
                if is_prime(p):
                    break
            while True:
                q = rng.randrange(p + 2, 2 * p) | 1
                if is_prime(q):
                    break
            assert hyperbola_factor(p * q) == p


def test_quadratic_floor_sum_is_a_class_number():
    for p in sieve(1500):
        if p > 3 and p % 4 == 3:
            h = class_number(-p)
            assert quadratic_floor_sum(p) == (p - 1) * (2 * p - 1) // 6 - (p - 1 - 2 * h) // 2


def test_divisor_is_a_vertex_and_its_edges_bracket_p_over_q():
    from fractions import Fraction

    from aksfactor.hyperbola import divisor_vertex_edges

    rng = random.Random(5)
    for _ in range(10):
        while True:
            p = rng.randrange(1 << 13, 1 << 14) | 1
            if is_prime(p):
                break
        while True:
            q = rng.randrange(p + 2, 3 * p) | 1
            if is_prime(q):
                break
        got = divisor_vertex_edges(p * q)
        assert got is not None and got[0] == p
        slopes = sorted(Fraction(dy, dx) for dx, dy, _, _ in got[2])
        assert len(slopes) == 2 and slopes[0] < Fraction(p, q) < slopes[1]
        for dx, dy, k, gap in got[2]:
            assert (dy * q + dx * p) ** 2 - 4 * k * p * q == gap * gap
