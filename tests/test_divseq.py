import random

from aksfactor.arith import is_prime, sieve
from aksfactor.divseq import (
    curve_eds,
    eds_at,
    eds_mod,
    elliptic_triangle_factor,
    fibonacci,
    mixed_radix_carries,
    naturals,
    nomial,
    point_order,
    q_integers,
    rank_of_apparition,
)

CURVES = [(0, 17, 2, 5), (0, -2, 3, 5), (1, 1, 0, 1), (-7, 10, 1, 2)]


def test_generalised_binomials_are_integers():
    seqs = [naturals(40), q_integers(3, 40), fibonacci(40)] + [curve_eds(*c, 30) for c in CURVES]
    for seq in seqs:
        rows = min(30, len(seq) - 1)
        for n in range(rows):
            for k in range(n + 1):
                nomial(seq, n, k)          # raises if not an integer


def test_carry_rule_and_rank_is_order():
    for c in CURVES:
        seq = curve_eds(*c, 50)
        for p in (7, 11, 13, 19, 23):
            if seq[2] % p == 0 or (4 * c[0] ** 3 + 27 * c[1] ** 2) % p == 0:
                continue
            r = rank_of_apparition(seq, p)
            if r is None or r > 30:
                continue
            assert r == point_order(*c, p)
            for n in range(35):
                for k in range(n + 1):
                    zero = nomial(seq, n, k) % p == 0
                    assert zero == (mixed_radix_carries(n, k, r, p) > 0)


def test_double_and_add_matches_recurrence():
    w = curve_eds(0, 17, 2, 5, 4)
    m = 1000003 * 999983
    ref = eds_mod(*w[1:5], 600, m)
    assert all(eds_at(*w[1:5], t, m) == ref[t] for t in range(1, 601))


def test_elliptic_triangle_factors():
    rng = random.Random(6)
    while True:
        p = rng.randrange(1 << 19, 1 << 20) | 1
        if is_prime(p):
            break
    q = 1000000007
    assert elliptic_triangle_factor(p * q, 1000, 400, rng) in (p, q)


def test_zero_pattern_matches_exact_binomials():
    from math import comb

    from aksfactor.figures import zero_pattern

    nat = list(range(121))
    for p in (5, 7, 11):
        z = zero_pattern(nat, p, 4, 120)
        assert all((comb(n, k) % p == 0) == z[n][k] for n in range(120) for k in range(n + 1))
    w = curve_eds(0, 17, 2, 5, 40)
    z = zero_pattern(eds_mod(*w[1:5], 41, 7 ** 5), 7, 5, 38)
    assert all((nomial(w, n, k) % 7 == 0) == z[n][k] for n in range(38) for k in range(n + 1))


def test_point_count_factors():
    from math import gcd

    from aksfactor.divseq import factor_with_point_count

    def count(a, b, p):
        total = 1
        for x in range(p):
            r = (x ** 3 + a * x + b) % p
            total += 1 if r == 0 else (2 if pow(r, (p - 1) // 2, p) == 1 else 0)
        return total

    rng = random.Random(9)
    p, q = 2003, 3001
    n = p * q
    done = 0
    for _ in range(10):
        while True:
            a, x, y = rng.randrange(n), rng.randrange(n), rng.randrange(n)
            b = (y * y - x ** 3 - a * x) % n
            if gcd((4 * a ** 3 + 27 * b * b) % n, n) == 1:
                break
        g, _ = factor_with_point_count(n, a, (x, y), count(a % p, b % p, p) * count(a % q, b % q, q))
        done += g in (p, q)
    assert done >= 9
