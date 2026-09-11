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
