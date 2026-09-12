import random
from math import isqrt

from aksfactor.arith import multiplicative_order
from aksfactor.ecm import is_b1b2_smooth
from aksfactor.smooth import (
    curve_order,
    dickman_rho,
    qr_table,
    smooth_rate,
    suyama_order,
    weierstrass_order,
)


def test_table_counter_agrees_with_the_exponentiating_one():
    rng = random.Random(3)
    for p in (101, 1009, 4001):
        table = qr_table(p)
        for _ in range(6):
            a, b = rng.randrange(p), rng.randrange(p)
            if (4 * a ** 3 + 27 * b * b) % p == 0:
                continue
            assert weierstrass_order(a, b, p, table) == curve_order(a, b, p)


def test_suyama_orders_are_in_the_hasse_interval_and_divisible_by_12():
    rng = random.Random(5)
    resolved = 0
    for p in (8009, 32771):
        table = qr_table(p)
        for _ in range(10):
            order = suyama_order(rng.randrange(6, p), p, table)
            if order is None:
                continue
            resolved += 1
            assert abs(order - p - 1) <= 2 * isqrt(p) + 2
            assert order % 12 == 0, (p, order)
    assert resolved >= 15


def test_powersmooth_predicate_matches_stage_one():
    assert is_b1b2_smooth(2 * 3 * 5 * 7, 200, 200)
    assert not is_b1b2_smooth(2 * 211, 200, 200)
    assert not is_b1b2_smooth(2 ** 9, 200, 200)     # 512 > 200: stage 1 stops at 2^7


def test_multiplicative_order_can_be_smoother_than_p_minus_one():
    assert multiplicative_order(2, 7) == 3
    assert multiplicative_order(3, 7) == 6
    for p in (1009, 65537):
        assert pow(2, multiplicative_order(2, p), p) == 1
    # 2 has order 22 modulo 683, far smoother than 683 - 1 = 2 * 11 * 31
    assert multiplicative_order(2, 683) == 22


def test_shifted_primes_beat_the_baseline_which_dickman_predicts():
    rng = random.Random(7)
    from aksfactor.arith import is_prime

    primes, randoms = [], []
    while len(primes) < 300:
        v = rng.randrange(1 << 14, 1 << 15) | 1
        if is_prime(v):
            primes.append(v)
        randoms.append(rng.randrange(1 << 14, 1 << 15))
    shifted = smooth_rate([p - 1 for p in primes], 100)
    baseline = smooth_rate(randoms[:300], 100)
    assert shifted > baseline
    assert 0 < dickman_rho(15 * 0.693 / 4.6) < 1
