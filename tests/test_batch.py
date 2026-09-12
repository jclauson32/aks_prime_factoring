import random

from aksfactor.arith import sieve
from aksfactor.batch import (
    batch_smooth,
    product_tree,
    remainder_tree,
    smooth_parts,
    trial_smooth,
)


def test_product_tree_multiplies_and_remainder_tree_reduces():
    values = [3, 5, 7, 11, 13]
    levels = product_tree(values)
    root = levels[-1][0]
    assert root == 3 * 5 * 7 * 11 * 13
    big = 1234567891011
    assert remainder_tree(big, levels) == [big % v for v in values]


def test_odd_sized_batches_still_reduce_correctly():
    rng = random.Random(1)
    for count in (1, 2, 3, 7, 9, 33):
        values = [rng.randrange(1 << 20, 1 << 21) for _ in range(count)]
        levels = product_tree(values)
        big = rng.randrange(1 << 200)
        assert remainder_tree(big, levels) == [big % v for v in values]


def test_batch_agrees_with_trial_division():
    rng = random.Random(2)
    primes = sieve(500)
    values = [rng.randrange(1 << 40, 1 << 41) for _ in range(200)]
    values += [2 ** 5 * 3 ** 7 * 5 ** 3, 7 ** 11, 499 ** 4, 503 * 499]
    assert batch_smooth(values, primes) == trial_smooth(values, primes)
    assert any(batch_smooth(values, primes))          # the planted ones show up


def test_smooth_part_is_the_largest_smooth_divisor():
    primes = sieve(100)
    values = [2 ** 8 * 3 * 101, 97 * 103, 2 * 3 * 5 * 7, 101 ** 3]
    parts = smooth_parts(values, primes)
    assert parts == [2 ** 8 * 3, 97, 2 * 3 * 5 * 7, 1]
    for v, part in zip(values, parts):
        rough = v // part
        assert v % part == 0
        assert all(rough % ell for ell in primes)


def test_empty_and_singleton_batches():
    primes = sieve(50)
    assert batch_smooth([], primes) == []
    assert smooth_parts([], primes) == []
    assert batch_smooth([49], primes) == [True]
    assert batch_smooth([53 * 59], primes) == [False]


def test_batch_gcd_finds_a_shared_prime():
    from math import gcd

    from aksfactor.batch import batch_gcd

    p, q, r, s = 1000003, 1000033, 1000037, 1000039
    moduli = [p * q, p * r, s * 1000081, 1000081 * 1000099]
    got = batch_gcd(moduli)
    assert got[0] == p and got[1] == p          # the shared prime, from one pass
    assert got[2] == 1000081 and got[3] == 1000081
    # and it agrees with the pairwise loop it replaces
    for i, v in enumerate(moduli):
        want = 1
        for j, w in enumerate(moduli):
            if i != j:
                want = want * w
        assert got[i] == gcd(v, want)


def test_batch_gcd_reports_nothing_when_nothing_is_shared():
    from aksfactor.batch import batch_gcd

    moduli = [1000003 * 1000033, 1000037 * 1000039, 1000081 * 1000099]
    assert batch_gcd(moduli) == [1, 1, 1]
    assert batch_gcd([15]) == [1]
    assert batch_gcd([]) == []
