"""The O~(n^1/4) multipoint search must agree with trial division exactly."""

import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import aksfactor.fast as fast
from aksfactor.arith import is_prime, spf_trial
from aksfactor.fast import multipoint_eval, poly_mul, product_of


def _horner(f, x, n):
    v = 0
    for c in reversed(f):
        v = (v * x + c) % n
    return v


def test_poly_mul_matches_schoolbook():
    rng = random.Random(1)
    for n in (10**9 + 7, 2**61 - 1, 1001 * 9973):
        for _ in range(60):
            a = [rng.randrange(n) for _ in range(rng.randrange(1, 45))]
            b = [rng.randrange(n) for _ in range(rng.randrange(1, 45))]
            ref = [0] * (len(a) + len(b) - 1)
            for i, x in enumerate(a):
                for j, y in enumerate(b):
                    ref[i + j] = (ref[i + j] + x * y) % n
            while ref and ref[-1] == 0:
                ref.pop()
            assert poly_mul(a, b, n) == ref


def test_multipoint_eval_matches_horner():
    rng = random.Random(2)
    for n in (10**9 + 7, 2**64, 6, 1001 * 9973):
        for _ in range(50):
            pts = [rng.randrange(n) for _ in range(rng.randrange(1, 33))]
            f = [rng.randrange(n) for _ in range(rng.randrange(1, 60))]
            got = multipoint_eval(f, pts, n)
            assert got == [_horner(f, x, n) for x in pts]


def test_product_of():
    n = 10**9 + 7
    f = product_of([[j, 1] for j in range(1, 9)], n)
    for x in (0, 1, 5, 12345):
        want = 1
        for j in range(1, 9):
            want = want * (x + j) % n
        assert _horner(f, x, n) == want


def test_fast_spf_exhaustive():
    saved = fast._CROSSOVER
    fast._CROSSOVER = 0  # force the polynomial path on small inputs too
    try:
        for n in range(4, 1500):
            assert fast.fast_spf(n) == spf_trial(n), n
    finally:
        fast._CROSSOVER = saved


def test_fast_spf_structured():
    saved = fast._CROSSOVER
    fast._CROSSOVER = 0
    try:
        for n in (1009 * 1013, 65537**2, 999983 * 999979, 104729 * 1299709, 2**31 - 1):
            assert fast.fast_spf(n) == spf_trial(n), n
        for p in (1000003, 999983):
            assert fast.fast_spf(p) == p
    finally:
        fast._CROSSOVER = saved


def test_fast_spf_bound():
    n = 999983 * 999979
    assert fast.fast_spf(n, bound=1000) is None
    assert fast.fast_spf(n, bound=10**6) == 999979


def test_fast_split():
    assert fast.fast_split(1009 * 1013) == (1009, 1013)
    assert fast.fast_split(1000003) is None


def test_covers_more_positions_than_gcds():
    """The whole point: many positions tested per gcd."""
    stats = {}
    n = 200210159 * 400420351
    fast.fast_spf(n, bound=300000, stats=stats)
    assert stats["c"] > 100
    assert stats["positions_covered"] >= 300000


def test_factorial_mod_matches_math_factorial():
    from math import factorial

    from aksfactor.fast import factorial_mod

    for n in range(2, 300):
        for m in range(0, 45):
            assert factorial_mod(m, n) == factorial(m) % n, (m, n)


def test_factorial_mod_large():
    from math import factorial

    from aksfactor.fast import factorial_mod

    n = 1009 * 1013
    for m in (500, 1000, 1010, 2000):
        assert factorial_mod(m, n) == factorial(m) % n, m


def test_threshold_spf_matches_trial_division():
    from aksfactor.arith import spf_trial
    from aksfactor.fast import threshold_spf

    for n in list(range(4, 500)) + [1009 * 1013, 10007 * 10009, 65537**2,
                                    999983 * 999979, 2**31 - 1, 1000003]:
        assert threshold_spf(n) == spf_trial(n), n


def test_factorial_threshold_is_monotone():
    """T20: the predicate binary search relies on."""
    from math import gcd

    from aksfactor.arith import spf_trial
    from aksfactor.fast import factorial_mod

    for n in (91, 143, 221, 1147, 2021, 1009 * 1013):
        s = spf_trial(n)
        for i in range(1, min(n, 80)):
            assert (gcd(factorial_mod(i, n), n) > 1) == (i >= s), (n, i)
