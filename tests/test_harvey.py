"""Harvey's exponent-one-fifth deterministic factoring algorithm (arXiv:2010.05450)."""

import random
import sys
from math import ceil, isqrt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from aksfactor.arith import is_prime, sieve
from aksfactor.harvey import (
    collision_search,
    harvey_factor,
    harvey_search,
    lehman_recover,
)


def test_lehman_recover_on_true_candidates():
    """Lemma 3.1: u = a*q + b*p must recover (p, q)."""
    for p, q in [(11, 13), (101, 103), (1009, 1013), (10007, 10009), (7, 101)]:
        n = p * q
        for a in range(1, 6):
            for b in range(1, 6):
                assert lehman_recover(n, a, b, a * q + b * p) == (p, q), (p, q, a, b)


def test_lehman_recover_only_returns_real_factorizations():
    """Any non-None answer must be a genuine factorization."""
    for p, q in [(101, 103), (1009, 1013)]:
        n = p * q
        for a in range(1, 4):
            for b in range(1, 4):
                for d in range(-30, 30):
                    got = lehman_recover(n, a, b, a * q + b * p + d)
                    if got is not None:
                        assert got[0] * got[1] == n and 1 < got[0] <= got[1]


def test_lehman_rejects_non_squares():
    assert lehman_recover(1009 * 1013, 1, 1, 3) is None
    assert lehman_recover(143, 1, 1, 0) is None


def test_algorithm_42_adjacent_primes():
    primes = [x for x in sieve(4000) if x > 50]
    ok = tot = 0
    for i in range(0, len(primes) - 1, 5):
        p, q = primes[i], primes[i + 1]
        n = p * q
        r = max(1, ceil(q / p) + 1)
        m = max(2, isqrt(isqrt(n)) + 2)
        tot += 1
        ok += harvey_search(n, r, m, 2) == (p, q)
    assert tot > 50
    assert ok == tot, f"{ok}/{tot}"


def test_algorithm_42_unbalanced():
    rng = random.Random(4)
    primes = [x for x in sieve(6000) if x > 50]
    ok = tot = 0
    for _ in range(40):
        p = rng.choice(primes[:150])
        q = rng.choice(primes)
        if p >= q:
            continue
        n = p * q
        r = max(1, ceil(q / p) + 1)
        m = max(2, isqrt(isqrt(n)) + 2)
        tot += 1
        ok += harvey_search(n, r, m, 2) == (p, q)
    assert ok == tot, f"{ok}/{tot}"


def test_algorithm_42_reports_primes():
    for x in [x for x in sieve(20000) if x > 10000][:15]:
        assert harvey_search(x, 3, 30, 2) is None


def test_algorithm_43_end_to_end():
    rng = random.Random(7)

    def randprime(lo, hi):
        while True:
            x = rng.randrange(lo, hi) | 1
            if is_prime(x):
                return x

    cases = [(101, 103), (1009, 1013), (10007, 10009), (7, 101), (11, 10007)]
    for lo, hi in [(100, 400), (1000, 4000), (10000, 40000)]:
        for _ in range(4):
            p = randprime(lo, hi)
            q = randprime(p + 2, 3 * p)
            cases.append((p, q))
    for p, q in cases:
        n = p * q
        assert harvey_factor(n) == (min(p, q), max(p, q)), n


def test_algorithm_43_primes_return_none():
    for x in (10007, 65537, 99991, 1000003, 2**31 - 1):
        assert harvey_factor(x) is None


def test_collision_search_finds_single_prime_collisions():
    n = 101 * 103
    powers = [pow(2, i, n) for i in range(20)]
    # a value congruent to 2**3 mod 101 but not mod 103
    target = 8 % 101
    v = None
    for cand in range(n):
        if cand % 101 == target and cand % 103 != pow(2, 3, 103):
            v = cand
            break
    got = collision_search(n, powers, [v])
    assert got == (101, 103)


def test_divisor_summatory_exact():
    from aksfactor.harvey import divisor_summatory

    for r in (1, 2, 10, 50, 200, 1000, 4321):
        brute = sum(1 for a in range(1, r + 1) for b in range(1, r // a + 1))
        assert divisor_summatory(r) == brute, r


def test_candidate_ratio_is_exactly_two_at_the_optimum():
    """At r = m = N^(1/5) the pair term is exactly twice the j term, for all N."""
    from aksfactor.harvey import candidate_counts

    for bits in (60, 100, 140, 180, 220):
        n = 1 << bits
        r = m = round(2.0 ** (bits * 0.2))
        c = candidate_counts(n, r, m)
        assert abs(c["pair_candidates"] / c["j_candidates"] - 2.0) < 1e-6, bits
        assert abs(c["pair_share"] - 2 / 3) < 1e-6


def test_exponent_runs_are_short():
    """No local arithmetic progressions in Lehman's range -- nothing to BSGS."""
    from aksfactor.harvey import exponent_runs, run_length_bound

    for bits in (30, 36, 40, 44):
        n = 1 << bits
        r = min(max(4, round(n ** (1 / 3))), 12000)
        e = exponent_runs(n, r)
        assert e["longest"] <= 4, (bits, e)
        # and the bound is the reason
        assert run_length_bound(r, n) < 10, (bits, run_length_bound(r, n))


def test_run_length_bound_shape():
    from aksfactor.harvey import run_length_bound

    n = 1 << 60
    # grows like b^1.5, so it only reaches sqrt(r) far beyond r = N^(1/3)
    assert run_length_bound(1000, n) < run_length_bound(10000, n)
    assert run_length_bound(round(n ** (1 / 3)), n) < 3
