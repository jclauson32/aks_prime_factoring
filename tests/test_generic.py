"""Generic ring programs achieve the zero-divisor density, and nothing more."""

import random
import sys
from math import gcd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from aksfactor.arith import is_prime
from aksfactor.generic import (
    alignment_gain,
    generic_rate,
    random_program,
    run_program,
    zero_divisor_density,
)


def test_run_program_stays_in_the_ring():
    rng = random.Random(1)
    n = 1009 * 1013
    for _ in range(50):
        ops = random_program(10, rng)
        v = run_program(ops, rng.randrange(2, n), n)
        assert 0 <= v < n


def test_program_is_deterministic():
    rng = random.Random(2)
    n = 10007 * 10009
    ops = random_program(12, rng)
    a = 12345
    assert run_program(ops, a, n) == run_program(ops, a, n)


def test_zero_divisor_density_matches_uniform_sampling():
    rng = random.Random(3)
    for p, q in [(113, 281), (353, 607)]:
        n = p * q
        hits = 0
        trials = 60000
        for _ in range(trials):
            if 1 < gcd(rng.randrange(1, n), n) < n:
                hits += 1
        predicted = zero_divisor_density(n, p, q)
        assert abs(hits / trials / predicted - 1) < 0.1, (p, q, hits / trials)


def test_generic_programs_hit_the_baseline_not_better():
    """Deep random programs match 1 - phi(N)/N; they never beat it."""
    for p, q in [(179, 419), (733, 761)]:
        n = p * q
        base = zero_divisor_density(n, p, q)
        rate = generic_rate(n, 24, 60000, random.Random(17))
        assert rate < 2 * base, (p, q, rate, base)
        assert rate > 0.4 * base, (p, q, rate, base)


def test_alignment_beats_the_baseline_by_orders_of_magnitude():
    """The same operations with a chosen exponent win by 100x or more."""
    p, q = 4643, 4733
    assert is_prime(p) and is_prime(q)
    n = p * q
    base = zero_divisor_density(n, p, q)
    got = alignment_gain(n, 200, 400, random.Random(1))
    assert got["rate"] / base > 100, got


def test_alignment_needs_the_order_to_be_smooth():
    """No smoothness, no gain -- the alignment is with |G|, not with arithmetic."""
    # p - 1 and q - 1 both have a large prime factor
    p, q = 1000003, 1000033
    assert is_prime(p) and is_prime(q)
    n = p * q
    got = alignment_gain(n, 30, 200, random.Random(2))
    assert got["rate"] == 0.0, got
