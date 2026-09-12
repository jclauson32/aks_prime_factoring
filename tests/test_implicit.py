import random

from aksfactor.arith import is_prime
from aksfactor.implicit import (
    implicit_factor,
    shared_low_primes,
    shared_lsb_lattice,
    threshold_bits,
)


def test_threshold_falls_towards_alpha_and_never_below():
    assert threshold_bits(2, 40) == 80
    assert threshold_bits(3, 40) == 60
    assert threshold_bits(10, 40) == 45
    for k in range(2, 40):
        assert threshold_bits(k, 40) > 40


def test_lattice_contains_the_cofactors():
    rng = random.Random(11)
    t = 90
    pairs = shared_low_primes(3, 180, t, 40, rng)
    moduli = [p * q for p, q in pairs]
    basis = shared_lsb_lattice(moduli, t)
    m = 1 << t
    for row in basis:                       # every basis row satisfies the congruence
        for i in range(1, 3):
            assert (row[0] * moduli[i] - row[i] * moduli[0]) % m == 0
    q0 = pairs[0][1]
    target = [q for _, q in pairs]
    for i in range(1, 3):                   # and so does the vector we are hunting
        assert (q0 * moduli[i] - target[i] * moduli[0]) % m == 0


def test_factors_above_the_threshold_and_fails_below():
    rng = random.Random(12)
    alpha = 40
    for k in (2, 3):
        pairs = shared_low_primes(k, 200, 100, alpha, rng)
        moduli = [p * q for p, q in pairs]
        want = sorted(p for p, _ in pairs)
        thr = threshold_bits(k, alpha)
        assert sorted(implicit_factor(moduli, thr + 10) or []) == want
        assert implicit_factor(moduli, thr - 12) is None


def test_unrelated_moduli_are_not_factored():
    rng = random.Random(13)
    moduli = []
    for _ in range(3):
        while True:
            p = rng.randrange(1 << 199, 1 << 200) | 1
            if is_prime(p):
                break
        while True:
            q = rng.randrange(1 << 39, 1 << 40) | 1
            if is_prime(q):
                break
        moduli.append(p * q)
    assert implicit_factor(moduli, 100) is None


def test_shared_primes_really_share_and_differ():
    rng = random.Random(14)
    pairs = shared_low_primes(4, 120, 60, 30, rng)
    primes = [p for p, _ in pairs]
    assert len(set(primes)) == 4                    # equal primes would be a gcd
    assert len({p % (1 << 60) for p in primes}) == 1
    for p, q in pairs:
        assert is_prime(p) and is_prime(q)
        assert p.bit_length() == 120 and q.bit_length() == 30


def test_generator_refuses_an_impossible_residue_class():
    """Asking for more primes than the class holds must fail, not spin."""
    rng = random.Random(15)
    try:
        shared_low_primes(8, 60, 54, 60, rng)     # 32 candidates, ~0.8 prime
    except ValueError as exc:
        assert "too few" in str(exc)
    else:
        raise AssertionError("expected a refusal")
