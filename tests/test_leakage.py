import random
from math import gcd

from aksfactor.arith import is_prime
from aksfactor.leakage import (
    bayes_sum_accuracy,
    feature_accuracy,
    r4,
    tau_prime_mod_691,
    tau_sum_candidates,
    tau_table,
    trace_divisor_term,
)

KNOWN_TAU = [0, 1, -24, 252, -1472, 4830, -6048, -16744, 84480, -113643, -115920]


def _sigma(k, n):
    return sum(d ** k for d in range(1, n + 1) if n % d == 0)


def test_tau_known_values():
    assert tau_table(len(KNOWN_TAU)) == KNOWN_TAU


def test_tau_ramanujan_congruence_and_multiplicativity():
    tau = tau_table(300)
    assert all((tau[n] - _sigma(11, n)) % 691 == 0 for n in range(1, 300))
    for m in range(2, 17):
        for n in range(2, 17):
            if m * n < 300 and gcd(m, n) == 1:
                assert tau[m * n] == tau[m] * tau[n]


def test_tau_parity_is_odd_squares():
    tau = tau_table(400)
    for n in range(1, 400):
        r = int(n ** 0.5 + 0.5)
        odd_square = r * r == n and r % 2 == 1
        assert (tau[n] % 2 == 1) == odd_square


def test_tau_pins_the_sum_mod_691():
    rng = random.Random(3)
    sizes = []
    for _ in range(40):
        while True:
            p = rng.randrange(10 ** 4, 10 ** 5) | 1
            if is_prime(p) and p % 691:
                break
        while True:
            q = rng.randrange(p + 2, 2 * p) | 1
            if is_prime(q) and q % 691:
                break
        t = tau_prime_mod_691(p) * tau_prime_mod_691(q) % 691
        got = tau_sum_candidates(p * q, t)
        assert (p + q) % 691 in got
        sizes.append(len(got))
    assert sum(sizes) / len(sizes) < 1.2


def test_trace_formula_names_the_smaller_factor():
    for p, q in [(3, 5), (5, 7), (11, 13), (29, 31)]:
        assert trace_divisor_term(p * q) == 2 + 2 * p ** 11


def test_four_squares_give_p_plus_q():
    for p, q in [(3, 5), (3, 7), (5, 7), (7, 11), (3, 13)]:
        n = p * q
        assert r4(n) == 8 * (n + 1 + p + q)


def test_bayes_baseline_matches_enumeration():
    for ell in (5, 7, 11, 13, 17):
        best = 0
        for r in range(1, ell):
            counts = {}
            for a in range(1, ell):
                s = (a + r * pow(a, -1, ell)) % ell
                counts[s] = counts.get(s, 0) + 1
            best += max(counts.values()) / (ell - 1)
        assert abs(best / (ell - 1) - bayes_sum_accuracy(ell)) < 1e-12


def test_feature_detector_positive_control():
    """Close primes: floor(sqrt N) nearly determines p+q, and the detector sees it."""
    from math import isqrt

    rng = random.Random(8)
    data = []
    while len(data) < 1200:
        p = rng.randrange(1 << 15, 1 << 16) | 1
        if not is_prime(p):
            continue
        q = p + 2
        while not is_prime(q):
            q += 2
        data.append((p * q, p, q))
    acc = feature_accuracy(data, 7, lambda n, ell: isqrt(n) % ell)
    assert acc > 0.9
    noise = feature_accuracy(data, 7, lambda n, ell: bin(n).count("1") % ell)
    assert noise < 0.9
