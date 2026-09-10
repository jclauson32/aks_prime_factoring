"""Norm-one cyclotomic factoring: the method that refutes round 2's domination claim."""

import sys
from math import gcd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from aksfactor.cyclo import (
    _norm,
    lucas_v,
    lucas_v_binomial,
    norm_one_search,
    pollard_pminus1,
    williams_pplus1,
)


def test_lucas_matches_recurrence():
    for n in (10**9 + 7, 1001 * 9973, 97):
        for a in (3, 5, 7, 11):
            v0, v1 = 2 % n, a % n
            for k in range(0, 60):
                assert lucas_v(k, a, n) == v0, (n, a, k)
                v0, v1 = v1, (a * v1 - v0) % n


def test_lucas_matches_binomial_sum():
    for n in (10**9 + 7, 1001 * 9973, 97):
        for a in (3, 5, 7, 11):
            for k in range(0, 45):
                assert lucas_v(k, a, n) == lucas_v_binomial(k, a, n), (n, a, k)


def test_norm_degree_two_closed_form():
    n = 10**9 + 7
    for a in (3, 5, 7):
        f = [1 % n, (-a) % n, 1 % n]  # x^2 - a x + 1
        for c0 in (0, 1, 2, 5, 123):
            for c1 in (0, 1, 3, 77):
                # roots satisfy z1+z2 = a, z1*z2 = 1
                want = (c0 * c0 + c0 * c1 * a + c1 * c1) % n
                assert _norm([c0, c1], f, n) == want, (a, c0, c1)


def test_williams_finds_smooth_p_plus_one():
    # p-1 is rough, p+1 is 60-smooth; q is rough on both sides
    q = 1419449461
    for p in (1286450213, 1219421579, 1820048873, 2938272193, 1277805479):
        n = p * q
        found = norm_one_search(n, bound=500, degree=2, bases=range(3, 60))
        assert found is not None, p
        assert found[0] == p, (p, found)


def test_pollard_fails_where_williams_succeeds():
    """The point of the correction: these are not reachable from the p-1 side."""
    q = 1419449461
    for p in (1286450213, 1820048873, 1277805479):
        n = p * q
        assert pollard_pminus1(n, bound=500) is None, p


def test_norm_one_returns_genuine_divisors():
    for n in (1009 * 1013, 10007 * 10009, 65537 * 65539):
        found = norm_one_search(n, bound=300, degree=2, bases=range(3, 30))
        if found:
            assert 1 < found[0] < n and n % found[0] == 0


def test_degree_three_splits():
    ok = 0
    cases = [(10007, 10009), (12007, 13009), (15013, 16007)]
    for p, q in cases:
        n = p * q
        found = norm_one_search(n, bound=300, degree=3, bases=range(3, 12))
        if found and 1 < found[0] < n and n % found[0] == 0:
            ok += 1
    assert ok >= 1


def test_exponent_bound_matters():
    """lcm(1..B) caps prime powers at B; p+1 with a high power of 2 needs more."""
    from aksfactor.cyclo import _prime_power_ladder

    small = list(_prime_power_ladder(500, exponent_bound=500))
    large = list(_prime_power_ladder(500))
    assert max(small) <= 500
    assert max(large) > 500
