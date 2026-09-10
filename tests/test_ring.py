"""The folded AKS object: identity, correctness, and factor leakage."""

import sys
from math import comb
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from aksfactor.arith import factorize
from aksfactor.ring import aks_pow, fold_attack, fold_coefficients, fold_identity


def test_aks_pow_matches_direct_fold():
    for n in (12, 21, 35, 45, 77):
        for r in (2, 3, 5, 7, 11):
            for a in (1, 2, 3):
                got = aks_pow(n, r, a)
                want = [0] * r
                for k in range(n + 1):
                    want[k % r] = (want[k % r] + comb(n, k) * pow(a, n - k)) % n
                assert got == want, (n, r, a)


def test_prime_satisfies_aks_identity():
    for n in (5, 7, 11, 13, 31, 101):
        for r in (3, 5, 8):
            for a in (1, 2):
                got = aks_pow(n, r, a)
                want = [0] * r
                want[n % r] = (want[n % r] + 1) % n
                want[0] = (want[0] + a) % n
                assert got == want, (n, r, a)


def test_fold_identity():
    for n in (35, 36, 77, 100, 120, 121):
        for p in factorize(n):
            for r in (3, 5, 7, 9):
                for a in (1, 2, 3):
                    ok, detail = fold_identity(n, p, r, a)
                    assert ok, detail


def test_fold_attack_finds_factor():
    for n in (1009 * 1013, 29651, 12317):
        res = fold_attack(n, rmax=120)
        assert res["factor"] is not None
        assert n % res["factor"] == 0
        assert 1 < res["factor"] < n


def test_fold_attack_reports_cost():
    res = fold_attack(1009 * 1013, rmax=120)
    assert res["pairs"] > 0


def test_cyclic_kronecker_path():
    # r > 24 forces the packed-multiplication branch of _cyclic_mul
    for n in (501, 1024, 2047):
        got = aks_pow(n, 64, 3)
        want = [0] * 64
        for k in range(n + 1):
            want[k % 64] = (want[k % 64] + comb(n, k) * pow(3, n - k)) % n
        assert got == want, n


def test_fold_coefficients_alias():
    assert fold_coefficients(35, 5, 1) == aks_pow(35, 5, 1)
