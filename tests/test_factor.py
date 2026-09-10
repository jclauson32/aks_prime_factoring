"""The Pascal factorizer must agree with a conventional one."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from aksfactor.arith import factorize, is_prime
from aksfactor.factor import certificate, factor, pascal_split, pascal_spf


def test_factor_matches_reference():
    for n in list(range(2, 400)) + [1234567, 999983 * 3, 2**16 + 1, 5**7 * 11]:
        assert factor(n) == factorize(n), n


def test_scan_mode_matches_certified():
    for n in range(4, 260):
        if is_prime(n):
            continue
        assert pascal_split(n, mode="scan") == pascal_split(n, mode="certified"), n


def test_certificate_is_exact_cofactor():
    for n, p in [(1234567, 127), (35, 5), (1155, 3), (1009 * 1013, 1009)]:
        cert = certificate(n, p)
        assert cert["valid"], cert
        assert cert["residue"] == n // p


def test_spf_of_prime_is_itself():
    for n in (2, 3, 97, 1009, 1000003):
        assert pascal_spf(n) == n


def test_bound_returns_none():
    n = 1000003 * 1000033
    assert pascal_spf(n, bound=1000) is None
    assert pascal_split(n, bound=1000) is None
