"""The three row engines must agree, and random access must be exact."""

import sys
from math import comb, gcd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from aksfactor.pascal import (
    first_nonzero,
    residue_shape,
    row_entry,
    row_exact,
    row_series,
    row_support,
)


def test_engines_agree():
    for n in range(2, 260):
        kmax = min(n, 60)
        exact = row_exact(n, kmax)
        assert row_series(n, kmax) == exact, n
        assert [row_entry(n, k) for k in range(kmax + 1)] == exact, n


def test_series_full_row():
    for n in (30, 64, 97, 105, 121):
        assert row_series(n, n) == [comb(n, k) % n for k in range(n + 1)]


def test_random_access_large_n():
    # A 67-digit modulus; the residue at k = p must be the cofactor exactly.
    p, q = 1000003, 2**200 - 75 - 12
    while True:
        from aksfactor.arith import is_prime

        if is_prime(q):
            break
        q -= 2
    n = p * q
    assert row_entry(n, p) == q
    assert row_entry(n, p - 1) == 0  # p-1 is coprime to n here
    assert row_entry(n, 2 * p) % (n // gcd(n, 2 * p)) == 0


def test_support_only_on_shared_factors():
    for n in range(4, 200):
        for k in row_support(n):
            assert gcd(n, k) > 1, (n, k)


def test_first_nonzero_prime_is_none():
    for n in (2, 3, 5, 7, 97, 101, 1009):
        assert first_nonzero(n) is None


def test_residue_shape():
    shape = residue_shape(35, 10)
    assert shape["residue"] == 21
    assert shape["x"] == 7 and shape["y"] == 3
    assert shape["x_divides_n"] and shape["x_divides_residue"]
    assert shape["gcd_residue_n"] == 7
