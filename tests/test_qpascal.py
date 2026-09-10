"""The q-deformed Pascal triangle: criterion, and what it does and doesn't buy."""

import sys
from math import gcd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from aksfactor.arith import spf_trial
from aksfactor.pascal import row_entry
from aksfactor.qpascal import (
    q_clause_one_position,
    q_first_hit,
    q_lucas_divides,
    q_pascal_row,
    q_period,
    shift_split,
)


def test_q_period():
    assert q_period(2, 7) == 3       # 2^3 = 1 mod 7
    assert q_period(3, 7) == 6       # 3 is a primitive root mod 7
    assert q_period(1, 97) == 1
    assert q_period(7, 7) is None


def test_q_pascal_row_matches_definition():
    """[n,k]_q evaluated from the polynomial definition."""
    from math import prod

    for n in (5, 6, 7, 9):
        for q in (2, 3, 5):
            row = q_pascal_row(n, n, q, 10**9 + 7)
            for k in range(n + 1):
                num = prod(q**(n - k + i) - 1 for i in range(1, k + 1))
                den = prod(q**i - 1 for i in range(1, k + 1))
                assert num % den == 0
                assert row[k] == (num // den) % (10**9 + 7), (n, k, q)


def test_q_lucas_criterion_exhaustive():
    """p | [n,k]_q  <=>  low-digit carry in base d, or p | C(n//d, k//d)."""
    bad = tot = 0
    for p in (7, 11, 13, 17, 19, 23):
        for base in (2, 3, 5, 6, 7, 10):
            if base % p == 0:
                continue
            d = q_period(base, p)
            for n in (p * 11, p * 13, p * 17):
                if n > 400:
                    continue
                row = q_pascal_row(n, min(40, n), base, p)
                for k in range(1, len(row)):
                    tot += 1
                    if (row[k] % p == 0) != q_lucas_divides(n, k, d, p):
                        bad += 1
    assert tot > 3000
    assert bad == 0, f"{bad}/{tot} mismatches"


def test_clause_one_position_predicts_divisibility():
    for p, r in [(7, 11), (11, 13), (13, 17), (19, 23)]:
        n = p * r
        for q in range(2, 30):
            if gcd(q, n) != 1:
                continue
            for prime in (p, r):
                d = q_period(q, prime)
                k = q_clause_one_position(n, d)
                if k is None:
                    continue
                row = q_pascal_row(n, min(k, 60), q, prime)
                if len(row) > k:
                    assert row[k] % prime == 0, (n, q, prime, k)


def test_shift_split_matches_order_condition():
    """shift_split finds the least j where exactly one order divides n-j."""
    for p, r in [(7, 11), (11, 13), (13, 17), (17, 19), (11, 23), (23, 29)]:
        n = p * r
        for q in range(2, 40):
            if gcd(q, n) != 1:
                continue
            dp, dr = q_period(q, p), q_period(q, r)
            want = None
            for j in range(0, 200):
                a = (n - j) % dp == 0
                b = (n - j) % dr == 0
                if a != b:
                    want = (j, p if a else r)
                    break
            assert shift_split(n, q, 200) == want, (n, q)


def test_q_row_fires_below_spf():
    """The deformation is real: the q-row can be non-trivial where the
    classical row is provably zero."""
    found = 0
    for p, r in [(11, 13), (13, 17), (17, 19), (19, 23), (23, 29)]:
        n = p * r
        s = spf_trial(n)
        for q in range(2, 30):
            if gcd(q, n) != 1:
                continue
            hit = q_first_hit(n, s - 1, q)
            if hit:
                # classical row is identically zero at every k < spf(n)
                assert all(row_entry(n, k) == 0 for k in range(1, s))
                assert hit[0] < s
                found += 1
    assert found > 0


def test_shift_split_returns_real_factors():
    for n in (1009 * 1013, 10007 * 10009, 77, 143):
        for q in (2, 3, 5, 7):
            got = shift_split(n, q, 200)
            if got:
                assert 1 < got[1] < n and n % got[1] == 0
