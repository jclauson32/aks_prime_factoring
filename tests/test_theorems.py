"""Exhaustive verification of T1-T5 and the original conjecture."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from aksfactor.theorems import (
    check_all,
    check_gist,
    check_t1_shape,
    check_t2_coprime_vanishing,
    check_t3_exact_value,
    check_t4_first_nonzero,
    check_t5_closed_form,
)

RANGE = range(4, 320)


def test_t1_shape():
    for n in RANGE:
        ok, detail = check_t1_shape(n)
        assert ok, detail


def test_t2_coprime_vanishing():
    for n in RANGE:
        ok, detail = check_t2_coprime_vanishing(n)
        assert ok, detail


def test_t3_exact_value():
    for n in RANGE:
        ok, detail = check_t3_exact_value(n)
        assert ok, detail


def test_t4_first_nonzero():
    for n in RANGE:
        ok, detail = check_t4_first_nonzero(n)
        assert ok, detail


def test_t5_closed_form():
    for n in range(4, 160):
        ok, detail = check_t5_closed_form(n)
        assert ok, detail


def test_gist_holds():
    for n in RANGE:
        ok, detail = check_gist(n)
        assert ok, detail


def test_check_all_bundles():
    for n in (12, 35, 97, 561, 1024, 2047):
        for name, (ok, detail) in check_all(n, kmax=200).items():
            assert ok, (n, name, detail)


def test_t7_aliasing():
    from math import gcd

    from aksfactor.theorems import check_t7_aliasing

    for n in (35, 77, 143, 221, 1001, 105):
        for r in range(2, 20):
            if gcd(r, n) != 1:
                continue
            ok, detail = check_t7_aliasing(n, r)
            assert ok, detail


def test_t9_second_digit():
    from aksfactor.theorems import check_t9_second_digit

    for n in (15, 35, 77, 143, 221, 1001, 2257, 105, 1155):
        ok, detail = check_t9_second_digit(n)
        assert ok, detail
        # the split rate is a constant, not a 1/p lottery
        if detail.get("rate") is not None and detail["trials"] > 50:
            assert detail["rate"] > 0.05, detail


def test_t20_factorial_threshold():
    from aksfactor.theorems import check_t20_factorial_threshold

    for n in (12, 60, 91, 143, 210, 875, 1155, 2310, 9009, 1024, 1009 * 1013):
        ok, detail = check_t20_factorial_threshold(n)
        assert ok, detail


def test_t21_stride_embedding():
    from aksfactor.arith import sieve
    from aksfactor.theorems import check_t21_stride_embedding

    primes = [x for x in sieve(120) if x > 2]
    for i, p in enumerate(primes):
        for q in primes[i + 1:i + 6]:
            ok, detail = check_t21_stride_embedding(p, q)
            assert ok, detail


def test_t21_support_count_follows_fine():
    """All multiples of the larger prime are non-zero; the smaller prime's
    multiples follow Fine's theorem -- few for close primes, many otherwise."""
    from aksfactor.pascal import row_entry

    for p, q, expect in [(1009, 1013, 8), (1009, 1511, 1004), (3, 11, 4)]:
        n = p * q
        at_q = sum(1 for m in range(1, p) if row_entry(n, m * q))
        at_p = sum(1 for m in range(1, q) if row_entry(n, m * p))
        assert at_q == p - 1
        assert at_p == expect
