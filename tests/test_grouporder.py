"""Group-order counting: the measurement behind Proposition 14."""

import random
import sys
from math import isqrt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from aksfactor.arith import is_prime
from aksfactor.grouporder import (
    degree_pattern,
    elliptic_order,
    partitions,
    poly_gcd,
    ring_unit_order,
)


def test_partitions():
    assert [partitions(d) for d in range(1, 8)] == [1, 2, 3, 5, 7, 11, 15]


def test_poly_gcd():
    p = 101
    # (x-1)(x-2) and (x-2)(x-3) share (x-2)
    a = [2, -3, 1]
    b = [6, -5, 1]
    g = poly_gcd([c % p for c in a], [c % p for c in b], p)
    assert g == [(-2) % p, 1]


def test_degree_pattern_known_cases():
    p = 101  # 101 = 1 mod 4, so x^2+1 splits
    assert degree_pattern([1, 0, 1], p) == [1, 1]
    p = 103  # 103 = 3 mod 4, so x^2+1 is irreducible
    assert degree_pattern([1, 0, 1], p) == [2]


def test_degree_pattern_rejects_non_squarefree():
    p = 97
    assert degree_pattern([1, 2, 1], p) is None  # (x+1)^2


def test_ring_unit_order_matches_formula():
    p = 101
    assert ring_unit_order([1, 0, 1], p) == (p - 1) ** 2
    p = 103
    assert ring_unit_order([1, 0, 1], p) == p * p - 1


def test_ring_orders_capped_by_partitions():
    """The heart of Proposition 14: the ring cannot exceed partitions(d) orders."""
    rng = random.Random(5)
    for p in (211, 503, 1009):
        for d in (2, 3, 4):
            seen = set()
            for _ in range(300):
                f = [rng.randrange(p) for _ in range(d)] + [1]
                order = ring_unit_order(f, p)
                if order is not None:
                    seen.add(order)
            assert len(seen) <= partitions(d), (p, d, len(seen))


def test_elliptic_order_respects_hasse():
    for p in (211, 503, 1009):
        assert is_prime(p)
        for a in range(0, 6):
            for b in range(0, 6):
                order = elliptic_order(a, b, p)
                if order is None:
                    continue
                assert abs(order - (p + 1)) <= 2 * isqrt(p) + 1, (p, a, b, order)


def test_elliptic_reaches_more_orders_than_ring():
    """Varying beats rigid, measurably."""
    rng = random.Random(6)
    p = 1009
    ring = set()
    for _ in range(300):
        f = [rng.randrange(p) for _ in range(4)] + [1]
        o = ring_unit_order(f, p)
        if o is not None:
            ring.add(o)
    ec = set()
    for _ in range(300):
        o = elliptic_order(rng.randrange(p), rng.randrange(p), p)
        if o is not None:
            ec.add(o)
    assert len(ring) <= partitions(4)
    assert len(ec) > 10 * len(ring)
