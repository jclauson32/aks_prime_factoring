"""LLL and Coppersmith: the one polynomial-time family in this project."""

import random
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from aksfactor.arith import is_prime
from aksfactor.lattice import (
    _integer_roots,
    coppersmith_small_root,
    factor_with_hint,
    hint_bits_needed,
    lll,
)


def _det2(m):
    return m[0][0] * m[1][1] - m[0][1] * m[1][0]


def _det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def test_lll_preserves_the_lattice():
    """A reduced basis spans the same lattice: |det| is invariant."""
    rng = random.Random(1)
    for _ in range(15):
        b = [[rng.randrange(-20, 20) for _ in range(3)] for _ in range(3)]
        if _det3(b) == 0:
            continue
        assert abs(_det3(lll(b))) == abs(_det3(b))


def test_lll_shortens():
    rng = random.Random(2)
    for _ in range(15):
        b = [[rng.randrange(-500, 500) for _ in range(3)] for _ in range(3)]
        if _det3(b) == 0:
            continue
        before = min(sum(x * x for x in r) for r in b)
        after = min(sum(x * x for x in r) for r in lll(b))
        assert after <= before


def test_lll_reduced_conditions():
    """Output satisfies the size-reduction condition |mu| <= 1/2."""
    from aksfactor.lattice import _dot

    rng = random.Random(3)
    b = [[rng.randrange(-50, 50) for _ in range(4)] for _ in range(4)]
    if _det2([[b[0][0], b[0][1]], [b[1][0], b[1][1]]]) == 0:
        return
    red = lll(b)
    ortho = []
    for i, row in enumerate(red):
        vec = [Fraction(x) for x in row]
        for j in range(i):
            d = _dot(ortho[j], ortho[j])
            if d:
                mu = Fraction(_dot(red[i], ortho[j])) / d
                assert abs(mu) <= Fraction(1, 2) + Fraction(1, 1000), (i, j, mu)
                vec = [vec[k] - mu * ortho[j][k] for k in range(len(vec))]
        ortho.append(vec)


def test_integer_roots_with_huge_constant_term():
    """The root finder must not enumerate divisors of an N**m-sized constant."""
    c = 123456789012345678901
    poly = [-7 * c, c, -7, 1]  # (x - 7)(x^2 + c)
    assert _integer_roots(poly, 1000) == [7]


def test_integer_roots_edge_cases():
    assert _integer_roots([0, 1], 10) == [0]
    assert _integer_roots([-4, 0, 1], 10) == sorted(_integer_roots([-4, 0, 1], 10))
    assert set(_integer_roots([-4, 0, 1], 10)) == {2, -2}
    assert _integer_roots([1], 10) == []


def test_coppersmith_finds_a_planted_small_root():
    rng = random.Random(5)
    for _ in range(4):
        p = rng.randrange(1 << 15, 1 << 16) | 1
        while not is_prime(p):
            p += 2
        q = rng.randrange(1 << 15, 1 << 16) | 1
        while not is_prime(q):
            q += 2
        if p == q:
            continue
        n = p * q
        unknown = 5
        approx = p & ~((1 << unknown) - 1)
        roots = coppersmith_small_root(n, [approx % n, 1], 1 << unknown, m=3)
        assert (p - approx) in roots, (p, approx, roots)


def test_factor_with_hint():
    rng = random.Random(7)
    ok = 0
    for _ in range(6):
        p = rng.randrange(1 << 13, 1 << 14) | 1
        while not is_prime(p):
            p += 2
        q = rng.randrange(1 << 13, 1 << 14) | 1
        while not is_prime(q):
            q += 2
        if p == q:
            continue
        n = p * q
        unknown = max(1, hint_bits_needed(n) - 3)
        got = factor_with_hint(n, p & ~((1 << unknown) - 1), bound=1 << unknown, m=4)
        if got:
            assert got[0] * got[1] == n
            ok += 1
    assert ok >= 4, ok


def test_hint_bits_needed_is_a_quarter():
    for bits in (32, 64, 128, 256):
        n = 1 << bits
        assert abs(hint_bits_needed(n) - bits / 4) <= 1
