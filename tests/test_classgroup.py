"""Class groups of imaginary quadratic orders: group law, and factoring."""

import random
import sys
from math import gcd, isqrt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from aksfactor.classgroup import (
    ambiguous_factor,
    class_number_bruteforce,
    compose,
    identity_form,
    pow_form,
    reduce_form,
    schnorr_lenstra_split,
    square_form,
)


def _reduced_forms(d):
    out = []
    for a in range(1, isqrt(-d // 3) + 2):
        for b in range(-a + 1, a + 1):
            num = b * b - d
            if num % (4 * a):
                continue
            c = num // (4 * a)
            if c < a or gcd(gcd(a, b), c) != 1 or (a == c and b < 0):
                continue
            out.append((a, b, c))
    return out


DISCS = [d for d in range(-3, -600, -1) if d % 4 in (0, 1)]


def test_reduce_is_idempotent_and_preserves_discriminant():
    rng = random.Random(1)
    for _ in range(400):
        a = rng.randrange(1, 60)
        b = rng.randrange(-80, 80)
        c = rng.randrange(1, 60)
        d = b * b - 4 * a * c
        if d >= 0:
            continue
        f = reduce_form(a, b, c)
        assert f[1] ** 2 - 4 * f[0] * f[2] == d
        assert reduce_form(*f) == f
        assert -f[0] < f[1] <= f[0] <= f[2]


def test_class_number_matches_enumeration():
    for d in DISCS[:80]:
        assert class_number_bruteforce(d) == len(_reduced_forms(d)), d


def test_identity_and_inverse():
    for d in DISCS[:80]:
        idf = identity_form(d)
        for f in _reduced_forms(d):
            assert compose(f, idf, d) == f, (d, f)
            inv = reduce_form(f[0], -f[1], f[2])
            assert compose(f, inv, d) == idf, (d, f)


def test_order_divides_class_number():
    for d in DISCS[:80]:
        h = class_number_bruteforce(d)
        idf = identity_form(d)
        for f in _reduced_forms(d):
            assert pow_form(f, h, d) == idf, (d, f, h)


def test_associativity():
    rng = random.Random(7)
    for d in DISCS[:80]:
        forms = _reduced_forms(d)
        for _ in range(12):
            x, y, z = (rng.choice(forms) for _ in range(3))
            assert compose(compose(x, y, d), z, d) == compose(x, compose(y, z, d), d)


def test_square_matches_compose():
    for d in DISCS[:40]:
        for f in _reduced_forms(d):
            assert square_form(f, d) == compose(f, f, d)


def test_improper_equivalence_regression():
    """An xgcd returning a negative gcd gives det = -1 -- an improper
    equivalence, which silently computes in the wrong class."""
    d = -23
    f = (2, -1, 3)
    inv = (2, 1, 3)
    assert compose(f, inv, d) == identity_form(d)


def test_ambiguous_factor():
    # b == 0 case: D = -4ac
    assert ambiguous_factor((2, 0, 13), -104, 26) in (2, 13)


def test_schnorr_lenstra_splits():
    for n in (91, 143, 187, 209, 221, 391, 667, 1147, 2021, 5183, 1009 * 1013):
        got = schnorr_lenstra_split(n, multipliers=range(1, 30), bound=200,
                                    forms_per_disc=3)
        assert got is not None, n
        assert 1 < got[0] < n and n % got[0] == 0, (n, got)
