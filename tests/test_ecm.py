import random

from aksfactor.arith import is_prime
from aksfactor.ecm import ecm, ladder, pminus1, suyama_curve


def _count(a, p):
    """#E(F_p) for y^2 = x^3 + a x^2 + x."""
    total = 1
    for x in range(p):
        r = (x ** 3 + a * x * x + x) % p
        total += 1 if r == 0 else (2 if pow(r, (p - 1) // 2, p) == 1 else 0)
    return total


def test_ladder_respects_the_group_order():
    p = 10007
    checked = 0
    for sigma in range(6, 40):
        a24, x, z = suyama_curve(sigma, p)
        a = (4 * a24 - 2) % p
        order = _count(a, p)
        twist = 2 * p + 2 - order
        # the x-line carries the curve and its quadratic twist
        assert ladder(order, x, z, a24, p)[1] % p == 0 or ladder(twist, x, z, a24, p)[1] % p == 0
        assert order % 12 == 0 or twist % 12 == 0
        checked += 1
    assert checked > 20


def test_ecm_factors():
    rng = random.Random(4)
    for pb, qb in ((20, 40), (26, 60)):
        while True:
            p = rng.randrange(1 << (pb - 1), 1 << pb) | 1
            if is_prime(p):
                break
        while True:
            q = rng.randrange(1 << (qb - 1), 1 << qb) | 1
            if is_prime(q):
                break
        assert ecm(p * q, b1=1000, b2=50000, curves=200) in (p, q)


def test_pminus1_needs_a_smooth_p_minus_1():
    p = next(k * 2 * 3 * 5 * 7 * 11 * 13 + 1 for k in range(1, 1000)
             if is_prime(k * 2 * 3 * 5 * 7 * 11 * 13 + 1) and k > 50)
    q = 1000000007                                                  # q - 1 = 2 * 500000003
    assert pminus1(p * q, 100) == p


def test_pminus1_stage_two_catches_one_large_prime():
    from aksfactor.ecm import is_b1b2_smooth

    # p - 1 = 2 * 3 * 5 * 7 * 11 * L with one prime L in (B1, B2]
    for big in (1009, 2003, 4001, 7919):
        p = 2 * 3 * 5 * 7 * 11 * big + 1
        if not is_prime(p):
            continue
        assert is_b1b2_smooth(p - 1, 100, 10000)
        assert pminus1(p * 1000000007, 100) is None
        assert pminus1(p * 1000000007, 100, 10000) == p
