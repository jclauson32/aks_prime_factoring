import random

from aksfactor.arith import is_prime
from aksfactor.collision import (
    binom_poly,
    conjugate_check,
    degenerate_start,
    fibre_kappa,
    pascal_rho,
    predicted_kappa,
)


def test_choose_two_is_pollard_rho():
    for p in (101, 1009, 65537):
        assert conjugate_check(p)


def test_fibre_statistic_follows_the_reflection():
    for k in range(2, 9):
        got = fibre_kappa(k, 10007)
        assert abs(got - predicted_kappa(k)) < 0.1, (k, got)


def test_small_starts_are_fixed_points_not_walks():
    """``0`` and ``k + 1`` are fixed, and starts below ``k`` fall into ``0``."""
    n = 1000003 * 1000033
    for k in (2, 3, 4, 5):
        f = binom_poly(k, n)
        assert f(0) == 0 and f(k + 1) == (k + 1) % n
        for x0 in range(0, k + 2):
            assert degenerate_start(k, x0)
            factor, steps = pascal_rho(n, k, x0)
            assert factor is None and steps < 10, (k, x0, factor, steps)
        assert not degenerate_start(k, k + 2)
    assert degenerate_start(1, 10 ** 6)      # C(x, 1) = x, nothing ever moves


def test_the_default_start_walks():
    """The default must be a start that can actually split something."""
    n = 1000003 * 1000033
    for k in (2, 3, 4, 5):
        g, steps = pascal_rho(n, k)
        assert g in (1000003, 1000033), (k, g)
        assert steps > 0


def test_pascal_rho_factors():
    rng = random.Random(4)
    for k in (2, 3, 4):
        for _ in range(4):
            while True:
                p = rng.randrange(10 ** 5, 10 ** 6) | 1
                if is_prime(p):
                    break
            while True:
                q = rng.randrange(10 ** 5, 10 ** 6) | 1
                if is_prime(q) and q != p:
                    break
            n = p * q
            for x0 in range(3, 30):
                g, _ = pascal_rho(n, k, x0)
                if g:
                    break
            assert g in (p, q)
