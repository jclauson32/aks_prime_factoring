import random

from aksfactor.arith import is_prime
from aksfactor.collision import (
    conjugate_check,
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
