import random

from aksfactor.arith import is_prime
from aksfactor.squfof import squfof


def test_squfof_factors_semiprimes():
    rng = random.Random(3)
    for bits in (24, 40, 56, 72):
        for _ in range(3):
            while True:
                p = rng.randrange(1 << (bits // 2 - 1), 1 << (bits // 2)) | 1
                if is_prime(p):
                    break
            while True:
                q = rng.randrange(1 << (bits // 2 - 1), 1 << (bits // 2)) | 1
                if is_prime(q) and q != p:
                    break
            assert squfof(p * q) in (p, q)


def test_squfof_handles_squares_and_even_numbers():
    assert squfof(2 * 7919) == 2
    assert squfof(7919 ** 2) == 7919
