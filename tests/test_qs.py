import random

from aksfactor.arith import is_prime
from aksfactor.qs import quadratic_sieve, sqrt_mod_prime


def test_tonelli_shanks():
    for p in (3, 5, 13, 17, 97, 257, 65537, 1000003):
        for a in range(1, min(p, 200)):
            if pow(a, (p - 1) // 2, p) == 1:
                r = sqrt_mod_prime(a, p)
                assert r * r % p == a


def test_quadratic_sieve_factors():
    rng = random.Random(7)
    for bits in (36, 48, 60):
        for _ in range(2):
            while True:
                p = rng.randrange(1 << (bits // 2 - 1), 1 << (bits // 2)) | 1
                if is_prime(p):
                    break
            while True:
                q = rng.randrange(1 << (bits // 2 - 1), 1 << (bits // 2)) | 1
                if is_prime(q) and q != p:
                    break
            assert quadratic_sieve(p * q) in (p, q)
