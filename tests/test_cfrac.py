import random
from math import isqrt

from aksfactor.arith import is_prime
from aksfactor.cfrac import cf_residues, cfrac


def test_cf_residues_satisfy_the_congruence():
    n = 1000003 * 999983
    for a, q, sign in list(cf_residues(n, 60)):
        assert q < 2 * isqrt(n) + 2
        assert (a * a - sign * q) % n == 0


def test_cfrac_factors():
    rng = random.Random(4)
    for bits in (32, 48, 64):
        for _ in range(2):
            while True:
                p = rng.randrange(1 << (bits // 2 - 1), 1 << (bits // 2)) | 1
                if is_prime(p):
                    break
            while True:
                q = rng.randrange(1 << (bits // 2 - 1), 1 << (bits // 2)) | 1
                if is_prime(q) and q != p:
                    break
            assert cfrac(p * q) in (p, q)
