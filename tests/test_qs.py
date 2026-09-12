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


def test_auto_factor_command():
    from aksfactor.cli import _auto_factor

    for n in (1234567, 2 ** 20 + 1, 1000003 * 999983, 101 ** 3 * 103, 6):
        facs = _auto_factor(n, verbose=lambda *a: None)
        product = 1
        for p, e in facs.items():
            product *= p ** e
            assert is_prime(p)
        assert product == n


def test_package_exports_resolve():
    import aksfactor

    missing = [name for name in aksfactor.__all__ if not hasattr(aksfactor, name)]
    assert not missing, missing
    n = 1000003 * 999983
    assert aksfactor.squfof(n) in (1000003, 999983)
    assert aksfactor.cfrac(n) in (1000003, 999983)
    assert aksfactor.quadratic_sieve(n) in (1000003, 999983)
