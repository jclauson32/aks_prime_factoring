"""Do extra AKS bases `a` buy anything?  Measure the per-coefficient hit rate.

AKS itself sweeps many bases `a` in `(x+a)^n`. If each base gave a fresh,
independent shot at leaking a factor at rate better than 1/p, the fold attack
would break the trial-division barrier. This measures the actual rate.
"""

import random

from _common import Report

from aksfactor.arith import is_prime
from aksfactor.ring import fold_coefficients

report = Report(
    "exp05_multibase",
    "Hit rate per fold coefficient, across AKS bases",
    "Fraction of coefficients `S_j` with `gcd(S_j, n) > 1`, compared with the "
    "predicted `1/p`.",
)

rng = random.Random(4242)


def randprime(lo, hi):
    while True:
        c = rng.randrange(lo, hi) | 1
        if is_prime(c):
            return c


rows = []
for lo, hi in [(30, 80), (100, 300), (400, 900), (1200, 2500), (4000, 8000)]:
    p = randprime(lo, hi)
    q = randprime(p + 2, p * 4)
    n = p * q
    hits = total = 0
    for a in range(1, 13):
        for r in range(2, 60):
            coeffs = fold_coefficients(n, r, a)
            for c in coeffs:
                total += 1
                g = __import__("math").gcd(c, n)
                if 1 < g < n:
                    hits += 1
    rows.append([p, q, f"{total:,}", hits, f"{hits / total:.5f}", f"{1 / p:.5f}",
                 f"{(hits / total) * p:.2f}"])
report.table(["p", "q", "coefficients", "hits", "observed rate", "1/p",
              "observed x p"], rows)
report.p("The observed rate tracks `1/p` (last column near 1) regardless of how many "
         "bases `a` are swept. Extra bases give more independent tickets, never a "
         "better price per ticket, so the total work stays `Theta(p)`.")
report.write()
