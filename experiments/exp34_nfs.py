"""Round 26: a toy number field sieve -- the sign mechanism at L[1/3].

Round 21's quadratic sieve needs numbers of size about sqrt(N) times the sieve
length to be smooth.  The number field sieve splits the job between two rings
and asks two numbers of size about N^(1/(d+1)) each to be smooth.  This is a
complete, working, pure-Python toy: base-m polynomial selection ranked by root
properties, a two-sided line sieve, quadratic characters, GF(2) algebra, and
the square root in Z[alpha] by p-adic Newton lifting.
"""

import math
import random
import statistics
import time
from math import isqrt, log2

from _common import Report

from aksfactor.arith import is_prime
from aksfactor.nfs import (
    collect_relations,
    factor_bases,
    norm,
    number_field_sieve,
    select_polynomial,
)
from aksfactor.qs import quadratic_sieve

report = Report(
    "exp34_nfs",
    "A toy number field sieve",
    "Round 26. The sign mechanism with the smallest numbers anyone knows how to make smooth.",
)

rng = random.Random(34)

PARAMS = {40: (400, 3000), 50: (700, 5000), 60: (1500, 8000), 70: (2500, 12000),
          80: (4000, 16000)}


def randprime(bits):
    while True:
        v = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if is_prime(v):
            return v


report.p("## 1. It works")
report.p()
rows = []
for bits, (bound, width) in PARAMS.items():
    times, qs_times, ok = [], [], 0
    reps = 3 if bits <= 70 else 2
    for _ in range(reps):
        p, q = randprime(bits // 2), randprime(bits // 2)
        while q == p:
            q = randprime(bits // 2)
        n = p * q
        st = {}
        t0 = time.perf_counter()
        g = number_field_sieve(n, 3, bound, bound, 24, width, 3000, st)
        times.append(time.perf_counter() - t0)
        ok += g in (p, q)
        t0 = time.perf_counter()
        assert quadratic_sieve(n) in (p, q)
        qs_times.append(time.perf_counter() - t0)
    rows.append([bits, f"{ok} of {reps}", f"{statistics.median(times):.1f} s",
                 f"{statistics.median(qs_times):.2f} s"])
report.table(["bits of N", "NFS factored", "NFS time (median)", "QS time (median)"], rows)
report.p("Every dependency that survives the rational check goes through the "
         "number-field square root, `beta^2 = f'(alpha)^2 prod (a - b alpha)` in "
         "`Z[alpha]`, lifted from an inert prime; `gcd(phi(beta) - f'(m) X, N)` then "
         "splits `N`. At these sizes the quadratic sieve is faster, as it is in "
         "practice until about a hundred decimal digits.")
report.p()

report.p("## 2. Why it wins eventually: the sizes that must be smooth")
report.p()
rows = []
gaps = []
for bits in (60, 80, 100, 120):
    p, q = randprime(bits // 2), randprime(bits // 2)
    n = p * q
    m, f = select_polynomial(n, 3, tries=12)[0]
    rat, alg, _ = factor_bases(f, 1500, 1500, 0)
    rels = collect_relations(n, m, f, rat, alg, 40, 8000, 3000)
    if not rels:
        continue
    nfs_bits = statistics.median(log2(abs(a - b * m)) + log2(abs(norm(f, a, b)))
                                 for a, b, _, _ in rels)
    qs_bits = log2(2 * 20000 * isqrt(n))
    rows.append([bits, f"{qs_bits:.0f}", f"{nfs_bits:.0f}", f"{nfs_bits - qs_bits:+.0f}"])
    gaps.append(nfs_bits - qs_bits)
report.table(["bits of N", "QS: bits of |Q(x)|, nominal sieve length 2 x 10^4",
              "NFS: bits of |a - bm| |N(a - b alpha)|, median of relations found",
              "NFS minus QS"], rows)
report.p(f"At degree 3 and these sizes the two numbers the NFS needs smooth are, "
         f"together, {min(gaps):.0f} to {max(gaps):.0f} bits *larger* than the one "
         f"the quadratic sieve needs -- so it finds relations more slowly and loses, "
         f"as the timings show. Its advantage is asymptotic: with the degree raised "
         f"as `d ~ (3 log N / log log N)^(1/3)`, the product shrinks to "
         f"`L[2/3]` in size against the sieve's `N^(1/2 + o(1))`, and the smoothness "
         f"probability of `L[2/3]`-sized numbers gives the `L[1/3]` running time. In "
         f"practice the crossover is near a hundred decimal digits.")
report.p()
report.p("That completes the list. Each of round 14's four mechanisms now has its "
         "best method implemented here and run on the repository's own numbers: "
         "ECM for order, Harvey for size, rho for collision, and the quadratic and "
         "number field sieves for sign. The fastest of them, the number field "
         "sieve, is sub-exponential and not polynomial, because smooth numbers of "
         "size `L[2/3]` are only `L[1/3]`-rare.")
report.write()
