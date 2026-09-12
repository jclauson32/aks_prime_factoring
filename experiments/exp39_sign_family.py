"""Round 33: the sign mechanism, five ways.

Every method here manufactures the same object -- a congruence of squares mod
`N` -- and they differ only in how they find one and how big the numbers
involved are.  Two of them use no smooth numbers at all and cost `N^(1/4)`;
three use smooth numbers and are sub-exponential.  This is the mechanism's
whole story on one page.
"""

import random
import statistics
import time
from math import isqrt, log2

from _common import Report

from aksfactor.arith import is_prime
from aksfactor.central import central_split
from aksfactor.cfrac import cf_residues, cfrac, factor_base
from aksfactor.nfs import number_field_sieve
from aksfactor.qs import default_parameters, quadratic_sieve, sqrt_mod_prime
from aksfactor.squfof import squfof

report = Report(
    "exp39_sign_family",
    "The sign mechanism, five ways",
    "Round 33. One mechanism, five relation-finders, and what each one buys.",
)

rng = random.Random(39)


def randprime(bits):
    while True:
        v = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if is_prime(v):
            return v


def timed(fn, n):
    t0 = time.perf_counter()
    out = fn(n)
    return out, time.perf_counter() - t0


def first(x):
    return x[0] if isinstance(x, tuple) else x


report.p("## The five")
report.p()
methods = [
    ("central column (round 15)", lambda n: first(central_split(n)), 50),
    ("SQUFOF", squfof, 90),
    ("continued fractions (CFRAC)", cfrac, 90),
    ("quadratic sieve", quadratic_sieve, 90),
    ("number field sieve", lambda n: number_field_sieve(n, 3, 1500, 1500, 24, 8000, 3000), 70),
]
sizes = (50, 60, 70, 80, 90)
numbers = {}
for bits in sizes:
    numbers[bits] = []
    while len(numbers[bits]) < 3:
        p, q = randprime(bits // 2), randprime(bits // 2)
        if p != q:
            numbers[bits].append((p * q, p, q))

rows = []
for name, fn, top in methods:
    row = [name]
    for bits in sizes:
        if bits > top:
            row.append("--")
            continue
        times, ok = [], 0
        for n, p, q in numbers[bits]:
            g, dt = timed(fn, n)
            times.append(dt)
            ok += g in (p, q)
        row.append(f"{statistics.median(times):.2f} s" + ("" if ok == 3 else f" ({ok}/3)"))
    rows.append(row)
report.table(["method"] + [f"{b} bits" for b in sizes], rows)
report.p("`--` means the method was not run at that size (the central column is "
         "`N^(1/4)` with a large constant; the toy number field sieve is a "
         "degree-3 implementation).")
report.p()

report.p("## The numbers that have to be smooth")
report.p()
rows = []
for bits in (60, 80):
    n = numbers[bits][0][0]
    cf = [q for _, q, _ in cf_residues(n, 4000)]
    cf_bits = statistics.median(log2(q) for q in cf if q > 1)
    b, width = default_parameters(n)
    m = isqrt(n) + 1
    qs_vals = [abs((x + m) * (x + m) - n) for x in range(1, width)]
    qs_bits = statistics.median(log2(v) for v in qs_vals if v > 1)
    rows.append([bits, f"{log2(n) / 2:.0f}", f"{cf_bits:.0f}", f"{qs_bits:.0f}",
                 f"{len(factor_base(n, b))}"])
report.table(["bits of N", "sqrt N", "CFRAC residue (median)",
              "QS residue over its sieve interval (median)", "factor base"], rows)
report.p("CFRAC's residues are the smallest of any method here -- below `2 sqrt N` "
         "by construction -- and smaller than the sieve's, which grow with the "
         "distance along the sieve interval. Smaller numbers are smooth more "
         "often, so CFRAC needs fewer candidates. It is still the slower method "
         "asymptotically, because the sieve finds its smooth values in bulk: one "
         "pass over an interval marks every candidate divisible by each factor "
         "base prime, while CFRAC must trial-divide every residue it generates.")
report.p()

report.p("## What each change buys")
report.p()
report.p("| method | how the congruence is found | smooth numbers? | cost |")
report.p("|---|---|---|---|")
report.p("| central column | Legendre symbol from a truncated hypergeometric sum | no | `N^(1/4)` |")
report.p("| SQUFOF | walk the cycle of forms of discriminant `4N` to a square form | no | `N^(1/4)` |")
report.p("| CFRAC | convergents of `sqrt N`, trial-divided | yes, residues `< 2 sqrt N` | `L[1/2]` |")
report.p("| quadratic sieve | sieve `(x+m)^2 - N` over an interval | yes, residues `~ sqrt N` times the interval | `L[1/2]`, faster constant |")
report.p("| number field sieve | two rings, `a - bm` and `N(a - b alpha)` | yes, two numbers, `L[2/3]` in total | `L[1/3]` |")
report.p()
report.p("The mechanism is constant across the table; only the relation-finder "
         "changes. Without smooth numbers the mechanism costs `N^(1/4)` -- which "
         "is exactly where Pascal's central column landed in round 15, and it is "
         "not a weakness of the triangle. With smooth numbers it becomes "
         "sub-exponential, and every subsequent improvement -- sieving instead of "
         "trial division, two rings instead of one -- makes the numbers that must "
         "be smooth smaller or the search for them faster. None of it makes the "
         "method polynomial, because `rho(u)` (round 31) falls too fast.")
report.write()
