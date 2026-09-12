"""Round 31: the smoothness wall, measured against Dickman's function.

Every sub-exponential method in this repository is priced by one quantity: how
often a number of a given size has only small prime factors.  This round
measures it -- for random integers, for `p - 1`, and for elliptic group orders
-- against Dickman's `rho`, and reads the three mechanisms off the numbers.
"""

import random
import statistics
from math import gcd, log

from _common import Report

from aksfactor.arith import is_prime, sieve
from aksfactor.smooth import curve_order, dickman_rho, smooth_rate

report = Report(
    "exp38_smoothness",
    "The smoothness wall, measured",
    "Round 31. Dickman's rho against the draws the order mechanism actually makes.",
)

rng = random.Random(38)

report.p("## Dickman's function")
report.p()
known = {1.0: 1.0, 2.0: 0.30685282, 3.0: 0.04860839, 4.0: 0.00491093, 5.0: 0.00035473}
rows = [[u, f"{dickman_rho(u):.6f}", f"{v:.6f}", f"{abs(dickman_rho(u) - v) / v:.2%}"]
        for u, v in known.items()]
report.table(["u", "computed rho(u)", "known", "relative error"], rows)
report.p("`rho(u)` is the density of `x^(1/u)`-smooth numbers near `x`: the "
         "probability that a random number of `u` times the bound's bit length "
         "factors entirely below it. Everything below is measured against it.")
report.p()

report.p("## The three draws")
report.p()
rows = []
for pbits, bound in ((13, 30), (13, 100), (16, 50), (16, 200), (20, 100), (20, 500)):
    lo, hi = 1 << (pbits - 1), 1 << pbits
    primes = []
    while len(primes) < 400:
        v = rng.randrange(lo, hi) | 1
        if is_prime(v):
            primes.append(v)
    randoms = [rng.randrange(lo, hi) for _ in range(400)]
    orders = []
    for p in primes[:120]:
        while True:
            a, x, y = rng.randrange(p), rng.randrange(p), rng.randrange(1, p)
            b = (y * y - x ** 3 - a * x) % p
            if (4 * a ** 3 + 27 * b * b) % p:
                break
        orders.append(curve_order(a, b, p))
    u = log(lo) / log(bound)
    rows.append([f"{pbits}-bit", bound, f"{u:.2f}", f"{dickman_rho(u):.4f}",
                 f"{smooth_rate(randoms, bound):.4f}",
                 f"{smooth_rate([p - 1 for p in primes], bound):.4f}",
                 f"{smooth_rate([p + 1 for p in primes], bound):.4f}",
                 f"{smooth_rate(orders, bound):.4f}"])
report.table(["size", "bound B", "u", "rho(u)", "random m", "p - 1", "p + 1",
              "#E(F_p), random curves"], rows)


def ratio(col):
    return statistics.median(float(r[col]) / max(float(r[4]), 1e-9) for r in rows)


report.p(f"Random integers track `rho(u)` -- that is the baseline. The shifted "
         f"primes beat it: `p - 1` is smooth {ratio(5):.1f} times as often as a "
         f"random number of its size, `p + 1` {ratio(6):.1f} times, because both "
         f"are even and carry small factors more often than a random integer does. "
         f"Elliptic group orders beat it too, by {ratio(7):.1f} times: the Hasse "
         f"interval is full of numbers, and the orders that land in it are biased "
         f"toward small factors in the same way.")
report.p()
report.p("The difference between the mechanisms is not the bias, which is a small "
         "constant, but how many draws each gets. `p - 1` and `p + 1` offer one "
         "number per prime: if it is not smooth, the method is finished. Every "
         "elliptic curve is a new draw from the same biased distribution, so ECM "
         "converts more budget into more draws, and its running time is set by how "
         "many draws are needed -- `1/rho(u)` of them -- against the cost of each. "
         "Optimising that trade-off is exactly the calculation that gives "
         "`L_p[1/2]`, and no amount of redrawing escapes it, because `rho(u)` "
         "falls like `u^(-u)`.")
report.p()

report.p("## What a polynomial-time smoothness method would need")
report.p()
rows = []
for bits in (64, 128, 256, 512):
    for bound_bits in (16, 24):
        u = bits / bound_bits
        rows.append([bits, f"2^{bound_bits}", f"{u:.1f}", f"{dickman_rho(u):.2e}",
                     f"{1 / max(dickman_rho(u), 1e-300):.1e}"])
report.table(["bits of p", "smoothness bound", "u", "rho(u)", "draws needed"], rows)
report.p("A polynomial-time method would need a polylogarithmic bound -- `2^16` or "
         "`2^24` here -- and the table gives the number of draws that costs. "
         "Smooth numbers of that shape are too rare by a factor that grows "
         "super-polynomially in the size of `p`; every known way of making the "
         "numbers smaller (the sieves' `N^(1/2)`, then `N^(2/(d+1))`, then "
         "`L[2/3]`) buys a better exponent inside `L`, never a polynomial.")
report.write()
