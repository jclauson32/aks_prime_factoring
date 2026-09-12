"""Round 31: the smoothness wall, measured against Dickman's function.

Every sub-exponential method in this repository is priced by one quantity: how
often a number of a given size has only small prime factors.  This round
measures it -- for random integers, for `p - 1`, and for elliptic group orders
-- against Dickman's `rho`, and reads the three mechanisms off the numbers.
"""

import random
import statistics
from math import isqrt, log, sqrt

from _common import Report

from aksfactor.arith import is_prime, sieve
from aksfactor.smooth import (
    dickman_rho,
    qr_table,
    smooth_rate,
    weierstrass_order,
)

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
def se(rate, n):
    return sqrt(max(rate * (1 - rate), 1e-12) / n)


def uniform_orders(primes, per_prime):
    """`#E` for curves drawn uniformly among the nonsingular `(a, b)`.

    Drawing a random *point* and solving for `b` -- which is what this round did
    first -- weights each curve by how many points it has, and the weighting is
    exactly what the table is trying to measure.
    """
    out = []
    for p in primes:
        table = qr_table(p)
        for _ in range(per_prime):
            while True:
                a, b = rng.randrange(p), rng.randrange(p)
                if (4 * a ** 3 + 27 * b * b) % p:
                    break
            out.append(weierstrass_order(a, b, p, table))
    return out


rows = []
for pbits, bound in ((13, 30), (13, 100), (16, 50), (16, 200), (20, 100), (20, 500)):
    lo, hi = 1 << (pbits - 1), 1 << pbits
    primes = []
    while len(primes) < 400:
        v = rng.randrange(lo, hi) | 1
        if is_prime(v):
            primes.append(v)
    randoms = [rng.randrange(lo, hi) for _ in range(400)]
    orders = uniform_orders(primes[:30 if pbits <= 16 else 12], 4)
    u = log(lo) / log(bound)
    rate_rand = smooth_rate(randoms, bound)
    rate_curve = smooth_rate(orders, bound)
    rows.append([f"{pbits}-bit", bound, f"{u:.2f}", f"{dickman_rho(u):.4f}",
                 f"{rate_rand:.4f} +- {se(rate_rand, len(randoms)):.4f}",
                 f"{smooth_rate([p - 1 for p in primes], bound):.4f}",
                 f"{smooth_rate([p + 1 for p in primes], bound):.4f}",
                 f"{rate_curve:.4f} +- {se(rate_curve, len(orders)):.4f} "
                 f"({len(orders)})"])
report.table(["size", "bound B", "u", "rho(u)", "random m (400)", "p - 1 (400)",
              "p + 1 (400)", "#E(F_p), uniform curves (samples)"], rows)


def ratio(col):
    return statistics.median(float(r[col].split(" +-")[0]) / max(float(r[4].split(" +-")[0]), 1e-9)
                             for r in rows)


report.p(f"Random integers track `rho(u)` -- that is the baseline. The shifted "
         f"primes beat it: `p - 1` is smooth {ratio(5):.1f} times as often as a "
         f"random number of its size, `p + 1` {ratio(6):.1f} times, because both "
         f"are even and carry small factors more often than a random integer "
         f"does. The curve column's median ratio is {ratio(7):.1f}, but at a few "
         f"dozen orders per cell that column cannot carry a claim, so the next block "
         f"settles it at one size with enough samples to have an error bar worth "
         f"reading.")
report.p()

report.p("### Are curve orders smoother than random integers?")
report.p()
settle_bits, settle_bound = 16, 200
lo, hi = 1 << (settle_bits - 1), 1 << settle_bits
settle_primes = []
while len(settle_primes) < 1500:
    v = rng.randrange(lo, hi) | 1
    if is_prime(v):
        settle_primes.append(v)
settle_orders = uniform_orders(settle_primes[:150], 8)
settle_random = [rng.randrange(lo, hi) for _ in range(1500)]
hasse = [rng.randrange(p + 1 - 2 * isqrt(p), p + 2 + 2 * isqrt(p))
         for p in settle_primes[:1200]]

settle = []
for name, sample in (("random integer of the same bit length", settle_random),
                     ("random integer from the Hasse interval", hasse),
                     ("`#E`, uniform curves", settle_orders),
                     ("`p - 1`", [p - 1 for p in settle_primes]),
                     ("`p + 1`", [p + 1 for p in settle_primes])):
    rate = smooth_rate(sample, settle_bound)
    settle.append((name, rate, se(rate, len(sample)), len(sample)))
report.table([f"what is tested for {settle_bound}-smoothness", "rate", "samples"],
             [[name, f"{rate:.4f} +- {err:.4f}", n] for name, rate, err, n in settle])

base = settle[0]
curve = settle[2]
shifted = settle[3]
gap_curve = (curve[1] - base[1]) / sqrt(curve[2] ** 2 + base[2] ** 2)
gap_shift = (shifted[1] - base[1]) / sqrt(shifted[2] ** 2 + base[2] ** 2)
report.p(f"`p - 1` is {shifted[1] / base[1]:.2f} times the baseline "
         f"({gap_shift:+.1f} sigma) -- a real bias, and the one the round's "
         f"argument uses. Uniform curve orders are {curve[1] / base[1]:.2f} times "
         f"the baseline ({gap_curve:+.1f} sigma), "
         + ("which is a real effect at this sample size."
            if abs(gap_curve) >= 2 else
            "which this sample cannot separate from the baseline: a curve order "
            "behaves like a random integer of its size, exactly as ECM's standard "
            "heuristic assumes. An earlier version of this round reported a "
            "factor of 1.5 for this row from 120 samples drawn by picking a "
            "random point -- too few, and weighted by the curve's own point "
            "count.") +
         " The mechanism's economics do not depend on which it is: what ECM buys "
         "is the *redraw*, not a better draw (round 35).")
report.p()
report.p("The difference between the mechanisms is not the bias, which is a small "
         "constant, but how many draws each gets. `p - 1` and `p + 1` offer one "
         "number per prime: if it is not smooth, the method is finished. Every "
         "elliptic curve is a new draw from the same distribution, so ECM "
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
