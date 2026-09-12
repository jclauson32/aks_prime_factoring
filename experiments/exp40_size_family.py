"""Round 34: the size mechanism, several ways -- measured as coverage per operation.

The size mechanism rules out candidate factors.  Trial division rules out one
per operation.  Pollard and Strassen's product tree rules out a whole batch per
operation, Lehman's geometry rules out an interval, and Coppersmith's lattice
rules out everything at once given enough bits of `p`.  This round measures the
exchange rate: candidates eliminated per unit of work.
"""

import random
import statistics
import time
from math import isqrt, log2

from _common import Report

from aksfactor.arith import is_prime, sieve
from aksfactor.fast import factorial_mod, fast_spf
from aksfactor.hyperbola import hull_walk
from aksfactor.lattice import cost_exponents, hint_bits_needed

report = Report(
    "exp40_size_family",
    "The size mechanism, by coverage per operation",
    "Round 34. How many candidate factors each method rules out per unit of work.",
)

rng = random.Random(40)


def randprime(bits):
    while True:
        v = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if is_prime(v):
            return v


report.p("## Coverage per second")
report.p()
report.p("Each method is given a range of candidates to rule out and timed on a "
         "number with no factor in that range, so nothing stops early. The "
         "measure is candidates ruled out per second.")
report.p()
rows = []
for bits in (48, 60):
    p, q = randprime(bits // 2), randprime(bits // 2)
    n = p * q
    limit = 1 << (bits // 2 - 2)              # below both primes: no early exit

    t0 = time.perf_counter()
    count = 0
    for d in sieve(min(limit, 200000)):
        count += 1
        if n % d == 0:
            break
    trial_rate = count / (time.perf_counter() - t0)

    t0 = time.perf_counter()
    factorial_mod(limit, n)
    strassen_rate = limit / (time.perf_counter() - t0)

    t0 = time.perf_counter()
    steps = 0
    top = isqrt(n - 1)
    stop = int(top * 0.97)
    for x, y, _ in hull_walk(n - 1, top):
        steps += 1
        if y < stop:
            break
    rows_covered = top - stop
    hull_rate = rows_covered / (time.perf_counter() - t0)

    rows.append([bits, f"{trial_rate:,.0f}", f"{strassen_rate:,.0f}", f"{hull_rate:,.0f}"])
report.table(["bits of N", "trial division (candidates/s)",
              "product tree, factorial mod n (candidates/s)",
              "hyperbola hull walk (rows/s)"], rows)
report.p("The product tree wins by batching: one modular multiplication absorbs a "
         "whole block of candidates, so its rate per second is far above trial "
         "division's even though both are `Theta(p)` in candidates. The hull walk "
         "covers rows of the hyperbola rather than candidate divisors -- each "
         "vertex settles a stretch of about `N^(1/6)` rows -- which is why it "
         "reaches `N^(1/3)` while the others are `N^(1/2)`.")
report.p()

report.p("## What each change buys")
report.p()
report.p("| method | what one operation covers | cost |")
report.p("|---|---|---|")
report.p("| trial division, Pascal row scan | one candidate | `Theta(p)` |")
report.p("| Pollard-Strassen product tree (T20) | a batch, by multipoint evaluation | `O~(N^(1/4))` |")
report.p("| Lehman / the hull walk (P23) | an interval of the hyperbola near a rational slope | `O~(N^(1/3))` |")
report.p("| Harvey's method | Lehman's intervals plus a Fermat congruence and a BSGS sweep | `N^(1/5)` |")
report.p("| Coppersmith with a hint | everything, once `N^(1/4)` of `p` is known | polynomial, given the bits |")
report.p()
bits_needed = hint_bits_needed(1 << 256)
report.p(f"The last row is the one that matters for a breakthrough, and round 12 "
         f"priced it: Coppersmith needs about a quarter of `p`'s bits "
         f"({bits_needed} of them for a 256-bit `N`), and guessing them costs "
         f"`N^(1/4)` -- exactly what the other size methods already cost. The "
         f"exponent `beta^2` in Coppersmith's bound is what would have to move.")
report.p()
report.table(["beta (p ~ N^beta)", "Coppersmith window", "guess + Coppersmith", "Strassen"],
             [[f"{b:.2f}"] + [f"N^{cost_exponents(b)[k]:.3f}" for k in
                              ("coppersmith_window", "guess_and_coppersmith", "strassen")]
              for b in (0.5, 0.4, 0.3)])
report.p("Every size method in this repository is a different way of paying the "
         "same `N^(1/4)`-to-`N^(1/3)` bill for locating `p` within a window "
         "Coppersmith can finish from. None of them is cheaper than the batching "
         "bound, and the one that is faster -- Harvey's `N^(1/5)` -- gets there by "
         "combining the geometry with a congruence, not by covering more "
         "candidates per operation.")
report.write()
