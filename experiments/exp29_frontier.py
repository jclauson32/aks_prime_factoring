"""Round 21: the frontier, measured in one codebase.

Twenty rounds put Pascal's triangle on all four mechanisms and found each at
its baseline -- N^(1/4) at best.  This experiment runs a small quadratic sieve
(the sign mechanism with smoothness behind it) against the repository's
N^(1/4) methods on the same semiprimes, to show where the curves cross and what
the crossing is made of.
"""

import math
import random
import time

from _common import Report

from aksfactor.arith import is_prime
from aksfactor.collision import pascal_rho
from aksfactor.fast import fast_spf
from aksfactor.qs import quadratic_sieve

report = Report(
    "exp29_frontier",
    "The frontier: a quadratic sieve against the N^(1/4) family",
    "Round 21. Where smoothness overtakes every Pascal-native method.",
)

rng = random.Random(29)


def randprime(bits):
    while True:
        x = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if is_prime(x):
            return x


def timed(fn):
    t0 = time.perf_counter()
    out = fn()
    return out, time.perf_counter() - t0


def run_rho(n):
    for x0 in range(3, 200):
        g, _ = pascal_rho(n, 2, x0)
        if g:
            return g
    return None


methods = [
    ("quadratic sieve (sign + smoothness)", quadratic_sieve, 5, 110),
    ("Pascal rho = Pollard rho (collision)", run_rho, 5, 100),
    ("factorial threshold, Strassen (size)", lambda n: fast_spf(n), 1, 60),
]
SIZES = (50, 60, 70, 80, 90, 100, 110)
medians = {name: {} for name, *_ in methods}
for bits in SIZES:
    numbers = []
    for _ in range(5):
        p, q = randprime(bits // 2), randprime(bits // 2)
        while q == p:
            q = randprime(bits // 2)
        numbers.append((p * q, p, q))
    for name, fn, reps, top in methods:
        if bits > top:
            continue
        times = []
        for n, p, q in numbers[:reps]:
            g, dt = timed(lambda: fn(n))
            assert g in (p, q), (name, bits, g)
            times.append(dt)
        times.sort()
        medians[name][bits] = times[len(times) // 2]

rows = []
for bits in SIZES:
    row = [bits]
    for name, _, reps, top in methods:
        t = medians[name].get(bits)
        row.append(f"{t:.2f} s" if t is not None else "--")
    rows.append(row)
report.p("Balanced semiprimes: five per size for the sieve and rho (median time "
         "shown, same numbers for both), one for Strassen's method, whose cost "
         "grows fastest in pure Python. `--` marks sizes not attempted.")
report.p()
report.table(["bits of N"] + [f"{name} ({reps} each)" for name, _, reps, _ in methods], rows)


def fit(name, lo=0.05):
    pts = [(b, math.log(t)) for b, t in medians[name].items() if t >= lo]
    if len(pts) < 2:
        return None
    mb = sum(b for b, _ in pts) / len(pts)
    mt = sum(t for _, t in pts) / len(pts)
    slope = (sum((b - mb) * (t - mt) for b, t in pts) /
             sum((b - mb) ** 2 for b, _ in pts))
    return math.exp(10 * slope), pts[0][0], pts[-1][0]


fits = {name: fit(name) for name, *_ in methods}
report.p("Least-squares growth of the median time per 10 bits of `N` (sizes "
         "taking at least 0.05 s):")
report.p()
report.table(["method", "x per 10 bits", "over", "N^(1/4) would be"],
             [[name, f"{f[0]:.1f}" if f else "--", f"{f[1]}-{f[2]} bits" if f else "--",
               "5.7"] for name, f in fits.items()])
qs_name, rho_name = methods[0][0], methods[1][0]
cross = next((b for b in SIZES if b in medians[rho_name]
              and all(medians[qs_name][c] < medians[rho_name][c]
                      for c in SIZES if c >= b and c in medians[rho_name])), None)
qs_top = max(medians[qs_name])
report.p(f"The sieve is faster than rho from {cross} bits on and factors the "
         f"{qs_top}-bit numbers in a median {medians[qs_name][qs_top]:.1f} s of pure "
         f"Python. Rho grows by {fits[rho_name][0]:.1f}x per ten bits "
         f"(`N^(1/4)` predicts 5.7); the sieve by {fits[qs_name][0]:.1f}x over "
         f"the same kind of range. That is the sub-exponential `L[1/2]` at work: "
         f"in theory its growth factor per ten bits keeps shrinking, and at these "
         f"sizes it is already the smaller one.")
report.p()
report.p("What the sieve has that the central column (round 15) lacks is the "
         "whole difference. Both manufacture a zero divisor from the *sign* of a "
         "square root. The central column asks for one Legendre symbol and pays "
         "`sqrt p` for it. The sieve never evaluates a Legendre symbol mod `p` at "
         "all: it collects `x^2 = Q(x) (mod N)` with `Q(x)` small enough -- about "
         "`sqrt N` times the sieve length -- to be smooth unusually often, and "
         "lets linear algebra over `F_2` find a product that is a square on both "
         "sides. Its cost is the density of smooth numbers, the same quantity "
         "that caps the order methods at `L[1/2]` (Proposition 17), and the "
         "number field sieve pushes it to `L[1/3]` by making the numbers smaller "
         "still.")
report.p()
report.p("Every mechanism a polynomial-time algorithm could use is on the table "
         "of round 14. Smoothness-based ones are sub-exponential and not "
         "polynomial, because smooth numbers of size `N^c` are too rare for any "
         "fixed `c > 0`. The others, measured here and in twenty rounds before, "
         "sit at `N^(1/4)` or worse. A polynomial-time method along smoothness "
         "lines would need relations carried by numbers of size `N^(o(1))`, and "
         "no known construction produces them.")
report.write()
