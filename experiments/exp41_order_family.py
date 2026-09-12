"""Round 35: the order mechanism, priced per draw.

Every order method buys the same thing: a group attached to `p` whose order it
hopes is smooth.  They differ in how many groups they can draw and what each
draw costs.  This round measures both, and multiplies them out into the expected
work for a split.
"""

import random
import statistics
import time
from math import log

from _common import Report

from aksfactor.arith import is_prime, sieve
from aksfactor.cyclo import pollard_pminus1, williams_pplus1
from aksfactor.divseq import curve_eds_mod, eds_at
from aksfactor.ecm import ecm_curve, is_b1b2_smooth
from aksfactor.smooth import curve_order, smooth_rate

report = Report(
    "exp41_order_family",
    "The order mechanism, priced per draw",
    "Round 35. How many groups each method can draw, and what each draw costs.",
)

rng = random.Random(41)
B1, B2 = 1000, 50000


def randprime(bits):
    while True:
        v = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if is_prime(v):
            return v


report.p("## How many draws, and how often each one wins")
report.p()
pbits = 20
primes = [randprime(pbits) for _ in range(200)]
pm1 = smooth_rate([p - 1 for p in primes], B1)
pp1 = smooth_rate([p + 1 for p in primes], B1)
orders = []
for p in primes[:80]:
    while True:
        a, x, y = rng.randrange(p), rng.randrange(p), rng.randrange(1, p)
        b = (y * y - x ** 3 - a * x) % p
        if (4 * a ** 3 + 27 * b * b) % p:
            break
    orders.append(curve_order(a, b, p))
ecm_draw = smooth_rate(orders, B1)
report.table(["family", "groups available per prime", f"chance one draw is {B1}-smooth"], [
    ["Pascal row, AKS fold, q-deformation", "1 (order divides `p^d - 1`)", f"{pm1:.2f} (as `p - 1`)"],
    ["Pollard p - 1", "1", f"{pm1:.2f}"],
    ["Williams p + 1", "2 (the two twists)", f"{pp1:.2f}"],
    ["class groups (Schnorr-Lenstra)", "one per multiplier `k`", "varies with `h(-kN)`"],
    ["ECM / the elliptic triangle", "one per curve, about `4 sqrt p` distinct orders",
     f"{ecm_draw:.2f} per curve"],
])
report.p(f"The rigid families get one ticket per prime: if `p - 1` is not smooth "
         f"to the bound, no amount of extra work helps. The elliptic families draw "
         f"again for the price of a curve, and each draw is about as likely to win "
         f"as `p - 1` was.")
report.p()

report.p("## What a draw costs")
report.p()
n = randprime(24) * randprime(40)
rows = []


def cost_of(fn, repeats=3):
    t0 = time.perf_counter()
    for _ in range(repeats):
        fn()
    return (time.perf_counter() - t0) / repeats


exp_bits = 1
for ell in sieve(B1):
    k = ell
    while k * ell <= B1:
        k *= ell
    exp_bits += k.bit_length()
rows.append(["Pollard p - 1", f"{cost_of(lambda: pollard_pminus1(n, B1)) * 1000:.1f} ms",
             "one exponentiation"])
rows.append(["Williams p + 1", f"{cost_of(lambda: williams_pplus1(n, B1)) * 1000:.1f} ms",
             "one Lucas chain"])
rows.append(["ECM, one curve (two stages)",
             f"{cost_of(lambda: ecm_curve(n, rng.randrange(6, n - 1), B1, B2)) * 1000:.1f} ms",
             "a Montgomery ladder and a baby-step giant-step stage 2"])
w = curve_eds_mod(0, 17, 2, 5, n)
m_index = 1
for ell in sieve(B1):
    k = ell
    while k * ell <= B1:
        k *= ell
    m_index *= k
rows.append(["elliptic triangle, row lcm(1..B1)",
             f"{cost_of(lambda: eds_at(*w, m_index, n)) * 1000:.1f} ms",
             "the sequence's own double-and-add"])
report.table(["one draw", "time", "what it is"], rows)
report.p("The elliptic draws cost more than the rigid ones -- a curve is a bigger "
         "object than a residue -- but they are draws, and the rigid families have "
         "only the one.")
report.p()

report.p("## Multiplying it out")
report.p()
report.p("Expected work for a split is the cost of a draw divided by the chance it "
         "wins. For the rigid families that quotient is meaningless when the single "
         "draw loses: the method simply fails, whatever the budget. For ECM it is a "
         "real number, and minimising it over `B1` is the calculation that produces "
         "`L_p[1/2]` -- larger `B1` makes each draw cost more and win more often, "
         "and the optimum balances them (round 31 measured the density that sets "
         "the balance).")
report.p()
report.p("Shor's algorithm belongs in this table too, in the row that does not "
         "exist classically: it uses one group, `(Z/N)*`, draws once, and needs no "
         "smoothness at all, because it reads the order instead of guessing a "
         "multiple of it (round 24). Every classical entry pays for not being able "
         "to do that.")
report.write()
