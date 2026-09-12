"""Round 30: which mechanism pays first, as a function of budget.

Round 27 measured time to factor. This measures the other axis: given a fixed
budget of modular multiplications, how often does each mechanism split `N` at
all?  That curve is what selection climbed in round 28, and it says why the
mechanisms overtake one another in the order they do.

Costs are counted in modular multiplications, with the usual estimates:
exponentiation `1.5 log2 E`, a Lucas chain `2 log2 E`, a Montgomery ladder
step 11 per bit of the stage-1 exponent, rho 3 per step (map, accumulator,
amortised gcd), trial division 1 per candidate.
"""

import math
import random
import statistics

from _common import Report

from aksfactor.arith import is_prime, sieve
from aksfactor.collision import pascal_rho
from aksfactor.cyclo import williams_pplus1
from aksfactor.ecm import ecm
from aksfactor.evolve import pminus1_success

report = Report(
    "exp37_budget",
    "Which mechanism pays first",
    "Round 30. Success probability against a budget of modular multiplications.",
)

rng = random.Random(37)
PBITS, QBITS = 24, 40
SAMPLES = 60
BUDGETS = (100, 1000, 10000, 100000)


def randprime(bits):
    while True:
        v = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if is_prime(v):
            return v


cases = []
while len(cases) < SAMPLES:
    p, q = randprime(PBITS), randprime(QBITS)
    n = p * q
    cases.append((n, p, q, rng.randrange(2, n - 1)))


def lcm_upto(b):
    e = 1
    for ell in sieve(b):
        k = ell
        while k * ell <= b:
            k *= ell
        e *= k
    return e


def best_smooth_exponent(budget, cost_per_bit):
    """The largest `lcm(1..B)` whose exponentiation fits the budget."""
    bits = budget / cost_per_bit
    b = 2
    best = 1
    while True:
        e = lcm_upto(b)
        if e.bit_length() > bits:
            return best, b - 1
        best, b = e, b + 1


def pminus1_rate(budget):
    e, _ = best_smooth_exponent(budget, 1.5)
    return pminus1_success(e, [(n, a) for n, _, _, a in cases])


def pplus1_rate(budget):
    _, b1 = best_smooth_exponent(budget, 2.0)
    hits = 0
    for n, p, q, _ in cases:
        g = williams_pplus1(n, max(b1, 2))
        hits += g in (p, q)
    return hits / len(cases)


def rho_rate(budget):
    steps = int(budget / 3)
    hits = 0
    for n, p, q, a in cases:
        g, _ = pascal_rho(n, 2, 2 + (a % 50), max_steps=steps)
        hits += g in (p, q)
    return hits / len(cases)


def ecm_rate(budget):
    best = 0.0
    for b1 in (50, 200, 1000, 5000, 20000):
        per_curve = 11 * lcm_upto(b1).bit_length() + 30 * (50 * b1) / 210
        curves = int(budget / per_curve)
        if curves < 1:
            continue
        hits = 0
        for n, p, q, _ in cases:
            g = ecm(n, b1=b1, b2=50 * b1, curves=curves, rng=random.Random(n))
            hits += g in (p, q)
        best = max(best, hits / len(cases))
    return best


def trial_rate(budget):
    return sum(p <= budget for _, p, _, _ in cases) / len(cases)


methods = [("trial division (size)", trial_rate),
           ("Pollard p - 1 (order, rigid)", pminus1_rate),
           ("Williams p + 1 (order, rigid)", pplus1_rate),
           ("ECM, best B1 for the budget (order, redrawn)", ecm_rate),
           ("Pollard rho (collision)", rho_rate)]

rows = []
table = {}
for name, fn in methods:
    row = [name]
    for b in BUDGETS:
        r = fn(b)
        table[name, b] = r
        row.append(f"{r:.2f}")
    rows.append(row)
report.p(f"{SAMPLES} semiprimes with a {PBITS}-bit `p` and a {QBITS}-bit `q`; each "
         f"method is given the budget and tuned to it (the smooth exponent or the "
         f"ECM bound that fits). Entries are the fraction of the {SAMPLES} numbers "
         f"split.")
report.p()
report.table(["method"] + [f"{b:,} mults" for b in BUDGETS], rows)

winners = []
for b in BUDGETS:
    winners.append((b, max(methods, key=lambda m: table[m[0], b])[0]))
report.p("Best at each budget: " + "; ".join(f"{b:,} -- {w}" for b, w in winners) + ".")
report.p()
rho_needs = math.sqrt(2 ** PBITS) * 3
report.p(f"The order of overtaking is set by how each mechanism's probability grows "
         f"with effort. `p - 1` pays immediately -- its first few small primes are "
         f"the ones most likely to divide `p - 1` -- and then saturates at the "
         f"fraction of primes whose `p - 1` is smooth enough to reach. Rho pays "
         f"nothing until its walk approaches `sqrt p` steps (about "
         f"{rho_needs:,.0f} multiplications here) and then pays everything. ECM "
         f"spends its budget on many cheap curves and climbs steadily, because "
         f"every curve is a fresh draw. Trial division is linear in the budget and "
         f"needs `p` of them.")
report.p()
report.p("Two methods are missing on purpose. The sieves spend their budget on "
         "additions in a sieve array, not on modular multiplications, so putting "
         "them on this axis would need a different unit; and Strassen's method "
         "covers about `(budget/c)^2` candidates for a constant `c` that depends on "
         "the polynomial arithmetic, which is an implementation number rather than "
         "a property of the mechanism.")
report.p()
fine = [3000, 6000, 12000, 24000, 48000]
rho_fine = [(b, rho_rate(b)) for b in fine]
pred = [(b, 1 - math.exp(-((b / 3) ** 2) / (2 * 2 ** PBITS))) for b in fine]
report.p("Rho against the birthday prediction `1 - exp(-steps^2 / 2p)`, with "
         "`steps = budget/3`:")
report.p()
report.table(["budget", "measured", "birthday prediction"],
             [[f"{b:,}", f"{m:.2f}", f"{q:.2f}"] for (b, m), (_, q) in zip(rho_fine, pred)])
report.p()
report.p("This is the curve selection climbed in round 28: at a budget of a few "
         "dozen multiplications, the only mechanism with a gradient is the order "
         "mechanism with a fixed group, so that is what evolution found -- and the "
         "reason it stopped there is visible in the same table, since the "
         "mechanisms that overtake it do so only after thousands of "
         "multiplications.")
report.write()
