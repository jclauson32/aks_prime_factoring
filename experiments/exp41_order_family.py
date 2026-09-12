"""Round 35: the order mechanism, priced per draw.

Every order method buys the same thing: a group attached to `p` whose order it
hopes is smooth.  They differ in how many groups they can draw, how likely a
draw is to win, and what each draw costs.  This round measures all three and
multiplies them out -- and then checks the model by running the methods.

Everything below uses one bound `B1` and one prime size, so the numbers compose.
"""

import random
import statistics
import time
from math import isqrt, sqrt

from _common import Report

from aksfactor.arith import is_prime, multiplicative_order, sieve
from aksfactor.cyclo import pollard_pminus1, williams_pplus1
from aksfactor.divseq import curve_eds_mod, eds_at
from aksfactor.ecm import ecm_curve, is_b1b2_smooth, pminus1
from aksfactor.smooth import qr_table, suyama_order, weierstrass_order

report = Report(
    "exp41_order_family",
    "The order mechanism, priced per draw",
    "Round 35. How many groups each method can draw, how often a draw wins, "
    "what a draw costs -- and whether the product predicts the run.",
)

rng = random.Random(41)
B1 = 200                      # one bound throughout: stage 1 only, no stage 2
PBITS, QBITS = 15, 45         # p small enough to count points on curves over F_p


def randprime(bits):
    while True:
        v = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if is_prime(v):
            return v


def powersmooth(m):
    """Exactly what stage 1 kills: every prime power dividing `m` is `<= B1`."""
    return is_b1b2_smooth(m, B1, B1)


def rate_se(hits, total):
    r = hits / total
    return r, sqrt(r * (1 - r) / total)


# ---------------------------------------------------------------- section 1
report.p("## How many draws each family gets")
report.p()
report.table(["family", "groups available per prime", "if the draw loses"], [
    ["Pascal row, AKS fold, q-deformation", "1 (the order divides `p^d - 1`)",
     "nothing to do"],
    ["Pollard p - 1", "1", "nothing to do"],
    ["Williams p + 1", "2 (the base's discriminant decides whether the "
     "order divides `p - 1` or `p + 1`)", "nothing to do"],
    ["class groups (Schnorr-Lenstra)", "one per multiplier `k`", "change `k`"],
    ["ECM / the elliptic triangle", f"one per curve, filling the Hasse interval "
     f"of about `4 sqrt p` orders", "change the curve"],
])
report.p("The split in the last column is the whole distinction between the "
         "rigid families and the elliptic ones. It is not about the group being "
         "elliptic -- class groups are not -- but about whether a fresh group is "
         "available at a price.")
report.p()

# ---------------------------------------------------------------- section 2
report.p("## How often one draw wins")
report.p()
report.p(f"All four samples are orders attached to primes of {PBITS} bits, and "
         f"the test is the one stage 1 actually applies: every prime power at "
         f"most `B1 = {B1}`. The last row is the control -- a random integer "
         f"from the Hasse interval, which is what a curve order would be if it "
         f"carried no bias of its own. It does carry one, so the control is what "
         f"measures it.")
report.p()

shift_sample = 3000
shift_primes = [randprime(PBITS) for _ in range(shift_sample)]
pm1_rate, pm1_se = rate_se(sum(powersmooth(p - 1) for p in shift_primes), shift_sample)
pp1_rate, pp1_se = rate_se(sum(powersmooth(p + 1) for p in shift_primes), shift_sample)

curve_primes = [randprime(PBITS) for _ in range(150)]
suyama_orders, uniform_orders, hasse_control, ambiguous = [], [], [], 0
for p in curve_primes:
    table = qr_table(p)
    lo, hi = p + 1 - 2 * isqrt(p), p + 2 + 2 * isqrt(p)
    for _ in range(8):
        order = suyama_order(rng.randrange(6, p), p, table)
        if order is None:
            ambiguous += 1
        else:
            suyama_orders.append(order)
        while True:
            a, b = rng.randrange(p), rng.randrange(p)
            if (4 * a ** 3 + 27 * b * b) % p:
                break
        uniform_orders.append(weierstrass_order(a, b, p, table))
        hasse_control.append(rng.randrange(lo, hi))

suy_rate, suy_se = rate_se(sum(map(powersmooth, suyama_orders)), len(suyama_orders))
uni_rate, uni_se = rate_se(sum(map(powersmooth, uniform_orders)), len(uniform_orders))
ctl_rate, ctl_se = rate_se(sum(map(powersmooth, hasse_control)), len(hasse_control))

measured = [
    ("`p - 1` (Pollard)", pm1_rate, pm1_se, shift_sample),
    ("`p + 1` (Williams)", pp1_rate, pp1_se, shift_sample),
    ("`#E` for a Suyama curve (ECM)", suy_rate, suy_se, len(suyama_orders)),
    ("`#E` for a uniform curve", uni_rate, uni_se, len(uniform_orders)),
    ("random integer in the Hasse interval (control)", ctl_rate, ctl_se, len(hasse_control)),
]
report.table(["what has to be smooth", f"chance it is {B1}-powersmooth", "samples"],
             [[name, f"{r:.3f} +- {se:.3f}", n] for name, r, se, n in measured])

ranking = sorted(measured, key=lambda row: -row[1])
report.p("Ranked by the measured rate: " +
         ", ".join(f"{name} ({r:.3f})" for name, r, _, _ in ranking) + ".")


def compare(a, b):
    """A sentence fragment stating the measured relation between two rows."""
    (na, ra, sa, _), (nb, rb, sb, _) = a, b
    sigma = sqrt(sa * sa + sb * sb)
    gap = (ra - rb) / sigma if sigma else 0.0
    if abs(gap) < 2:
        return f"{na} and {nb} are level ({ra:.3f} against {rb:.3f}, {abs(gap):.1f} sigma)"
    word = "beats" if ra > rb else "loses to"
    return (f"{na} {word} {nb} ({ra:.3f} against {rb:.3f}, "
            f"{abs(gap):.1f} sigma)")


report.p()
report.p("The comparisons that matter, each computed from the samples above:")
report.p()
report.p(f"- {compare(measured[0], measured[2])}. The two draws are not "
         f"interchangeable, and the direction is worth pinning down: the curve "
         f"ECM actually runs wins more often than `p - 1` does, because Suyama's "
         f"parametrisation hands it a factor of 12 for free (next section).")
report.p(f"- {compare(measured[3], measured[4])}. A uniform curve order is *not* "
         f"quite a random integer of its size -- it carries small factors more "
         f"often, the bias Galbraith and McKee describe. Round 31 measures the "
         f"same gap against a plain random integer at higher power "
         f"(+5.4 sigma there).")
suy_vs_uni = abs(suy_rate - uni_rate) / sqrt(suy_se ** 2 + uni_se ** 2)
report.p(f"- {compare(measured[2], measured[3])}. "
         + ("Suyama's parametrisation is not neutral -- see the next section."
            if suy_vs_uni >= 2 else
            "At this sample size the parametrisation's bias does not separate "
            "from a uniform curve; the next section measures the bias directly "
            "instead of through its effect on smoothness."))
report.p()

# ---------------------------------------------------------------- section 3
report.p("## What Suyama's parametrisation is worth")
report.p()
div12 = sum(o % 12 == 0 for o in suyama_orders)
div12_uniform = sum(o % 12 == 0 for o in uniform_orders)
report.table(["curves", "orders resolved", "divisible by 12"],
             [["Suyama (what ECM runs)", len(suyama_orders),
               f"{div12}/{len(suyama_orders)}"],
              ["uniform `y^2 = x^3 + ax + b`", len(uniform_orders),
               f"{div12_uniform}/{len(uniform_orders)}"]])
report.p(f"The order was counted over `F_p` and the twist resolved by the "
         f"starting point ({ambiguous} draws where the point's order divided "
         f"both were discarded). Suyama's curves carry a rational 12-torsion "
         f"structure, so a factor of 12 is free on every draw and the cofactor "
         f"that still has to be smooth is smaller by that much, against "
         f"{div12_uniform / len(uniform_orders):.0%} of uniform curves getting "
         f"the same discount. "
         + (f"On the smoothness rate that discount was worth "
            f"{suy_rate - uni_rate:+.3f} here ({suy_vs_uni:.1f} sigma)."
            if suy_vs_uni >= 2 else
            f"On the smoothness rate the discount is worth "
            f"{suy_rate - uni_rate:+.3f} here, which this sample "
            f"({suy_vs_uni:.1f} sigma) cannot separate from zero -- the bias is "
            f"a constant factor, and a constant is what it looks like.")
         + " Either way it is bought once and never improved on: no "
           "parametrisation makes the remaining cofactor smooth more than a "
           "constant more often.")
report.p()

# ---------------------------------------------------------------- section 4
report.p("## What a draw costs")
report.p()
n_cost = randprime(PBITS) * randprime(QBITS)
w = curve_eds_mod(0, 17, 2, 5, n_cost)
m_index = 1
for ell in sieve(B1):
    k = ell
    while k * ell <= B1:
        k *= ell
    m_index *= k


def cost_of(fn, repeats=20):
    t0 = time.perf_counter()
    for _ in range(repeats):
        fn()
    return (time.perf_counter() - t0) / repeats


costs = [
    ["Pollard p - 1", cost_of(lambda: pollard_pminus1(n_cost, B1)), "one exponentiation"],
    ["Williams p + 1", cost_of(lambda: williams_pplus1(n_cost, B1)), "one Lucas chain"],
    ["ECM, one curve", cost_of(lambda: ecm_curve(n_cost, rng.randrange(6, n_cost - 1), B1)),
     "a Montgomery ladder"],
    ["elliptic triangle, row lcm(1..B1)", cost_of(lambda: eds_at(*w, m_index, n_cost)),
     "the sequence's own double-and-add"],
]
report.table(["one draw", "time", "what it is"],
             [[name, f"{t * 1000:.2f} ms", what] for name, t, what in costs])
ecm_cost = dict((row[0], row[1]) for row in costs)["ECM, one curve"]
report.p("The elliptic draw costs more than the rigid ones -- a curve is a "
         "bigger object than a residue -- but it is a draw, and the rigid "
         "families have only the one.")
report.p()

# ---------------------------------------------------------------- section 5
report.p("## Does the model predict the run?")
report.p()
report.p("The table above says a curve wins with probability "
         f"{suy_rate:.3f} and `p - 1` with probability {pm1_rate:.3f}. Both are "
         "now run on fresh semiprimes, one draw each.")
report.p()

trials = 400
ecm_wins = pm1_wins = pm1_agree = 0
won_without = lost_with = order_smooth = both_smooth = 0
for _ in range(trials):
    p, q = randprime(PBITS), randprime(QBITS)
    n = p * q
    if ecm_curve(n, rng.randrange(6, n - 1), B1) is not None:
        ecm_wins += 1
    got = pminus1(n, B1) is not None
    predicted = powersmooth(p - 1)
    pm1_wins += got
    pm1_agree += (got == predicted)
    if got and not predicted:             # won anyway: is the base's order smooth?
        won_without += 1
        order_smooth += powersmooth(multiplicative_order(2, p))
    if predicted and not got:             # lost anyway: did the gcd come out n?
        lost_with += 1
        both_smooth += powersmooth(q - 1)

ecm_run, ecm_run_se = rate_se(ecm_wins, trials)
pm1_run, pm1_run_se = rate_se(pm1_wins, trials)
report.table(["method, one draw", "predicted from the order sample", "measured on "
              f"{trials} semiprimes", "agree?"],
             [["ECM, one Suyama curve", f"{suy_rate:.3f} +- {suy_se:.3f}",
               f"{ecm_run:.3f} +- {ecm_run_se:.3f}",
               "yes" if abs(suy_rate - ecm_run) < 2 * sqrt(suy_se ** 2 + ecm_run_se ** 2)
               else "no"],
              ["Pollard p - 1", f"{pm1_rate:.3f} +- {pm1_se:.3f}",
               f"{pm1_run:.3f} +- {pm1_run_se:.3f}",
               "yes" if abs(pm1_rate - pm1_run) < 2 * sqrt(pm1_se ** 2 + pm1_run_se ** 2)
               else "no"]])
report.p(f"Per instance, `p - 1` succeeded exactly when `p - 1` was "
         f"{B1}-powersmooth on {pm1_agree} of {trials} numbers: the method is "
         f"barely probabilistic, it is close to a lookup on a property of `p` "
         f"that was decided before the method started. ECM's draw is the "
         f"probabilistic one.")
report.p()
if trials - pm1_agree:
    parts = [f"The {trials - pm1_agree} disagreements were chased down rather "
             f"than waved at."]
    if won_without:
        parts.append(f"{won_without} were wins without `p - 1` being "
                     f"powersmooth, and in {order_smooth} of those the order of "
                     f"the base `2` modulo `p` was powersmooth even though "
                     f"`p - 1` was not -- what stage 1 needs is the order of the "
                     f"base, and a divisor can be smoother than the number it "
                     f"divides.")
    if lost_with:
        parts.append(f"{lost_with} were losses with `p - 1` powersmooth, and in "
                     f"{both_smooth} of those `q - 1` was powersmooth too, so "
                     f"the gcd came out `n` and the split was lost to the second "
                     f"prime being as weak as the first.")
    report.p(" ".join(parts))
else:
    report.p("There were no disagreements at all: on these 400 numbers the "
             "method won exactly when `p - 1` was powersmooth.")
report.p()

# ---------------------------------------------------------------- section 6
report.p("## Multiplying it out")
report.p()
expected = ecm_cost / suy_rate if suy_rate else float("inf")
report.table(["family", "chance a draw wins", "cost per draw", "expected cost of a split"],
             [["Pollard p - 1", f"{pm1_rate:.3f}",
               f"{dict((r[0], r[1]) for r in costs)['Pollard p - 1'] * 1000:.2f} ms",
               f"succeeds on {pm1_rate:.0%} of primes; on the rest, never"],
              ["ECM (redrawn per curve)", f"{suy_rate:.3f}", f"{ecm_cost * 1000:.2f} ms",
               f"{1 / suy_rate:.1f} curves, {expected * 1000:.1f} ms"]])
report.p("That quotient is the whole of ECM's analysis. Raising `B1` raises the "
         "cost of a draw roughly linearly and raises the chance it wins by the "
         "Dickman factor (round 31); minimising the quotient over `B1` is what "
         "produces `L_p[1/2]`, and it is a minimisation over a function whose "
         "smooth-number term falls like `u^(-u)`. Nothing in the order mechanism "
         "is sub-exponential except through that term.")
report.p()
report.p("Shor's algorithm belongs in this table too, in the row that does not "
         "exist classically: one group, one draw, no smoothness at all, because "
         "it reads the order instead of guessing a multiple of it (round 24). "
         "Every classical entry here is paying for not being able to do that.")
report.write()
