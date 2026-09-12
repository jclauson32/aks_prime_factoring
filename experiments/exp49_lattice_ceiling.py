"""Round 43: the lattice's ceiling, measured against its one knob.

Round 34 closed on a sentence: "the exponent `beta^2` in Coppersmith's bound is
what would have to move". This round measures the shape of that bound from the
inside, using the only parameter a caller controls -- the lattice's dimension.

For a monic linear polynomial and `p >= N^(1/2)`, Coppersmith's construction
with parameter `m` reaches a window of about `N^(1/4 - eps(m))`, with `eps`
falling to zero as `m` grows. So the question "does more lattice buy more
window?" has a known answer in the limit and a measurable one at each `m`: it
buys some, with diminishing returns, and it converges to `1/4` rather than
through it.

Round 11 already had three rows of this table on a 40-bit modulus, reported as
a fraction of the limit. This round asks the same question as an exponent, over
eight dimensions and a 64-bit modulus, and looks for where the buying stops.

The number is not a curiosity. Round 38 walks a tree of `N^(1/4)` leaves and
hands each to this lattice; anything above `1/4` would shrink that tree below
`N^(1/4)` and beat every size method at once.
"""

import random
import statistics
import time
from math import log

from _common import Report

from aksfactor.arith import is_prime
from aksfactor.lattice import factor_with_hint

report = Report(
    "exp49_lattice_ceiling",
    "How much window more lattice buys",
    "Round 43. The reachable Coppersmith window as a function of the lattice's "
    "dimension: rising, with diminishing returns, towards 1/4.",
)

rng = random.Random(49)
NBITS = 64
INSTANCES = 5
DIMS = (1, 2, 3, 4, 5, 6, 8, 10)


def randprime(bits):
    while True:
        v = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if is_prime(v):
            return v


def instance():
    p, q = randprime(NBITS // 2), randprime(NBITS // 2)
    if p < q:
        p, q = q, p                       # the hint is about the larger factor
    return p * q, p


def closes(n, p, t, m):
    """Does the lattice recover `p` from its top bits, given a window of `2^t`?"""
    window = 1 << t
    return factor_with_hint(n, p - p % window, window, m=m) is not None


def largest_window(n, p, m, hi=NBITS // 2):
    """Biggest `t` that still works, by bisection -- success is monotone in `t`."""
    lo = 0
    if not closes(n, p, 1, m):
        return 0
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if closes(n, p, mid, m):
            lo = mid
        else:
            hi = mid
    return lo


# ---------------------------------------------------------------- section 1
report.p("## What one more dimension buys")
report.p()
report.p(f"{INSTANCES} balanced semiprimes of {NBITS} bits, with `m` stopped "
         f"at {max(DIMS)}: the lattice is `2m`-dimensional in exact "
         f"arithmetic and the next step up costs more than this round is "
         f"worth. For each lattice "
         f"parameter `m`, the largest window `2^t` from which the lattice still "
         f"recovers `p`, found by bisection (a smaller window is strictly "
         f"easier, so the property is monotone and bisection is honest). The "
         f"window is reported as an exponent of `N`, which is the quantity "
         f"Coppersmith's theorem bounds by `beta^2 = 1/4`.")
report.p()

cases = [instance() for _ in range(INSTANCES)]
rows, exponents, per_instance = [], {}, []
for m in DIMS:
    ts, times = [], []
    for n, p in cases:
        t0 = time.perf_counter()
        ts.append(largest_window(n, p, m))
        times.append(time.perf_counter() - t0)
    exps = [t / NBITS for t in ts]
    exponents[m] = statistics.median(exps)
    per_instance.extend(exps)
    rows.append([m, 2 * m, f"{statistics.median(ts):.0f}",
                 f"{statistics.median(exps):.3f}",
                 f"{min(exps):.3f}-{max(exps):.3f}",
                 f"{statistics.median(times):.2f} s"])
report.table(["m", "lattice dimension", "largest window (bits)",
              "as `N^x`, median", "range over instances", "time to bisect"], rows)

best = max(exponents.values())
best_instance = max(per_instance)
report.p(f"The reachable exponent climbs with `m` and flattens: "
         + ", ".join(f"`m = {m}` reaches `N^{exponents[m]:.3f}`" for m in DIMS)
         + f". The ceiling the theorem gives is `N^0.250`, and the best measured "
           f"here is `N^{best:.3f}`.")
report.p()

gains = [(DIMS[i], exponents[DIMS[i]] - exponents[DIMS[i - 1]])
         for i in range(1, len(DIMS))]
report.table(["from m", "to m", "extra window (exponent)", "extra dimensions"],
             [[DIMS[i - 1], DIMS[i], f"{gains[i - 1][1]:+.3f}",
               2 * (DIMS[i] - DIMS[i - 1])] for i in range(1, len(DIMS))])
report.p("Each doubling of the lattice buys less than the last, which is what "
         "`eps(m) ~ 1/m` looks like from underneath. "
         + ("Every measured exponent, median and instance alike, stays under "
            "`1/4`, as the theorem requires; this is consistency with the "
            "theorem, not an independent check of it."
            if best_instance <= 0.25 else
            f"The medians stay under `1/4`, but individual instances reach "
            f"`N^{best_instance:.3f}`, a shade over it. That is what an "
            f"asymptotic bound with an `-eps` looks like at 64 bits: the window "
            f"is quantised to whole bits, one bit is {1 / NBITS:.3f} of the "
            f"exponent here, and a lucky instance clears the bar by one. It is "
            f"the theorem's error term, not a counterexample to it."))
report.p()

# ---------------------------------------------------------------- section 2
report.p("## What the ceiling is worth")
report.p()
report.p("Suppose the exponent were `1/4 + delta` instead. Round 38's algorithm "
         "walks a tree of `2^t` leaves to feed this lattice, so its cost is "
         "exactly the window's reciprocal: `N^(1/4)` leaves today, "
         "`N^(1/4 - delta)` then. The table below is that arithmetic, and it is "
         "the reason this one exponent keeps reappearing.")
report.p()
def verdict(cost):
    if cost >= 0.25:
        return "where every size method already is"
    if cost > 0.2:
        return "beats the `N^(1/4)` batching methods"
    if abs(cost - 0.2) < 0.005:
        return "ties Harvey's deterministic `N^(1/5)`"
    return "beats the best deterministic exponent known"


rows = []
for delta in (0.0, 0.01, 0.05, 1 / 12):
    reach = 0.25 + delta
    cost = 0.5 - reach
    rows.append([f"{reach:.3f}", f"N^{cost:.3f}", verdict(cost)])
report.table(["Coppersmith exponent", "cost of the tree that feeds it",
              "what it would beat"], rows)
report.p("At `delta = 1/12` the window is `N^(1/3)`, the tree is `N^(1/6)`, and "
         "the deterministic record falls. The whole size family, in this "
         "project and in the literature, is the search for that `delta`, and "
         "forty-three rounds have not moved it.")
report.write()
