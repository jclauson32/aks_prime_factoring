"""Round 11: the polynomial-time algorithm -- and the exact price of its input.

Ten rounds of searching produced barriers, never an algorithm faster than
`N^(1/5)`.  There is, however, a genuinely *polynomial-time* factoring
algorithm, and it has been sitting outside this project's scope the whole time:

    Coppersmith's method.  Given `N = pq` and an approximation to `p` accurate to
    within about `N^(1/4)`, it recovers `p` in time polynomial in `log N`.

No search, no smoothness, no group order.  It builds a lattice whose short
vectors are polynomials with the *same* small root over the integers as `f` has
modulo `p`, LLL-reduces, and reads the root off.

So factoring *is* polynomial time -- conditional on `N^(1/4)`-precision
information about `p`.  This round implements it, measures how much information
it really needs, and then shows why that condition cannot be met cheaply.
"""

import random
import time
from math import isqrt, log2, sqrt

from _common import Report

from aksfactor.arith import is_prime
from aksfactor.harvey import lehman_intervals
from aksfactor.lattice import factor_with_hint, hint_bits_needed, lll

report = Report(
    "exp19_coppersmith",
    "Coppersmith: polynomial time, at a price of exactly N^(1/4)",
    "Round 11. The algorithm exists. The information it needs is the whole problem.",
)

rng = random.Random(11)


def randprime(bits):
    while True:
        x = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if is_prime(x):
            return x


report.p("## It works, and it is polynomial time")
report.p()
rows = []
for pb in (13, 15, 17, 19):
    p = randprime(pb)
    q = randprime(pb)
    if p == q:
        continue
    n = p * q
    unknown = max(1, hint_bits_needed(n) - 3)
    t0 = time.time()
    got = factor_with_hint(n, p & ~((1 << unknown) - 1), bound=1 << unknown, m=4)
    dt = time.time() - t0
    rows.append([n.bit_length(), f"{p:,}", unknown, "yes" if got else "no",
                 f"{dt:.2f}"])
report.table(["bits of N", "p", "unknown low bits of p", "factored", "seconds"], rows)
report.p("Given the high bits of `p`, the factorisation falls out with no search "
         "at all.")
report.p()

report.p("## How much information does it need?")
report.p()
report.p("Coppersmith's bound is `X < N**(beta^2/d)` asymptotically in the lattice "
         "dimension; at finite dimension it falls short by about `1/m`. Measured "
         "live, sweeping the number of unknown bits until recovery "
         "fails:")
report.p()
sweep_p = randprime(20)
sweep_q = randprime(20)
while sweep_q == sweep_p:
    sweep_q = randprime(20)
sweep_n = sweep_p * sweep_q
sweep_need = hint_bits_needed(sweep_n)
rows = []
for m in (2, 3, 5):
    t0 = time.time()
    best = 0
    for unknown in range(1, sweep_need + 2):
        hint = sweep_p & ~((1 << unknown) - 1)
        if factor_with_hint(sweep_n, hint, bound=1 << unknown, m=m):
            best = unknown
    rows.append([m, 2 * m, best, sweep_need, f"{best / sweep_need:.2f}",
                 f"{time.time() - t0:.1f}"])
report.table(
    ["m", "lattice dimension", "max unknown bits recovered", "(1/4)·log2 N",
     "fraction of the limit", "seconds"],
    rows,
)
report.p("The achievable bound climbs toward `N^(1/4)` as the lattice grows, "
         "exactly as the theory says, and the cost of LLL climbs with it. The "
         "limit is `N^(1/4)`: **a quarter of the bits of `N`, or equivalently half "
         "the bits of `p`.**")
report.p()

report.p("## And why that cannot be bootstrapped")
report.p()
report.p("The obvious move is to guess. Cover the range of possible `p` with "
         "intervals of width `N^(1/4)` and run Coppersmith on each centre. That "
         "converts a polynomial-time conditional algorithm into an unconditional "
         "one -- at a cost fixed by pure counting:")
report.p()
report.p("```")
report.p("covering an interval of length L with intervals of width w")
report.p("needs at least L/w of them.")
report.p("```")
report.p()
rows = []
for bits in (64, 128, 256, 512, 1024, 2048):
    # exact, in exponent form -- 2**2048 overflows a float
    rows.append([bits, f"2^{bits // 2}", f"2^{bits // 4}", f"2^{bits // 4}"])
report.table(["bits of N", "range of p", "Coppersmith window", "guesses needed"],
             rows)
report.p("`sqrt(N) / N^(1/4) = N^(1/4)`, always. Guess-and-Coppersmith is "
         "`Theta(N^(1/4))` -- the same exponent Strassen reached in 1977, and it "
         "**cannot** reach `N^(1/5)`, let alone `N^(1/6)`. This is not a heuristic "
         "or a measured tendency; it is a counting bound with no escape.")
report.p()
report.p("Nor does any cleverer organisation of the guesses help, and it is worth "
         "being precise about why. The counting bound depends on **only two "
         "numbers**: the length of the range and the width of the window. No "
         "structure imposed on the guesses changes either. Lehman's fan, for "
         "instance, has windows that are wider or narrower than Coppersmith's "
         "depending on `r`:")
report.p()
rows = []
for bits, r in [(32, 100), (40, 220), (48, 420)]:
    n = float(1 << bits)
    iv, lo, hi = lehman_intervals(n, r)
    widths = sorted(h - l for l, h, _, _ in iv)
    median = widths[len(widths) // 2]
    rows.append([bits, r, f"{len(iv):,}", f"2^{log2(median):.1f}",
                 f"2^{bits / 4:.1f}", f"{median / (n**0.25):.4f}"])
report.table(["bits of N", "r", "pairs", "median interval width",
              "N^(1/4)", "ratio"], rows)
report.p("At these `r` the Lehman windows are *wider* than `N^(1/4)`, so each needs "
         "several Coppersmith calls; at `r ~ N^(1/3)` they are narrower (width "
         "`~N^(1/6)`), so one Coppersmith call spans many of them. Either way the "
         "total is unchanged: covering a range of length `sqrt(N)` with windows of "
         "width `N^(1/4)` takes `N^(1/4)` windows, however they are grouped. The "
         "fan reorganises the guesses; it does not reduce them.")
report.p()

report.p("## Where this leaves the search")
report.p()
report.p("The result is oddly clarifying. Factoring **is** polynomial time given "
         "`N^(1/4)`-precision knowledge of `p`; Coppersmith converts that "
         "information into a factorisation for free. So the entire difficulty of "
         "the problem is the *cost of obtaining a quarter of the bits*, and every "
         "route this project examined pays at least `N^(1/4)` for them:")
report.p()
report.p("| route | what it buys | cost |")
report.p("|---|---|---|")
report.p("| trial division / Pascal row scan | `p` exactly | `Theta(p)` |")
report.p("| Strassen / the gasket threshold | `p` exactly | `O~(N^(1/4))` |")
report.p("| Harvey's BSGS sweep | `p` exactly | `O~(N^(1/5))` |")
report.p("| guess + Coppersmith | `p` exactly | `Theta(N^(1/4))`, by counting |")
report.p("| Coppersmith alone | `p` exactly | `poly(log N)` -- **given `N^(1/4)` of it** |")
report.p()
report.p("Two ways out, both famous open problems. Improve Coppersmith's exponent "
         "`beta^2/d` -- widening the window past `N^(1/4)` would immediately beat "
         "everything above, and is a well-known barrier in lattice cryptanalysis. "
         "Or find a source of high-order bits of `p` costing less than `N^(1/4)`, "
         "which is what all eleven rounds here failed to do.")
report.p()
report.p("The honest summary of the search: there is a polynomial-time algorithm, "
         "its input is a quarter of the bits, and nothing classical is known that "
         "produces those bits for less than the cost of just finding `p` outright.")
report.write()
