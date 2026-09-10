"""Round 12: the search has been reduced to an equivalent problem, not an easier one.

Round 11 showed factoring is polynomial time *given* `p` to within `N^(1/4)`, and
that guessing costs `N^(1/4)` by counting.  The obvious hope is that partial
information is cheaper to obtain than `p` itself.  This round shows it is not,
in two steps:

1. Coppersmith's window is `N^(beta^2)` where `p ~ N^beta`.  Since `p` is the
   smaller factor, `beta <= 1/2`, so the window never exceeds `N^(1/4)` -- and
   the balanced semiprime is the case where it is widest.
2. Guessing over that window costs `N^(beta(1-beta))`, while Strassen finds `p`
   outright in `N^(beta/2)`.  The first is never smaller.  So the lattice route
   never beats the algebraic one, and they coincide at exactly `N^(1/4)`.

Which gives the honest end state: for balanced semiprimes, polynomial-time
factoring and polynomial-time approximation of `p` are *equivalent*.
"""

import random
from math import log2

from _common import Report

from aksfactor.arith import is_prime
from aksfactor.lattice import cost_exponents, factor_with_hint

report = Report(
    "exp20_equivalence",
    "Partial information is exactly as hard as the factorisation",
    "Round 12. The reduction closes into an equivalence.",
)

rng = random.Random(4)


def randprime(bits):
    while True:
        x = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if is_prime(x):
            return x


report.p("## 1. The window is `N^(beta^2)`, and `beta <= 1/2`")
report.p()
report.p("Coppersmith recovers a root smaller than `N^(beta^2/d)`, and for "
         "factoring the polynomial is linear (`d = 1`). Writing `p ~ N^beta`:")
report.p()
rows = []
for pb, qb in [(14, 14), (12, 18), (10, 22), (9, 26), (8, 30)]:
    p = randprime(pb)
    q = randprime(qb)
    n = p * q
    nb = n.bit_length()
    beta = log2(p) / log2(n)
    predicted = beta * beta * nb
    best = 0
    for unknown in range(1, int(predicted) + 4):
        if factor_with_hint(n, p & ~((1 << unknown) - 1), bound=1 << unknown,
                            m=4, beta=beta):
            best = unknown
    rows.append([nb, f"{beta:.3f}", best, f"{predicted:.1f}",
                 f"{best / predicted:.2f}" if predicted else "-"])
report.table(["bits of N", "beta", "bits actually recovered",
              "beta^2 · log2 N", "ratio"], rows)
report.p("The window tracks `beta^2 log N`. And `p` is the *smaller* factor, so "
         "`beta <= 1/2` always: **the widest possible Coppersmith window is "
         "`N^(1/4)`, attained exactly at the balanced semiprime.** There is no "
         "regime where the lattice sees further.")
report.p()

report.p("## 2. Guessing the window never beats Strassen")
report.p()
report.p("Covering a range of length `N^beta` with windows of width `N^(beta^2)` "
         "takes `N^(beta - beta^2)` guesses. Strassen finds `p` outright in "
         "`O~(sqrt p) = O~(N^(beta/2))`.")
report.p()
rows = []
for beta in (0.50, 0.45, 0.40, 0.33, 0.25, 0.20, 0.10):
    c = cost_exponents(beta)
    if abs(c["guess_and_coppersmith"] - c["strassen"]) < 1e-9:
        winner = "tie"
    elif c["strassen"] < c["guess_and_coppersmith"]:
        winner = "Strassen"
    else:
        winner = "Coppersmith"
    rows.append([f"{beta:.2f}", f"{c['coppersmith_window']:.4f}",
                 f"{c['guess_and_coppersmith']:.4f}", f"{c['strassen']:.4f}",
                 winner])
report.table(["beta", "window exponent beta^2", "guess+Coppersmith beta(1-beta)",
              "Strassen beta/2", "winner"], rows)
report.p("`beta(1-beta) <= beta/2` iff `beta >= 1/2`, and `beta <= 1/2` always. So "
         "**guessing plus Coppersmith is never strictly better than Strassen**, "
         "and is strictly worse for every unbalanced semiprime. The two meet at "
         "`N^(1/4)` precisely at `beta = 1/2` -- the balanced case, which is also "
         "where `beta(1-beta)` is maximised. The hardest input for the lattice "
         "route is the hardest input for everything else.")
report.p()

report.p("## 3. The equivalence")
report.p()
report.p("Putting the two together gives a statement rather than another barrier.")
report.p()
report.p("> **Theorem.** For balanced semiprimes `N = pq` with `p <= q <= 2p`, the "
         "following are equivalent:")
report.p("> ")
report.p("> **(a)** a deterministic `poly(log N)` algorithm that factors `N`;")
report.p("> ")
report.p("> **(b)** a deterministic `poly(log N)` algorithm that outputs an integer "
         "`p~` with `|p - p~| <= N^(1/4)`.")
report.p()
report.p("*Proof.* (a) implies (b): factor and output `p`. (b) implies (a): feed "
         "`p~` to Coppersmith, whose window at `beta = 1/2` is exactly `N^(1/4)`. "
         "∎")
report.p()
report.p("The implication that matters is the second one, and it is the reason "
         "this round exists. Eleven rounds hunted for *partial* information about "
         "`p` on the theory that a quarter of the bits might be cheaper than all "
         "of them. They are the same problem. Any algorithm producing a quarter "
         "of the bits in polynomial time **is** a polynomial-time factoring "
         "algorithm, and every barrier measured in this repository applies to it "
         "unchanged.")
report.p()

report.p("## What this settles, and what it does not")
report.p()
report.p("It settles the strategy. Looking for cheap partial information is not a "
         "way around the barrier; it is the barrier restated. The search space of "
         "*approaches* has collapsed to two doors, both long known and both shut:")
report.p()
report.p("- **widen the window**: improve Coppersmith's exponent past `beta^2/d`. "
         "A well-known barrier in lattice cryptanalysis; the bound is tight for "
         "this lattice family, so it needs a genuinely different construction.")
report.p("- **shrink the search**: find structure in the *location* of `p` that no "
         "algebraic, geometric, group-theoretic or fractal property examined over "
         "twelve rounds provided.")
report.p()
report.p("It does not settle the question. Factoring is not known to be hard -- it "
         "is in NP and co-NP, so it is very unlikely to be NP-complete, and a "
         "polynomial-time algorithm is not ruled out by anything proved. What this "
         "project can report is a map: eleven distinct routes, each closed by a "
         "measurement rather than an intuition, and a precise statement of what "
         "any twelfth route must do that none of them did.")
report.write()
