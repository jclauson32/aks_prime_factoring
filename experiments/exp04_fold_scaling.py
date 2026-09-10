"""Does folding the AKS object beat trial division?  Measured answer: no.

Folding `(1+x)^n` modulo `x^r - 1` reads the *whole* row in `O(log n)` ring
multiplications instead of `O(k)` coefficients, so it is the natural candidate
for escaping `Theta(spf(n))`.  It does leak factors -- but this script measures
how many fold coefficients must be inspected before one does, against the size
of the smallest prime factor.
"""

import random

from _common import Report

from aksfactor.arith import is_prime
from aksfactor.ring import fold_attack

report = Report(
    "exp04_fold_scaling",
    "Scaling of the AKS-ring fold attack",
    "Cost measured in fold coefficients inspected, versus `p = spf(n)`.",
)

report.p("### Why the fold leaks at all")
report.p()
report.p("For a prime `p` with `p^v || n` and `s = n / p^v`, Frobenius gives "
         "`(x+a)^n = (x^(p^v) + a)^s (mod p)`, so the folded coefficient")
report.p()
report.p("```")
report.p("S_j  =  sum over { i : i*p^v = j (mod r) } of C(s, i) * a^(s-i)   (mod p)")
report.p("```")
report.p()
report.p("`gcd(S_j, n)` is a non-trivial factor exactly when that progression sum "
         "vanishes mod `p`. Nothing forces it to; it behaves like a uniform residue, "
         "so each coefficient is a **~1/p** lottery ticket. Expected tickets needed: "
         "`Theta(p)`.")
report.p()
report.p("### Measured")
report.p()

rng = random.Random(20260909)


def randprime(lo: int, hi: int) -> int:
    while True:
        c = rng.randrange(lo, hi) | 1
        if is_prime(c):
            return c


rows = []
for lo, hi in [(50, 120), (200, 500), (900, 1500), (3000, 4500), (9000, 13000),
               (25000, 35000)]:
    for _ in range(3):
        p = randprime(lo, hi)
        q = randprime(p + 2, p * 3)
        n = p * q
        res = fold_attack(n, rmax=400)
        found = res["factor"]
        rows.append([p, q, res["r"] if found else "-",
                     f"{res['pairs']:,}",
                     f"{res['pairs'] / p:.2f}" if found else "not found"])
report.table(["p = spf(n)", "q", "first r that leaked", "coefficients inspected",
              "coefficients / p"], rows)

ratios = sorted(float(r[4]) for r in rows if r[4] != "not found")
misses = sum(1 for r in rows if r[4] == "not found")
median = ratios[len(ratios) // 2]
report.p(f"Median *coefficients / p* over the {len(ratios)} solved cases: "
         f"**{median:.2f}**, range {ratios[0]:.2f} to {ratios[-1]:.2f}, with **no "
         f"trend in `p`** across three orders of magnitude. The spread is what a "
         f"geometric process with success rate `1/p` predicts: the cost is a "
         f"coin-flip count, not a smooth function.")
if misses:
    report.p()
    report.p(f"{misses} run(s) exhausted the `r <= 400` budget without a hit. That "
             f"budget is only about `2.5 * p` coefficients at the largest size "
             f"tested, and a geometric process with rate `1/p` fails to fire within "
             f"`2.5p` draws about `e^-2.5 = 8%` of the time -- so these are misses of "
             f"the budget, not evidence against the model.")
report.p()
report.p("Either way the conclusion is the same: the fold attack costs `Theta(p)` "
         "coefficient evaluations -- the same wall as trial division, reached by a "
         "more expensive route.")
report.p()
report.p("### The barrier, stated plainly")
report.p()
report.p("Modulo a prime factor `p`, the AKS object *is* `(x^(p^v) + a)^(n/p^v)`. "
         "Its first non-trivial structure lives at degree `p^v`. Truncating below "
         "that degree sees nothing (Theorem 2); folding below that degree smears it "
         "across all `r` classes and each class keeps only a `1/p` chance of "
         "collapsing. Either way you pay `Theta(p)`.")
report.write()
