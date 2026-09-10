"""Round 9: an anatomy of Harvey's open question (Remark 3.4).

Harvey's paper closes its Remark 3.4 with an explicit open problem:

    "An interesting question is whether it is possible to obtain a fully
     square-root speedup for Lehman's original choice r ~ N^(1/3). This would
     presumably lead to a factoring algorithm with complexity N^(1/6+o(1))."

This experiment does not resolve it.  It locates it: measures which term of the
cost actually binds at the optimum, and measures whether the structure a
baby-step/giant-step sweep would need is present in the search space at all.
"""

from math import isqrt, log, sqrt

from _common import Report

from aksfactor.harvey import (
    candidate_counts,
    divisor_summatory,
    exponent_runs,
    run_length_bound,
)

report = Report(
    "exp17_remark34",
    "Where the N^(1/6) question binds",
    "Round 9. Not a resolution of Harvey's Remark 3.4 -- a measurement of it.",
)

report.p("## The cost, and its two terms")
report.p()
report.p("Harvey's candidate count is")
report.p()
report.p("```")
report.p("s  ~  sqrt(N) lg r / (2 sqrt(r) m)   +   r lg r")
report.p("      \\_______ j sweep ________/       \\_ one per (a,b) pair _/")
report.p("```")
report.p()
report.p("and the running time is `O(max(s, m, r))` up to log factors. The second "
         "term counts pairs with `ab <= r`, one candidate each at `j = 0`.")
report.p()
rows = []
for r in (10**2, 10**3, 10**4, 10**5, 10**6):
    exact = divisor_summatory(r)
    approx = r * (log(r) + 2 * 0.5772156649 - 1)
    rows.append([f"{r:,}", f"{exact:,}", f"{approx:,.0f}", f"{exact / approx:.4f}"])
report.table(["r", "exact #{(a,b) : ab <= r}", "r(ln r + 2γ − 1)", "ratio"], rows)
report.p("So the pair term is `Theta(r log r)` exactly, as claimed.")
report.p()

report.p("## Which term binds?")
report.p()
rows = []
for bits in (60, 100, 140, 180, 220, 260):
    n = 1 << bits
    r = m = round(2.0 ** (bits * 0.2))
    c = candidate_counts(n, r, m)
    rows.append([bits, f"{c['j_candidates']:.3e}", f"{c['pair_candidates']:.3e}",
                 f"{c['pair_candidates'] / c['j_candidates']:.4f}",
                 f"{c['pair_share']:.1%}"])
report.table(["bits of N", "j-candidates", "pair-candidates", "ratio", "pair share"],
             rows)
report.p("**The ratio is exactly 2, at every size.** At the optimum "
         "`r = m = N^(1/5)` the two terms are `N^(1/5) ln r / 2` and "
         "`N^(1/5) ln r`, so the pair enumeration is *always* two thirds of the "
         "candidate count. It is not an artifact of a particular `N`; it is what "
         "pins the exponent.")
report.p()
report.p("The arithmetic: `s >= r` forces `cost >= r`, while balancing the `j` "
         "term against the table size gives `m ~ N^(1/4)/r^(1/4)`. Setting "
         "`N^(1/4)/r^(1/4) = r` yields `r = N^(1/5)`. Reaching `N^(1/6)` with "
         "Lehman's `r = N^(1/3)` therefore requires sweeping the `Theta(r)` pairs "
         "themselves in `O(sqrt(r))` -- exactly the 'fully square-root speedup' "
         "Remark 3.4 asks for.")
report.p()

report.p("## Is there structure to sweep?")
report.p()
report.p("A baby-step/giant-step sweep needs the candidates to form a geometric "
         "progression, i.e. the exponents to form an arithmetic one. The exponents "
         "are")
report.p()
report.p("```")
report.p("e(a,b) = a*N + b - ceil(2 sqrt(abN))")
report.p("```")
report.p()
report.p("and the obstruction is visible: `sqrt(ab)` is not additive in `(a, b)`. "
         "But it could still be *locally* additive. Taking the `a = 1` slice, "
         "`e(1,b) = N + b - ceil(2 sqrt(bN))`, its first difference stays constant "
         "only while `f'(b) = sqrt(N/b)` moves by less than one, and "
         "`f''(b) = -sqrt(N)/(2 b^1.5)`, so a run has length at most about")
report.p()
report.p("```")
report.p("L(b)  <=  2 b^1.5 / sqrt(N).")
report.p("```")
report.p()
rows = []
for bits in (30, 36, 40, 44, 50, 60):
    n = 1 << bits
    r = min(max(4, round(n ** (1 / 3))), 20000)
    e = exponent_runs(n, r)
    rows.append([bits, f"{r:,}", e["longest"], f"{e['mean']:.2f}",
                 f"{e['predicted_bound']:.2f}", f"{sqrt(r):.0f}"])
report.table(["bits of N", "b up to r", "longest run measured", "mean run",
              "bound 2b^1.5/sqrt(N)", "run length BSGS needs (sqrt r)"], rows)
report.p("The measured runs track the bound, and both are `O(1)`. The sweep would "
         "need runs of length `sqrt(r) = N^(1/6)`; a run of length `L` requires "
         "`b ~ (L sqrt(N))^(2/3)`, so `L = N^(1/6)` needs `b ~ N^(4/9)`. But "
         "`b <= r = N^(1/3)`, and `1/3 < 4/9`. **The runs are never long enough, "
         "by a fixed margin in the exponent.**")
report.p()

report.p("## What this says about Remark 3.4")
report.p()
report.p("Not a resolution. A localisation, in three parts:")
report.p()
report.p("1. The binding constraint is the pair enumeration, and it binds by a "
         "factor of exactly 2 over the `j` sweep at every size. Any `N^(1/6)` "
         "algorithm must attack *that* term.")
report.p("2. The obvious attack -- BSGS over `(a,b)` -- needs local arithmetic "
         "progressions in `e(a,b)`, and the second-difference bound "
         "`L(b) <= 2 b^1.5/sqrt(N)` rules them out throughout Lehman's range. This "
         "is not a shortage of cleverness; the structure is measurably absent.")
report.p("3. So a square-root speedup over the pairs, if one exists, cannot come "
         "from the group structure of the candidates. It would have to come from "
         "the *arithmetic* of the pairs themselves -- the set `{(a,b) : ab <= r}` "
         "with `a/b` a convergent to `p/q`, which is a Farey/Stern-Brocot object, "
         "not a geometric progression.")
report.p()
report.p("That last observation is where this project's own thread ends up "
         "pointing. Lehman's candidates are good rational approximations to "
         "`p/q`; the search space is a Farey fan, not an interval. Whether a fan "
         "admits a square-root sweep is a question about continued fractions, not "
         "about binomial coefficients -- and nothing in Pascal's triangle speaks "
         "to it.")
report.write()
