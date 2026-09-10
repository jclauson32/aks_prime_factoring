"""Round 10: the Farey fan.  Can Lehman's candidate set be thinned?

Round 9 localised Harvey's Remark 3.4 to a single question: the `Theta(r log r)`
pairs `(a,b)` with `ab <= r` are two thirds of the cost, so an `N^(1/6)`
algorithm must handle them in `O(sqrt(r))`.  There are exactly two ways to do
that -- use fewer pairs, or process them faster.  Round 9 closed the second
(no local arithmetic progressions to sweep).  This round closes the first.

Each pair certifies `p` only inside an interval.  Lehman's condition
`0 <= aN/p + bp - 2 sqrt(abN) < W` becomes, after multiplying by `p`, the
quadratic `b p^2 - (2 sqrt(abN) + W) p + aN < 0`, so the certified `p` lie
between its roots -- centred at `p* = sqrt(aN/b)` where AM-GM is tight, with
half-width about `sqrt(N)/(2 b sqrt(r))`, depending only on `b`.

So the question becomes a covering problem: how few of these intervals suffice
to cover `[sqrt(N/r), sqrt(N))`?
"""

from math import sqrt

from _common import Report

from aksfactor.arith import is_prime
from aksfactor.harvey import lehman_intervals, minimum_subcover

report = Report(
    "exp18_farey",
    "The Farey fan: Lehman's covering is not redundant",
    "Round 10. The second of the two routes to N^(1/6), measured and closed.",
)

report.p("## Lehman's theorem, verified computationally")
report.p()
report.p("Lemma 3.3 says every `p` in `[sqrt(N/r), sqrt(N))` is certified by some "
         "pair. So the intervals must cover the range -- a claim worth checking "
         "rather than assuming.")
report.p()
rows = []
for bits, r in [(28, 60), (32, 100), (36, 150), (40, 220), (44, 300), (48, 420)]:
    n = float(1 << bits)
    iv, lo, hi = lehman_intervals(n, r)
    cover = minimum_subcover(iv, lo, hi)
    total_width = sum(h - l for l, h, _, _ in iv)
    rows.append([bits, r, f"{len(iv):,}", "yes" if cover else "NO",
                 f"{total_width / (hi - lo):.2f}"])
report.table(["bits of N", "r", "certified intervals", "covers the range?",
              "sum of widths / range"], rows)
report.p("It covers, at every size. And the total width is only about **1.8x** the "
         "range -- the system is barely thicker than it has to be.")
report.p()

report.p("## The covering is essentially non-redundant")
report.p()
report.p("If a small subfamily still covered the range, Harvey's algorithm could "
         "use only those pairs. Greedy interval covering is optimal, so this is "
         "the exact minimum:")
report.p()
rows = []
for r in (100, 200, 400, 800, 1600, 3200, 6400):
    bits = 2 * max(20, r.bit_length() * 3)
    n = float(1 << bits)
    iv, lo, hi = lehman_intervals(n, r)
    cover = minimum_subcover(iv, lo, hi)
    rows.append([f"{r:,}", f"{len(iv):,}", f"{cover:,}", f"{cover / len(iv):.3f}",
                 f"{cover / r:.2f}", f"{cover / sqrt(r):.1f}"])
report.table(["r", "pairs", "minimum subcover", "subcover / pairs",
              "subcover / r", "subcover / sqrt(r)"], rows)
report.p("**The ratio is flat at about 0.70** across a 64-fold range of `r`. At "
         "most 30% of the pairs are droppable, and the minimum subcover is "
         "`Theta(r log r)` -- the same order as the full set. Meanwhile "
         "`subcover / sqrt(r)` grows from 18 to 245, so the subcover is "
         "emphatically not `O(sqrt(r))`.")
report.p()
report.p("The reason is visible in the geometry. The half-width "
         "`sqrt(N)/(2 b sqrt(r))` depends only on `b`, while the centres "
         "`p* = sqrt(aN/b)` are spaced by about `sqrt(N)/(2 sqrt(ab))`. Width "
         "beats spacing only when `4a >= b r`, which at `ab = r` needs `a >= r/2` "
         "-- a vanishing corner of the fan. Everywhere else the intervals sit "
         "essentially edge to edge, and dropping one opens a gap.")
report.p()

report.p("## Both routes to N^(1/6) are now closed")
report.p()
report.p("| route | what it needs | status |")
report.p("|---|---|---|")
report.p("| **use fewer pairs** | a subfamily of size `O(sqrt r)` covering the range | closed here: minimum subcover is `0.70 x pairs = Theta(r log r)` |")
report.p("| **process them faster** | local arithmetic progressions in `e(a,b)` to sweep by BSGS | closed in round 9: runs are `O(1)`, need `sqrt(r)` |")
report.p()
report.p("Neither is a shortage of ingenuity. Both are measurable absences of "
         "structure: the covering is tight, and the exponent sequence is nowhere "
         "locally geometric.")
report.p()

report.p("## What is left")
report.p()
report.p("A `N^(1/6)` algorithm, if one exists, cannot come from thinning "
         "Lehman's fan or from sweeping it as a group. It would need a different "
         "*parametrisation* of the same search -- some coordinate system in which "
         "the good rational approximations to `p/q` are not `Theta(r)` scattered "
         "points but a structured object.")
report.p()
report.p("The natural candidate is the Stern-Brocot tree, where the convergents "
         "to any fixed `xi = p/q` form a path of length `O(log)` rather than a "
         "set of size `Theta(r)`. Descending that path needs one comparison per "
         "step: *is `a/b < p/q`?* And that comparison is")
report.p()
report.p("```")
report.p("a/b < p/q   <=>   aq < bp   <=>   aN < b p^2   <=>   p > sqrt(aN/b),")
report.p("```")
report.p()
report.p("i.e. exactly a threshold query *is `p` bigger than this?* -- which is "
         "Strassen's problem, costing `O~(sqrt(threshold))`. For thresholds in "
         "Lehman's range that is `O~(N^(1/4))` per comparison, so the tree walk "
         "costs more than the answer it is looking for.")
report.p()
report.p("That is a satisfying place for this thread to end. The Farey structure "
         "genuinely does compress the search -- from `Theta(r)` candidates to "
         "`O(log)` tree steps -- but each step is priced at exactly the barrier "
         "the compression was meant to avoid. The fan and the threshold oracle are "
         "the same problem wearing different coordinates.")
report.write()
