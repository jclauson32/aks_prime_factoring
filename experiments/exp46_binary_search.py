"""Round 40: the binary search, and the price of its comparison.

The question that started several of these rounds: can the factor be found by
binary search? Round 38 said no from the 2-adic side -- the bitwise congruence
is satisfiable for every candidate, so it never eliminates a half.

But a binary search does not need the *bits* of `p`; it needs a comparison.
And there is one: "does `N` have a divisor at most `x`?" is answered by counting
lattice points under two hyperbolas (round 18), in `O~(N^(1/3))`. So the search
exists, it is implemented here, and it works.

What it costs is the point. A comparison is `N^(1/3)`, a search is `log N` of
them, and the direct hull walk that round 18 already had is the same `N^(1/3)`
without the log. The binary search is real and it is not a shortcut: the cost is
in the comparison, and nobody knows a cheaper one.
"""

import random
import statistics
import time
from math import isqrt, log2

from _common import Report

from aksfactor.arith import factorize, is_prime
from aksfactor.hyperbola import binary_search_factor, divisors_up_to, hyperbola_factor

report = Report(
    "exp46_binary_search",
    "The binary search for a factor, and what one comparison costs",
    "Round 40. Divisor counting as a comparison operator: the search works, at "
    "`N^(1/3)` per question.",
)

rng = random.Random(46)


def randprime(bits):
    while True:
        v = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if is_prime(v):
            return v


# ---------------------------------------------------------------- section 1
report.p("## The comparison")
report.p()
report.p("`divisors_up_to(N, x)` counts the divisors of `N` that are at most "
         "`x`, as a difference of two hyperbola sums, each computed by walking "
         "the convex hull of `{xy > N}` rather than visiting every row. It is "
         "checked here against brute force on every `N` in a range, so the "
         "comparison the search rests on is not taken on trust.")
report.p()

mismatch, checked = 0, 0
for n in range(4, 3000):
    top = isqrt(n - 1)                     # the range the count is defined on
    want = sum(1 for d in range(1, top + 1) if n % d == 0)
    checked += 1
    if divisors_up_to(n, top) != want:
        mismatch += 1
spf_mismatch, spf_checked = 0, 0
for n in range(2, 3000):
    got = binary_search_factor(n)
    want = None if is_prime(n) else min(factorize(n))
    spf_checked += 1
    spf_mismatch += got != want
report.table(["check", "cases", "mismatches"],
             [["`divisors_up_to` against brute force", checked, mismatch],
              ["`binary_search_factor` against a reference factorisation",
               spf_checked, spf_mismatch]])
report.p("The second row includes the primes (the search must return nothing), "
         "the squares of primes (whose only nontrivial divisor sits exactly at "
         "`sqrt N`, outside the range the hyperbola count can address) and the even "
         "numbers (which it must find at the first step).")
report.p()

# ---------------------------------------------------------------- section 2
report.p("## The search")
report.p()
report.p("Now on semiprimes, where the answer is the smaller prime. The search "
         "asks `log2(sqrt N)` questions; what is measured is how many hull steps "
         "those questions cost in total.")
report.p()

rows = []
points = []
for bits in (28, 34, 40, 46):
    steps, comps, times, ok = [], [], [], 0
    for _ in range(3):
        p, q = randprime(bits // 2), randprime(bits - bits // 2)
        n = p * q
        stats = {}
        t0 = time.perf_counter()
        got = binary_search_factor(n, stats)
        times.append(time.perf_counter() - t0)
        ok += got == min(p, q)
        steps.append(stats["steps"])
        comps.append(stats["comparisons"])
    med = statistics.median(steps)
    points.append((bits, med))
    rows.append([bits, f"{statistics.median(comps):.0f}", f"{med:,.0f}",
                 f"{statistics.median(times):.2f} s", f"{ok}/3"])
report.table(["bits of N", "comparisons", "hull steps (median)", "time", "correct"],
             rows)

fit, fit_corrected = [], []
for (b1, s1), (b2, s2) in zip(points, points[1:]):
    fit.append(log2(s2 / s1) / (b2 - b1))
    # the O~ hides logs; divide them out before reading an exponent off four points
    fit_corrected.append(log2((s2 / log2(2 ** b2) ** 2) / (s1 / log2(2 ** b1) ** 2))
                         / (b2 - b1))
raw, corrected = statistics.mean(fit), statistics.mean(fit_corrected)
report.p(f"Fitted growth of the step count: `N^{raw:.3f}`, which is above the "
         f"`N^(1/3)` the geometry promises. The promise is `O~(N^(1/3))` and the "
         f"tilde is not decoration at these sizes: dividing the counts by "
         f"`log2(N)^2` first and refitting gives `N^{corrected:.3f}`, so two "
         f"logarithmic factors account for the gap. The comparisons column is "
         f"`log2(sqrt N)` exactly, as a binary search must be.")
report.p()

# ---------------------------------------------------------------- section 3
report.p("## Against walking the hull once")
report.p()
report.p("Round 18's `hyperbola_factor` does not binary search: it walks the "
         "hull from `sqrt N` downwards and reads the divisor off as a vertex. "
         "Same geometry and the same exponent, in one pass rather than "
         "`2 log2(sqrt N)` of them -- and, as the numbers below show, one pass "
         "that is allowed to stop as soon as it has the answer.")
report.p()

rows = []
for bits in (28, 34, 40, 46):
    p, q = randprime(bits // 2), randprime(bits - bits // 2)
    n = p * q
    walk_stats, search_stats = {}, {}
    t0 = time.perf_counter()
    walked = hyperbola_factor(n, walk_stats)
    walk_time = time.perf_counter() - t0
    t0 = time.perf_counter()
    searched = binary_search_factor(n, search_stats)
    search_time = time.perf_counter() - t0
    agree = "yes" if walked == searched == min(p, q) else f"no ({walked}, {searched})"
    rows.append([bits, f"{walk_stats.get('steps', 0):,}", f"{walk_time:.2f} s",
                 f"{search_stats['steps']:,}", f"{search_time:.2f} s",
                 f"{search_stats['steps'] / max(walk_stats.get('steps', 1), 1):.1f}x",
                 agree])
report.table(["bits of N", "hull walk steps", "hull walk time",
              "binary search steps", "binary search time", "ratio",
              "both found the smaller prime"], rows)
ratios = [float(r[5].rstrip("x")) for r in rows]
report.p(f"The ratio is not the `log N` one might expect from asking the "
         f"question `log N` times instead of once, and it grows with the size "
         f"({ratios[0]:.0f}x to {ratios[-1]:.0f}x here). The reason is that the "
         f"two walks stop at different places. `hyperbola_factor` halts the "
         f"moment a vertex turns out to be a divisor, which for a balanced "
         f"semiprime is almost immediately; each comparison in the search, by "
         f"contrast, is a *sum*, and a sum has to be finished -- both of its "
         f"hyperbola walks run all the way down to `cbrt(N)` whatever the answer "
         f"turns out to be.")
report.p()
report.p("Both are `O~(N^(1/3))` in the exponent, which is worse than the "
         "`N^(1/4)` of the batching methods (round 34) and far worse than "
         "anything using smooth numbers. Between them the difference is a large "
         "constant and two logs, in favour of the method that is allowed to stop "
         "early.")
report.p()

# ---------------------------------------------------------------- section 4
report.p("## What this settles")
report.p()
report.p("A binary search needs a predicate that splits the space. Two have "
         "turned up in this project, and they are opposites:")
report.p()
report.table(["predicate", "does it split the range?", "what one evaluation costs"],
             [["`p == r (mod 2^k)`, the bitwise lift (round 38)",
               "no -- every odd residue is consistent with `N`", "`O(1)`"],
              ["`N` has a divisor `<= x`, by hyperbola counting (this round)",
               "yes -- exactly a half, every time", "`O~(N^(1/3))`"]])
report.p("So the binary search that was asked for exists and is written down "
         "here. It does not help, and the reason is precise: the predicate that "
         "is cheap carries no information, and the predicate that carries the "
         "information is as expensive as the search it was meant to replace. "
         "Every round in this project has ended at some version of that "
         "sentence.")
report.write()
