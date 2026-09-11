"""Round 18: factoring by counting lattice points under a hyperbola.

"Is there any way to do some sort of binary search?"  The predicate

    D(X) = #{d | N : d <= X}  =  sum_{y<=X} floor(N/y) - sum_{y<=X} floor((N-1)/y)

is monotone, and each sum is an exact lattice-point count computable without
visiting every row: walk the convex hull of the points above xy = N.
"""

import math
import random
import statistics
import time
from math import isqrt

from _common import Report

from aksfactor.arith import is_prime, pollard_rho, sieve
from aksfactor.hyperbola import (
    class_number,
    divisors_up_to,
    hyperbola_factor,
    hyperbola_piece_exponent,
    quadratic_floor_sum,
    row_sum,
    row_sum_naive,
)

report = Report(
    "exp26_hyperbola",
    "Factoring by counting lattice points under a hyperbola",
    "Round 18. A binary-searchable predicate with no factorial in it, and where it stops.",
)

rng = random.Random(18)


def randprime(lo, hi):
    while True:
        x = rng.randrange(lo, hi) | 1
        if is_prime(x):
            return x


report.p("## 1. The count, exactly")
report.p()
report.p("The lattice points above `xy = N` form a convex set, so none lies between "
         "its hull and the curve, and the number of points under the curve in row "
         "`y` is read off the hull edge crossing that row. A Stern-Brocot stack "
         "finds the next hull vertex in amortised `O(1)` steps.")
report.p()
checked = 0
for trial in range(400):
    n = rng.randrange(2, 10 ** 10)
    top = rng.randrange(1, isqrt(n) + 1)
    assert row_sum(n, top)[0] == row_sum_naive(n, top)
    checked += 1
rows = []
for e in (8, 10, 12, 14, 16):
    n = rng.randrange(10 ** e, 2 * 10 ** e)
    st = {}
    t0 = time.perf_counter()
    row_sum(n, isqrt(n), st)
    dt = time.perf_counter() - t0
    rows.append([f"~1e{e}", st["vertices"], st["steps"], f"{st['steps'] / n ** (1 / 3):.2f}",
                 f"{st['steps'] / (n ** (1 / 3) * math.log(n)):.3f}", f"{isqrt(n):,}",
                 f"{dt:.3f} s"])
report.p(f"Against an independent `O(sqrt N)` reference: **{checked} of {checked}** "
         f"random `(N, X)` pairs agree.")
report.p()
report.table(["N", "hull vertices", "steps", "steps / N^(1/3)",
              "steps / (N^(1/3) ln N)", "rows covered", "time"], rows)
drift = float(rows[-1][4]) / float(rows[0][4]) - 1
report.p(f"The step count grows like `N^(1/3) log N`: the fifth column moves by "
         f"{drift:+.0%} across eight orders of magnitude. Approximating the "
         f"hyperbola by segments of rational slope is also how Voronoi (1903) got "
         f"his `x^(1/3)` error term in the divisor problem; here the segments are "
         f"exact rather than approximate.")
report.p()

report.p("## 2. The predicate, and one pass instead of a binary search")
report.p()
p, q = 1009, 2003
n = p * q
rows = [[x, divisors_up_to(n, x)] for x in (p - 2, p - 1, p, p + 1, isqrt(n - 1))]
report.table(["X", f"#{{d | {n} : d <= X}}"], rows)
report.p("Binary search on `X` would cost `log N` pairs of walks. It is not needed: "
         "walking the hulls for `N` and `N - 1` side by side, their cumulative row "
         "counts agree exactly until the row of the first divisor, so one pair of "
         "walks locates `p`.")
report.p()

rows = []
for bits in (32, 38, 44, 50):
    steps, times, rho_times, found = [], [], [], 0
    trials = 4
    for _ in range(trials):
        pp = randprime(1 << (bits // 2 - 1), 1 << (bits // 2))
        qq = randprime(pp + 2, 2 * pp)
        nn = pp * qq
        st = {}
        t0 = time.perf_counter()
        d = hyperbola_factor(nn, st)
        times.append(time.perf_counter() - t0)
        found += d == pp
        steps.append(st["steps"] / nn ** (1 / 3))
        t0 = time.perf_counter()
        g = pollard_rho(nn, random.Random(1))
        rho_times.append(time.perf_counter() - t0)
        assert g in (pp, qq)
    rows.append([bits, f"{found} of {trials}", f"{statistics.mean(steps):.1f}",
                 f"{statistics.mean(times):.3f} s", f"{statistics.mean(rho_times):.4f} s"])
report.table(["bits of N", "p found", "steps / N^(1/3)", "hyperbola walk",
              "Pollard rho"], rows)
report.p("It works, exactly and deterministically, at `O~(N^(1/3))` -- Lehman's "
         "exponent, and worse than the `N^(1/4)` of Strassen and rho, which pull "
         "away in the last column. It is a size method: the hull's edge "
         "directions are Farey fractions, the same fan round 10 found in Lehman.")
report.p()

report.p("## 3. Where it stops: curved pieces and class numbers")
report.p()
report.p("Each hull edge is an exact *linear* piece of `floor(N/y)`, summed in "
         "`O(1)`. If exact degree-`d` pieces could be summed as cheaply, the piece "
         "count per dyadic block would fall from `N^(1/3)` to `N^(1/(d+2))`: a "
         "degree-`d` Taylor piece of length `h` at `x` misses `N/y` by "
         "`N h^(d+1) / x^(d+2)`, and `O(1)` stray lattice points per piece allows "
         "`h ~ x N^(-1/(d+2))`.")
report.p()
known = {1: "Lehman; this walk", 2: "Strassen, Pollard rho", 3: "Harvey (deterministic)",
         4: "nothing known -- below every deterministic method", 5: "nothing known"}
rows = [[d, f"N^(1/{d + 2})", f"{hyperbola_piece_exponent(d):.3f}", known[d]]
        for d in range(1, 6)]
report.table(["piece degree d", "pieces", "exponent", "same exponent as"], rows)
primes = [x for x in sieve(2000) if x > 3 and x % 4 == 3]
ok = sum(quadratic_floor_sum(x) == (x - 1) * (2 * x - 1) // 6 - (x - 1 - 2 * class_number(-x)) // 2
         for x in primes)
report.p(f"The obstruction arrives at `d = 2`. By Dirichlet's class number formula, "
         f"the full-period quadratic floor sum is a class number in disguise:")
report.p()
report.p("```")
report.p("sum_{k<p} floor(k^2/p)  =  (p-1)(2p-1)/6  -  (p-1-2 h(-p))/2      (p = 3 mod 4)")
report.p("```")
report.p()
report.p(f"(checked for **{ok} of {len(primes)}** primes below 2000). No "
         f"polynomial-time algorithm for `h(-p)` is known, so exact quadratic floor "
         f"sums in polylogarithmic time would be news in their own right. The pieces "
         f"the hyperbola needs are short arcs, not full periods, so this is a "
         f"warning rather than a proof -- but it puts the next step in a known "
         f"hard neighbourhood.")
report.p()
report.p("A second requirement hides in \"`O(1)` stray points per piece\": the "
         "linear walk never has strays, because hull edges are chosen so that no "
         "lattice point lies between edge and curve. A curved analogue needs the "
         "lattice points *near* a short parabolic arc located exactly -- small "
         "fractional parts of a quadratic sequence, a Diophantine problem with no "
         "known Euclid-like algorithm.")
report.write()
