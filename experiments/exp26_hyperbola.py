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
report.p("Binary search on `X` would cost `log N` pairs of walks. It is not needed. "
         "The hull of `{xy > N - 1}` contains every lattice point of `{xy > N}` "
         "plus the divisor points `(N/d, d)`, and those lie on the strictly convex "
         "curve `xy = N`, so each is an extreme point -- a *vertex* of the hull. "
         "One walk of that hull, testing `x * y = N` at each point, finds `p`.")
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
report.p("The walk stops at the divisor, so its cost is the number of hull "
         "vertices between rows `p` and `sqrt N` -- a constant times `N^(1/3)` per "
         "factor of two in `sqrt(N)/p`. Balanced semiprimes (`q < 2p` above) are its best "
         "case; the worst is `p` near `N^(1/3)`:")
report.p()
rows = []
for bits, pbits in ((44, 22), (44, 18), (44, 16), (50, 25), (50, 20), (50, 18)):
    st_all = []
    for _ in range(3):
        pp = randprime(1 << (pbits - 1), 1 << pbits)
        qq = randprime((1 << (bits - 1)) // pp + 1, (1 << bits) // pp)
        nn = pp * qq
        small = min(pp, qq)
        st = {}
        assert hyperbola_factor(nn, st) == small
        st_all.append((st["steps"] / nn ** (1 / 3), math.log2(math.sqrt(nn) / small)))
    rows.append([bits, pbits, f"{statistics.mean(a for a, _ in st_all):.2f}",
                 f"{statistics.mean(b for _, b in st_all):.1f}"])
report.table(["bits of N", "bits of p (target)", "steps / N^(1/3)", "log2(sqrt(N)/p)"], rows)
per_doubling = statistics.mean(float(r[2]) / float(r[3]) for r in rows if float(r[3]) > 2)
report.p(f"Away from `sqrt N` the cost is {per_doubling:.1f} `N^(1/3)` steps per "
         f"factor of two in `sqrt(N)/p`, as the curvature count predicts.")
report.p()
report.p("Exactly and deterministically, then, at `O~(N^(1/3))` in the worst "
         "case -- Lehman's exponent, and worse than the `N^(1/4)` of Strassen and "
         "rho, which pull away in the last column of the table before.")
report.p()

report.p("### The divisor vertex and Lehman")
report.p()
from fractions import Fraction

from aksfactor.hyperbola import divisor_vertex_edges

total = vertex = bracket = edges = in_box = 0
for _ in range(60):
    pp = randprime(1 << 15, 1 << 16)
    qq = randprime(pp + 2, 3 * pp)
    nn = pp * qq
    total += 1
    got = divisor_vertex_edges(nn)
    if not got or got[0] != pp:
        continue
    vertex += 1
    slopes = sorted(Fraction(dy, dx) for dx, dy, _, _ in got[2])
    bracket += len(slopes) == 2 and slopes[0] < Fraction(pp, qq) < slopes[1]
    for dx, dy, k, gap in got[2]:
        edges += 1
        c = dy * qq + dx * pp
        assert c * c - 4 * k * nn == gap * gap
        in_box += k <= nn ** (1 / 3) and c - 2 * math.sqrt(k * nn) <= nn ** (1 / 6) / (4 * math.sqrt(k))
report.table(["check", "count"], [
    ["divisor point is a hull vertex", f"{vertex} of {total}"],
    ["its two edges bracket the slope p/q", f"{bracket} of {vertex}"],
    ["edge (dx, dy) satisfies (dy q + dx p)^2 - 4(dx dy)N = (dy q - dx p)^2", f"{edges} of {edges}"],
    ["... with (k, c) = (dx dy, dy q + dx p) inside Lehman's search box", f"{in_box} of {edges}"],
])
report.p("The edges at the divisor are rational approximations of `p/q` from both "
         "sides, and each is literally a Fermat-Lehman certificate. But the hull "
         "and Lehman enumerate differently: Lehman walks `k <= N^(1/3)` and a "
         "short `c`-window for each; the hull walks the curve's own best "
         "approximations at the local curvature, and its certificates fall "
         "outside Lehman's box "
         f"{edges - in_box} times in {edges}. Same objects -- Farey fractions near "
         "`p/q`, the fan round 10 found in Lehman -- and the same exponent.")
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
