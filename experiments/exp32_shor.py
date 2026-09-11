"""Round 24: Shor's algorithm, simulated -- the order mechanism without smoothness.

Every classical order method in this repository exponentiates blindly by a
number divisible by everything small and hopes the order r is smooth.  Shor's
algorithm finds r itself, by interference over Q >= N^2 values of a^x at once.
Simulating it classically is exact and easy; it is also exponential, and the
place the cost lands is the whole story.
"""

import random
import time

from _common import Report

from aksfactor.arith import factorize
from aksfactor.shor import order_from_outcome, shor, shor_distribution

report = Report(
    "exp32_shor",
    "Shor's algorithm, simulated classically",
    "Round 24. The order mechanism without the smoothness condition, and what it costs to fake.",
)

rng = random.Random(32)

report.p("## 1. The interference pattern")
report.p()
n, a = 21, 2
t, probs, c = shor_distribution(n, a, random.Random(3))
q = 1 << t
r = next(k for k in range(1, n) if pow(a, k, n) == 1)
peaks = sorted(range(q), key=lambda y: -probs[y])[:r]
report.p(f"`N = {n}`, `a = {a}`, order `r = {r}`, first register `t = {t}` qubits "
         f"(`Q = {q}`). After measuring `a^x mod N = {c}`, the Fourier transform puts "
         f"the first register's probability on multiples of `Q/r = {q / r:.1f}`:")
report.p()
report.table(["outcome y", "probability", "y r / Q", "order read by continued fractions"],
             [[y, f"{probs[y]:.3f}", f"{y * r / q:.3f}", order_from_outcome(y, t, n, a)]
              for y in sorted(peaks)])
report.p(f"The {r} peaks hold {sum(probs[y] for y in peaks):.0%} of the probability. "
         f"Outcomes `y` with `gcd(y r / Q, r) = 1` give `r` directly; the others give "
         f"a divisor, recovered by trying small multiples.")
report.p()

report.p("## 2. Factoring with it")
report.p()
rows = []
series = []
for n in (15, 21, 33, 35, 51, 77, 91, 143, 221, 323, 437, 667, 899, 1147, 1517):
    st = {}
    t0 = time.perf_counter()
    g = shor(n, rng, stats=st, lucky_gcd=False)
    dt = time.perf_counter() - t0
    assert g and n % g == 0 and 1 < g < n
    per_run = dt / st["runs"]
    series.append((n, per_run, st["t"]))
    rows.append([n, f"{g} x {n // g}", st["t"], 1 << st["t"], st["a"], st["r"],
                 max(factorize(st["r"])), st["runs"], f"{per_run:.2f} s"])
report.table(["N", "factors", "qubits t", "amplitudes Q = 2^t", "base a", "order r",
              "largest prime of r", "runs", "time per run"], rows)
(n0, t0_, q0), (n1, t1_, q1) = series[-4], series[-1]
report.p("Every base here was coprime to `N`, so every factor came out of the "
         "period finding, not a lucky gcd. At these sizes the orders are small, "
         "hence smooth, and Pollard's `p - 1` would have found them too; the "
         "simulation cannot reach the sizes where the difference shows. What it "
         "does show is that nothing in the procedure asks whether `r` is smooth: "
         "the Fourier transform reads `r` whatever its factorisation.")
report.p()

report.p("## 3. Where the cost went")
report.p()
report.p(f"The simulation's time per run follows `Q = 2^t ~ N^2`: from `N = {n0}` to "
         f"`N = {n1}` the register grew {2 ** (q1 - q0)}x, from `2^{q0}` to `2^{q1}` "
         f"amplitudes, and the time per run {t1_ / t0_:.1f}x. That is `N^2` -- worse "
         f"than trial division's `N^(1/2)` -- because a classical machine has to write "
         f"down every amplitude.")
report.p()
report.table(["resource", "quantum computer", "this simulation"], [
    ["qubits / memory", "`t ~ 2 log2 N` qubits", "`2^t ~ N^2` complex amplitudes"],
    ["modular exponentiation", "`O(t)` multiplications, applied to all `x` at once",
     "`2^t` evaluations of `a^x mod N`"],
    ["Fourier transform", "`O(t^2)` gates", "an FFT of length `2^t`"],
    ["total", "polynomial in `log N`", "`Theta(N^2 log N)`"],
])
report.p("That is the exact shape of the gap. The order mechanism needs `r`; "
         "classically `r` can only be used blind, as an exponent divisible by it, "
         "which is why smoothness decides everything. Interference reads `r` off a "
         "superposition of all `Q` exponents. Nothing in twenty-three classical "
         "rounds found a way to read it without writing the `Q` values down, and "
         "by Miller's reduction, reading it any other way would factor `N`.")
report.write()
