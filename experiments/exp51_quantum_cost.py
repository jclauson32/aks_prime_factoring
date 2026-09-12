"""Round 45: the other side of the wall, priced in qubits.

Round 41 priced the classical wall: 2048-bit `N` at about 10^15 core-years
along a curve calibrated on RSA-768 and RSA-250. This round asks what the
quantum side costs, because the one algorithm that is not stuck behind that
wall has been known since 1994 and the reason it is not used is not
mathematical.

Two kinds of number appear below and they must not be mixed up. The first
section is *measured*: what it costs this machine to simulate Shor's algorithm,
which is the only part of the quantum story this repository can actually run.
The second is *quoted*: published resource estimates for running it, which are
engineering forecasts and not measurements of anything here.

The gap between them is the whole point. Simulating the algorithm costs
`2^(2n)`; running it would cost a polynomial. That difference is not a detail
of the simulator -- it is what "quantum" means here, and round 25 measured its
consequence from the other side when every classically computable function of
`a^x mod N` came back with a flat spectrum.
"""

import random
import time
from math import log, log2

from _common import Report

from aksfactor.arith import factorize, is_prime
from aksfactor.shor import shor, shor_distribution

report = Report(
    "exp51_quantum_cost",
    "The other side of the wall, in qubits",
    "Round 45. What simulating Shor costs here (measured), what running it "
    "would cost (quoted), and why the two differ by an exponential.",
)

rng = random.Random(51)

# ---------------------------------------------------------------- section 1
report.p("## Measured: simulating the algorithm")
report.p()
report.p("The simulator holds an amplitude for every value of the first "
         "register, and that register needs `2 log2 N` qubits for the period to "
         "be readable -- so the state vector has about `N^2` entries. Below, "
         "the cost of one run, on composites small enough to fit.")
report.p()

rows = []
points = []
for bits in (5, 6, 7, 8, 9, 10):
    while True:
        n = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        f = factorize(n)
        if len(f) >= 2 and not is_prime(n):        # odd, composite, not a power
            break
    a = 2
    while True:
        a = rng.randrange(2, n)
        if all(a % p for p in f):
            break
    t0 = time.perf_counter()
    t, probs, _ = shor_distribution(n, a, rng)
    dt = time.perf_counter() - t0
    entries = 1 << t
    points.append((bits, dt))
    rows.append([n, bits, t, f"{entries:,}", f"{entries * 16 / 1e6:.1f} MB",
                 f"{dt:.2f} s"])
report.table(["N", "bits of N", "qubits in the first register",
              "amplitudes held", "at 16 bytes each", "one run"], rows)

growth = []
for (b1, t1), (b2, t2) in zip(points, points[1:]):
    if t1 > 0:
        growth.append(log2(t2 / t1) / (b2 - b1))
per_bit = sum(growth) / len(growth) if growth else float("nan")
report.p(f"Each extra bit of `N` multiplies the simulation's cost by "
         f"{2 ** per_bit:.1f} -- the state vector is `N^2` amplitudes, so the "
         f"factor to expect is 4 per bit, and the excess is the FFT's own"
         f" logarithm. No "
         f"extrapolation is needed to see where this ends: a 30-bit `N` wants "
         f"`2^60` amplitudes, {(1 << 60) * 16 / 1e18:.0f} exabytes at 16 bytes "
         f"each, which is why round 24 stops where it stops.")
report.p()

n_small = 21
t0 = time.perf_counter()
stats = {}
got = shor(n_small, random.Random(7), stats=stats, lucky_gcd=False)
elapsed = time.perf_counter() - t0
report.p(f"The full algorithm still runs, on numbers this size: `{n_small}` "
         + (f"factors as `{min(got, n_small // got)} x "
            f"{max(got, n_small // got)}` after "
            f"{stats.get('runs', 0)} period-finding runs in {elapsed:.2f} s"
            if got else
            f"was not split within {stats.get('runs', 0)} runs this time")
         + ", with the lucky-gcd shortcut disabled so every success goes "
           "through the interference rather than around it.")
report.p()

# ---------------------------------------------------------------- section 2
report.p("## Quoted: running the algorithm")
report.p()
report.p("Nothing in this section is measured here. It is what the published "
         "resource estimates say, and it is included because round 41's classical "
         "table is meaningless without it.")
report.p()
report.table(["quantity", "how it scales", "for `N` of 2048 bits"],
             [["logical qubits", "`O(n)` -- around `3n` in the modern layouts",
               "a few thousand"],
              ["Toffoli or `T` gates", "`O(n^3)` with schoolbook multiplication, "
               "`O(n^2 log n)` with fast", "about `10^10`"],
              ["physical qubits, surface code", "`O(n log n)` times the code "
               "distance squared", "about 20 million (Gidney-Ekerå, 2021)"],
              ["wall clock", "hours, at the same reference",
               "about 8 hours"]])
report.p("Set against round 41's classical figure for the same modulus -- about "
         "`10^15` core-years, or a billion years on a million cores -- the "
         "comparison is not close, and it is not a comparison between two "
         "algorithms on one machine. It is a comparison between a machine that "
         "exists and one that does not.")
report.p()

# ---------------------------------------------------------------- section 3
report.p("## Why the simulation cannot stand in for the machine")
report.p()
report.p("Shor's algorithm is the order mechanism (round 35's table) with one "
         "step no classical method has: it *reads* the order out of an "
         "interference pattern instead of guessing a multiple of it. The reading "
         "is where the exponential sits, and the measurement above is what it "
         "costs to fake: `N^2` amplitudes to hold the superposition the machine "
         "would hold for free.")
report.p()
report.p("Round 25 tested the obvious escape -- that the interference might be "
         "computing something a classical algorithm could compute directly -- by "
         "taking the spectrum of every classically cheap function of "
         "`a^x mod N` it could construct. All of them came back flat. The "
         "quantum step is not a classical function of anything cheap that this "
         "project could find, which is a negative result about this project and "
         "not a theorem.")
report.p()
report.p("So the honest summary of forty-five rounds is two sentences. "
         "Classically, every mechanism lands at `N^(1/4)` without smooth numbers "
         "and at `L[1/3]` with them, and nothing here moved either exponent. "
         "Quantumly the problem is polynomial and has been since 1994, and what "
         "stands between that and a factored 2048-bit modulus is twenty million "
         "physical qubits, which is an engineering number rather than a "
         "mathematical one.")
report.write()
