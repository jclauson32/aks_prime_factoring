"""Round 4: deforming the period.  A parameter that finally varies -- over the wrong set.

Round 3 closed with a requirement: find a construction whose relevant parameter
*varies* at fixed `p`.  The classical Pascal row has none; Theorem 7 pins its
period to `p`.  The Gaussian binomial `[n,k]_q` has one.

Divisibility of `[n,k]_q` by `p` is governed by `d = ord_p(q)`, and `d` moves as
the base `q` moves.  That is a genuine deformation: the q-row can be non-trivial
at positions where the classical row is *provably* zero.

This measures whether it escapes.
"""

import random
from math import gcd

from _common import Report

from aksfactor.arith import is_prime, spf_trial
from aksfactor.pascal import row_entry
from aksfactor.qpascal import (
    q_clause_one_position,
    q_first_hit,
    q_lucas_divides,
    q_pascal_row,
    q_period,
    shift_split,
)

report = Report(
    "exp12_qdeformation",
    "The q-deformation: a period that varies, over the divisors of p-1",
    "Round 4. The first construction here with a tunable period -- and where it lands.",
)

rng = random.Random(23)


def randprime(lo, hi):
    while True:
        x = rng.randrange(lo, hi) | 1
        if is_prime(x):
            return x


def largest_pf(m):
    d, worst = 2, 1
    while d * d <= m:
        while m % d == 0:
            m //= d
            worst = max(worst, d)
        d += 1
    return max(worst, m)


def _prime_power_ladder(bound, cap=1 << 40):
    from aksfactor.arith import sieve

    for q in sieve(bound):
        pw = q
        while pw * q <= cap:
            pw *= q
        yield pw


def _smooth_prime():
    while True:
        v = 2
        while v < 10**6:
            v *= rng.choice([2, 3, 5, 7, 11, 13, 17, 19, 23])
        if is_prime(v + 1):
            return v + 1


def _rough_prime():
    while True:
        t = randprime(10**5, 10**6)
        if is_prime(2 * t + 1):
            return 2 * t + 1

# ---------------------------------------------------------------- criterion
report.p("## The q-analogue of Kummer's theorem")
report.p()
report.p("```")
report.p("p | [n,k]_q   <=>   (k mod d) > (n mod d)          [clause 1]")
report.p("               or    p | C(floor(n/d), floor(k/d))  [clause 2]")
report.p("```")
report.p()
report.p("with `d = ord_p(q)`. Clause 1 is a low-digit carry in base `d`; clause 2 is "
         "a classical base-`p` carry one level up.")
report.p()
bad = tot = 0
for p in (7, 11, 13, 17, 19, 23):
    for base in (2, 3, 5, 6, 7, 10):
        if base % p == 0:
            continue
        d = q_period(base, p)
        for n in (p * 11, p * 13, p * 17):
            if n > 400:
                continue
            row = q_pascal_row(n, min(40, n), base, p)
            for k in range(1, len(row)):
                tot += 1
                if (row[k] % p == 0) != q_lucas_divides(n, k, d, p):
                    bad += 1
report.p(f"Verified on **{tot:,}** cases: **{bad}** mismatches.")
report.p()

# ---------------------------------------------------- the deformation is real
report.p("## The deformation is real: firing below `spf(n)`")
report.p()
report.p("Theorem 2 says the classical row is *identically zero* at every "
         "`k < spf(n)`. The q-row is not.")
report.p()
rows = []
for p, r in [(11, 13), (13, 17), (17, 19), (19, 23), (23, 29), (29, 31)]:
    n = p * r
    s = spf_trial(n)
    assert all(row_entry(n, k) == 0 for k in range(1, s))
    best = None
    for q in range(2, 60):
        # exclude bases that give the factor away for free: gcd(q-1, n) > 1 means
        # ord_p(q) = 1, and a single gcd would already have split n.
        if gcd(q, n) != 1 or gcd(q - 1, n) != 1:
            continue
        hit = q_first_hit(n, s - 1, q)
        if hit and (best is None or hit[0] < best[0]):
            best = (hit[0], hit[1], q)
    if best:
        rows.append([n, f"{p} x {r}", s, "0 (all of them)", best[0], best[2], best[1]])
report.table(["n", "factors", "spf(n)", "classical row for k < spf", "q-row first hit k",
              "base q", "factor found"], rows)
report.p("Bases with `gcd(q-1, n) > 1` are excluded throughout: those have "
         "`ord_p(q) = 1`, and a single gcd would already have split `n`, so they "
         "prove nothing. What remains is genuine -- the q-deformation reaches "
         "inside the region Theorem 2 seals off, and round 3's requirement (a "
         "parameter that varies at fixed `p`) is met: `d = ord_p(q)` moves with "
         "`q`.")
report.p()

# ---------------------------------------------------------------- scaling
report.p("## But where does it land?")
report.p()
report.p("`shift_split` runs the same search at scale: clause 1 fires at "
         "`k = (n mod d) + 1`, so a hit means `ord_p(q) | n - j` for some small "
         "`j`, which is one modular exponentiation plus `j` gcds -- no row needed, "
         "any size `n`.")
report.p()
rows = []
for lo, hi in [(200, 600), (2000, 6000), (20000, 60000), (200000, 600000),
               (2 * 10**6, 6 * 10**6)]:
    p = randprime(lo, hi)
    r = randprime(p + 2, 3 * p)
    n = p * r
    budget = 64
    work = 0
    found = None
    for q in range(2, 4000):
        if gcd(q, n) != 1:
            continue
        if gcd(q - 1, n) != 1:
            continue
        work += budget
        got = shift_split(n, q, budget)
        if got:
            found = got
            break
    rows.append([f"{p:,}", f"{r:,}", f"{work:,}" if found else f">{work:,}",
                 f"{work / p:.2f}" if found else "-",
                 f"{largest_pf(p - 1):,}"])
report.table(["p", "r", "gcds until split", "work / p", "largest prime factor of p-1"],
             rows)
report.p("Work stays proportional to `p`. The deformation moves *where* the "
         "structure appears; it does not reduce how much of it you must read.")
report.p()

# ------------------------------------------------------- why: the dichotomy
report.p("## Why: the tuned period is 1 or huge, with nothing in between")
report.p()
report.p("A first guess -- mine -- was that this reduces to Pollard `p-1`, since "
         "`ord_p(q)` divides `p - 1`. That is wrong, and the data says so: with "
         "`p - 1` *smooth*, a random base still has order close to `p - 1`.")
report.p()
rows = []
for _ in range(3):
    p = _smooth_prime()
    ords = [q_period(q, p) for q in range(2, 9)]
    rows.append([f"{p:,}", f"{largest_pf(p - 1):,}",
                 ", ".join(f"{o / (p - 1):.2f}" for o in ords)])
report.table(["p with smooth p-1", "largest prime factor of p-1",
              "ord_p(q)/(p-1) for q = 2..8"], rows)
report.p("Smooth `p - 1` supplies many small *divisors*, but almost no elements of "
         "small *order*. Clause 1 needs a small order, not a smooth one, so raw "
         "q-deformation is **not** Pollard `p-1` -- it is weaker.")
report.p()
report.p("The natural fix is to tune the base: replace `q` by `Q = q^M` where `M` "
         "is a smooth prime-power ladder to bound `B`. Then "
         "`ord_p(Q) = ord_p(q) / gcd(ord_p(q), M)` is the `B`-rough part of the "
         "order. Measured over 60 random `(p, q)` per ladder:")
report.p()
rows = []
for B in (50, 200):
    M = 1
    for pw in _prime_power_ladder(B):
        M *= pw
    ones = large = mid = 0
    smallest_large = None
    for _ in range(60):
        p = randprime(10**5, 10**7)
        q = rng.randrange(2, p - 1)
        Q = pow(q, M, p)
        d = 1 if Q == 1 else q_period(Q, p)
        if d == 1:
            ones += 1
        elif d > B:
            large += 1
            smallest_large = d if smallest_large is None else min(smallest_large, d)
        else:
            mid += 1
    rows.append([B, ones, mid, large, f"{smallest_large:,}"])
report.table(["ladder bound B", "tuned period = 1", "1 < period <= B",
              "period > B", "smallest large period seen"], rows)
report.p("**The middle column is empty.** The tuned period is either `1` or larger "
         "than `B`, never usefully small. And `period = 1` means `Q = 1 (mod p)`, so "
         "`gcd(Q - 1, n)` has *already* split `n` -- that is Pollard `p-1` verbatim, "
         "with the q-row contributing nothing on top.")
report.p()
report.p("So the dichotomy is: tune the base and you either land exactly on Pollard "
         "`p-1`, or you are left with a period above `B` and clause 1 becomes a "
         "`1/B` lottery. There is no regime where the tunable period pays for "
         "itself.")
report.p()

report.p("## The hierarchy, updated")
report.p()
report.p("| construction | governing quantity | varies at fixed `p`? | over what set |")
report.p("|---|---|---|---|")
report.p("| classical Pascal row mod `n` | period `p` | no | -- |")
report.p("| **q-Pascal row mod `n`** | period `ord_p(q)` | **yes** | divisors of `p-1` |")
report.p("| norm-one subgroup, degree `d` | order `Phi_d(p)` | no | -- |")
report.p("| elliptic curve | order `p + 1 - t` | yes | an interval of width `4 sqrt(p)` |")
report.p()
report.p("The q-row is the first construction in this project whose parameter "
         "genuinely varies. It is still not enough, and the reason is now sharp: "
         "**the achievable periods are the divisors of `p-1`, and that set is "
         "sparse where it needs to be dense.** Elements of small order are almost "
         "nonexistent; tuning the base collapses the period to `1` (Pollard `p-1`) "
         "or leaves it above the ladder bound. An elliptic curve's order ranges "
         "over an *interval* of width `4 sqrt(p)` -- dense, so a fresh draw is "
         "always a genuinely fresh number. A divisor lattice is not.")
report.p()
report.p("That is the round-4 lesson, and it refines round 3's: it is not enough "
         "for the parameter to vary. **It has to vary over a dense set.**")
report.write()
