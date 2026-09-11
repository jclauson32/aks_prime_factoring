"""Round 15: the central column of Pascal's triangle computes Legendre symbols.

Round 14's taxonomy (as corrected) lists four mechanisms by which factoring
methods manufacture zero divisors, and every earlier round of this project used
one of two -- order or size.  This round builds a Pascal-native method on the
third that had not been tried here: the **sign** of a square root.

    sum_k C(2k,k) x^k  =  (1 - 4x)^(-1/2)

and modulo p the truncated sum is the Legendre symbol ((1-4a)/p), for every
truncation point in a window spanning a factor of two.
"""

import random
import time
from math import comb, gcd

from _common import Report

from aksfactor.arith import is_prime
from aksfactor.central import (
    central_split,
    central_sum_naive,
    jacobi,
    power_series_naive,
)
from aksfactor.fast import fast_spf

report = Report(
    "exp23_central_column",
    "The central column: Legendre symbols from Pascal's triangle",
    "Round 15. A Pascal-native method on the sign mechanism -- and its cost.",
)

rng = random.Random(8)


def randprime(lo, hi):
    while True:
        x = rng.randrange(lo, hi) | 1
        if is_prime(x):
            return x


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


report.p("## The identity")
report.p()
report.p("For `0 <= k <= (p-1)/2`, `C(2k,k) = C((p-1)/2, k) (-4)^k (mod p)`; for "
         "`(p+1)/2 <= k <= p-1` Lucas makes `C(2k,k) = 0 (mod p)`. So for every "
         "`T` in `[(p-1)/2, p-1]`")
report.p()
report.p("```")
report.p("P_T(a) = sum_{k<=T} C(2k,k) a^k  ==  (1-4a)^((p-1)/2)  ==  ((1-4a)/p)   (mod p)")
report.p("```")
report.p()
bad = tot = 0
for p in (11, 13, 17, 19, 23, 29, 31, 37, 41, 43):
    for a in range(1, 25):
        leg = legendre(1 - 4 * a, p)
        if leg == 0:
            continue
        for t in range((p - 1) // 2, p):
            tot += 1
            bad += central_sum_naive(t, a, p) != leg % p
report.p(f"Checked on **{tot:,}** `(p, a, T)` triples: **{bad}** mismatches.")
report.p()

report.p("## Inside the window it always splits; below it, never")
report.p()
report.p("A first test of this method reported 12 of 12 semiprimes split -- all at "
         "`T = 2` or `4`, far below the window. That was round 1's trap again: with "
         "`p <= 120` and bases up to 40, the tiny polynomial `P_2(a) - 1 = "
         "2a(1+3a)` simply hits `p`, which is trial division. With `p` large enough "
         "that such coincidences are rare:")
report.p()
buckets = {"below the window": [0, 0], "inside the window": [0, 0]}
for _ in range(10):
    p = randprime(1500, 4000)
    q = randprime(p + 2, int(1.6 * p))
    n = p * q
    a = next(a for a in range(2, 500) if jacobi(1 - 4 * a, n) == -1)
    lo, hi = (p - 1) // 2, p - 1
    grid = list(range(10, lo, max(1, lo // 25))) + list(range(lo, hi + 1, max(1, (hi - lo) // 25)))
    for t in grid:
        s = central_sum_naive(t, a, n)
        split = any(1 < gcd((s - e) % n, n) < n for e in (1, n - 1))
        key = "below the window" if t < lo else "inside the window"
        buckets[key][0] += split
        buckets[key][1] += 1
report.table(["truncation T", "splits", "of", "rate"],
             [[k, h, t, f"{h / t:.1%}"] for k, (h, t) in buckets.items()])
report.p("So the mechanism is real. It is also, notably, a detection **below `p`**: "
         "the window `[(p-1)/2, p-1]` sits entirely under `p`, where "
         "`gcd(T!, N) = 1` and Strassen's factorial test (Theorem 20) sees nothing "
         "at all.")
report.p()

report.p("## End to end, against Strassen")
report.p()
report.p("`central_split` tries `T = 2, 4, 8, ...`; because the window spans a factor "
         "of two, one power of two always lands in it, so there is no scan. Each "
         "`P_T(a)` is a hypergeometric partial sum, evaluated in `O~(sqrt(T))` by "
         "baby-step/giant-step on 2x2 polynomial matrices -- the factorial "
         "algorithm, one dimension up.")
report.p()
rng2 = random.Random(21)


def randprime2(lo, hi):
    while True:
        x = rng2.randrange(lo, hi) | 1
        if is_prime(x):
            return x


rows = []
mech = {}
for lo, hi in [(10**3, 4 * 10**3), (10**4, 4 * 10**4), (10**5, 4 * 10**5),
               (10**6, 4 * 10**6), (10**7, 4 * 10**7)]:
    p = randprime2(lo, hi)
    q = randprime2(p + 2, int(1.8 * p))
    n = p * q
    stats = {}
    t0 = time.time()
    got = central_split(n, stats=stats)
    t_c = time.time() - t0
    t0 = time.time()
    fast_spf(n)
    t_s = time.time() - t0
    assert got == (p, q)
    mech[stats["mechanism"]] = mech.get(stats["mechanism"], 0) + 1
    rows.append([f"{p:,}", f"{stats['T']:,}", f"[{(p - 1) // 2:,}, {p - 1:,}]",
                 stats["mechanism"], f"{t_c:.2f}", f"{t_s:.2f}",
                 f"{t_c / max(t_s, 1e-9):.0f}x"])
report.table(["p", "T that split", "window", "mechanism", "central (s)",
              "Strassen (s)", "ratio"], rows)
report.p("Same exponent, `O~(sqrt(p))`, and a constant roughly an order of magnitude "
         "worse: each tree node multiplies 2x2 polynomial matrices, eight products "
         "where the factorial needs one. Splits the classifier labels "
         "*coincidence* are the partial sum hitting `+-1 (mod p)` by chance below "
         "the window -- legitimate, but trial-division-like, and not credited to "
         "the mechanism.")
report.p()

report.p("## The generalisation: d-th power characters")
report.p()
report.p("The central column is the `d = 2` case of `(1 - d^d x)^(-1/d)`, whose "
         "coefficients are integers for every `d`. When `d | p-1`, the same "
         "argument makes the truncated series equal to the **`d`-th power residue "
         "character** of `1 - d^d a` for every `T` in `[(p-1)/d, p-1]` -- so the "
         "detection threshold drops from `p/2` to `p/d`, and `gcd(P^d - 1, N)` "
         "splits. Measured, scanning `T` upward for the first split:")
report.p()
from aksfactor.arith import sieve  # noqa: E402

pool = [x for x in sieve(9000) if x > 2000]
rows = []
exact = 0
for d in (2, 3, 4, 6):
    picked = 0
    for p in pool:
        if (p - 1) % d:
            continue
        # q just above p keeps T < (q-1)/d, so mod q the series is not yet a
        # character and only p can fire -- no need to restrict q - 1 mod d
        q = next(x for x in pool if x > p + 50)
        n = p * q
        a = 2
        if (1 - d**d * a) % p == 0:
            continue
        first = None
        for t in range(max(1, (p - 1) // d - 40), p - 1):
            s = power_series_naive(t, a, n, d)
            if gcd((pow(s, d, n) - 1) % n, n) == p:
                first = t
                break
        if first is None:
            continue
        exact += first == (p - 1) // d
        rows.append([d, f"{p:,}", f"{(p - 1) // d:,}", f"{first:,}",
                     f"{first / p:.3f}", f"{1 / d:.3f}"])
        picked += 1
        if picked >= 2:
            break
report.table(["d", "p (d divides p-1)", "(p-1)/d", "first T that splits",
              "T / p", "1/d"], rows)
below = sum(1 for r in rows if int(r[3].replace(",", "")) < int(r[2].replace(",", "")))
if rows:
    msg = (f"The first split lands exactly at `T = (p-1)/d` in **{exact} of "
           f"{len(rows)}** cases.")
    if exact < len(rows):
        msg += (f" The {len(rows) - exact} exception(s) split *below* the window "
                f"({below} of them) -- the same ~2/p chance hit the classifier flags "
                "earlier, not the mechanism.")
    msg += (" So the cost falls from `sqrt(p/2)` to `sqrt(p/d)` -- a factor "
            "`sqrt(d/2)` -- *provided* `d` divides `p-1`.")
    report.p(msg)
report.p()
report.p("That proviso is the whole story. For a generic `p` one does not know "
         "which `d` divide `p-1`, and trying `D` candidates in parallel costs "
         "`D` series evaluations at each `T`, i.e. about `sqrt(p/d_max) * D`, which "
         "beats `sqrt(p)` only if `D < sqrt(d_max)` -- impossible when "
         "`d_max <= D`. The family therefore interpolates between the Legendre "
         "method (`d = 2`, works for every odd `p`) and Pollard `p-1` (large "
         "`d | p-1`, needs `p-1` to have a large known divisor), and at every point "
         "in between it pays for a large `d` with exactly the smoothness condition "
         "the order mechanism already needed.")
report.p()

report.p("## Verdict")
report.p()
report.p("A genuine addition, and not a breakthrough. The central column is the "
         "first method in this repository built on the **sign** mechanism, it is "
         "native to Pascal's triangle, and it detects `p` from below -- at `T ~ p/2`, "
         "where the factorial test used by every size method is blind. Its cost is "
         "Strassen's exponent with a worse constant. Its natural generalisation "
         "trades threshold for a divisibility condition on `p-1`, which walks it "
         "back into the order mechanism. The sign mechanism's strength, in the "
         "sieves, comes from finding congruences of squares through smooth "
         "relations; a single character evaluation per `T` does not tap that.")
report.write()
