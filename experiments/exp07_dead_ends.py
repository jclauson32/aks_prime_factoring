"""A scorecard of theories tried against the Theta(spf(n)) barrier, and how each died.

Every entry here was a genuine candidate for breaking out of trial-division
complexity. Recording *why* each failed is the useful part: the failures are
not accidents, they share a mechanism.
"""

from math import gcd

from _common import Report

from aksfactor.arith import is_prime
from aksfactor.ring import aks_pow, fold_norm, fold_norm_expected
from aksfactor.theorems import check_t7_aliasing

report = Report(
    "exp07_dead_ends",
    "Dead ends, and the mechanism they share",
    "Four attempts to beat `Theta(spf(n))`, each tested and each killed.",
)

# ---------------------------------------------------------------- theory A
report.p("## A. Symmetric functions of the fold (norms, resultants)")
report.p()
report.p("A single fold coefficient leaks a factor with probability `~1/p`. The "
         "*norm* `prod_t F(w^t)` aggregates all of them at once, so it looked like "
         "a way to buy every ticket with one purchase.")
report.p()
rows = []
ok = True
for n in (15, 35, 77, 143, 105, 1155, 27, 49, 221):
    for r in (3, 5, 7):
        for a in (1, 2, 3):
            got, want = fold_norm(n, r, a), fold_norm_expected(n, r, a)
            if got != want:
                ok = False
            if r == 5 and a == 2:
                rows.append([n, r, a, got, want, got == want])
report.table(["n", "r", "a", "Norm((x+a)^n)", "(a^r + 1)^n mod n", "equal"], rows)
report.p(f"Verified for every `(n, r, a)` tested (odd `r`): **{ok}**.")
report.p()
report.p("**Why it dies.** The norm equals `(a^r + 1)^n mod n` -- a function of "
         "`n`, `r`, `a` alone. It takes the same value modulo *every* prime factor "
         "of `n`, so `gcd(Norm - expected, n) = n` identically. Aggregating the "
         "tickets symmetrically destroys exactly the asymmetry that a factor "
         "consists of.")
report.p()
report.p("**The general lesson.** Any quantity invariant under the Galois action "
         "that permutes the `r`-th roots of unity is Frobenius-invariant mod every "
         "`p | n`, hence determined by `n`. Factoring needs a statistic that "
         "distinguishes `p` from `q`; symmetric ones cannot.")
report.p()

# ---------------------------------------------------------------- theory B
report.p("## B. Twisted moduli `x^r - c`")
report.p()
report.p("If the plain fold gives each coefficient a `1/p` chance, maybe weighting "
         "the wrapped terms by powers of `c` biases some coefficient toward "
         "vanishing.")
report.p()


def twist_pow(n, r, c, a):
    def mul(u, v):
        out = [0] * r
        for i, ui in enumerate(u):
            if ui:
                for j, vj in enumerate(v):
                    if vj:
                        s = i + j
                        if s < r:
                            out[s] = (out[s] + ui * vj) % n
                        else:
                            out[s - r] = (out[s - r] + ui * vj * c) % n
        return out

    base = [0] * r
    base[0] = a % n
    base[1 % r] = (base[1 % r] + 1) % n
    res = [0] * r
    res[0] = 1 % n
    e = n
    while e:
        if e & 1:
            res = mul(res, base)
        e >>= 1
        if e:
            base = mul(base, base)
    return res


import random  # noqa: E402

rng = random.Random(11)


def randprime(lo, hi):
    while True:
        x = rng.randrange(lo, hi) | 1
        if is_prime(x):
            return x


rows = []
for lo, hi in [(100, 300), (400, 900), (1200, 2500), (4000, 8000)]:
    p = randprime(lo, hi)
    q = randprime(p + 2, p * 4)
    n = p * q
    h1 = t1 = h2 = t2 = 0
    for r in range(2, 26):
        for a in (1, 2, 3):
            for c in (1, 2, 3, 5, 7):
                co = aks_pow(n, r, a) if c == 1 else twist_pow(n, r, c, a)
                for v in co:
                    hit = 1 < gcd(v, n) < n
                    if c == 1:
                        t1 += 1
                        h1 += hit
                    else:
                        t2 += 1
                        h2 += hit
    rows.append([p, q, f"{h1 / t1:.5f}", f"{h2 / t2:.5f}", f"{1 / p:.5f}"])
report.table(["p", "q", "plain `x^r-1` rate", "twisted `x^r-c` rate", "1/p"], rows)
report.p("Both track `1/p`. The twist changes which tickets you hold, not the "
         "price. Dead.")
report.p()

# ---------------------------------------------------------------- theory C
report.p("## C. The AKS-ring generalisation of Pollard `p-1`  *(partly retracted)*")
report.p()
report.p("In `F_p[x]/(x^r - 1) = prod_i F_{p^{k_i}}`, an element's order divides "
         "`lcm_i (p^{k_i} - 1)`. Raising `x + a` to `lcm(1..B)` and taking a gcd "
         "generalises Pollard `p-1`, which is the `k = 1` case. More `k` values, "
         "more chances -- so it seemed strictly better.")
report.p()


def largest_prime_factor(m):
    d, worst = 2, 1
    while d * d <= m:
        while m % d == 0:
            m //= d
            worst = max(worst, d)
        d += 1
    return max(worst, m)


rows = []
for p in (101, 1009, 10007, 65537, 99991):
    bs = [largest_prime_factor(p**k - 1) for k in (1, 2, 3)]
    rows.append([p, bs[0], bs[1], bs[2], "no" if min(bs) >= bs[0] else "yes"])
report.table(["p", "smoothness of p-1", "of p^2-1", "of p^3-1", "larger k helps?"],
             rows)
report.p("**Why it dies:** `p - 1` divides `p^k - 1` for every `k`, so `p^k - 1` is "
         "`B`-smooth only if `p - 1` already was. Working with `x + a` in the full "
         "unit group is dominated by the classical `p-1` method it generalises.")
report.p()
report.p("> **Correction (round 3).** The sentence that used to stand here said "
         "*strictly dominated*, full stop. That was too strong, and "
         "[exp10](exp10_norm_one.md) refutes it. The argument above is about the "
         "**full unit group**. The norm-one subgroup of `F_{p^d}*` has order "
         "`(p^d - 1)/(p - 1)`, which `p - 1` does **not** divide -- and you can "
         "land in it without knowing `p`, by choosing a monic polynomial whose "
         "roots multiply to `1`. For `d = 2` that is `x^2 - a x + 1`, i.e. Lucas "
         "sequences, i.e. Williams `p+1`. Primes with `p-1` rough and `p+1` smooth "
         "are plentiful, and the norm-one method factors them while `p-1` cannot. "
         "This entry is a dead end only for the unconstrained element `x + a`.")
report.p()

# ---------------------------------------------------------------- theory D
report.p("## D. Reading position out of the fold")
report.p()
report.p("If some residue class mod `r` were *occupied differently* by the support "
         "of the row, the class index would leak `p mod r`, and CRT over several "
         "`r` would reconstruct `p` in polylog time. This is the one that would "
         "actually have been polynomial.")
report.p()
rows = []
for n in (35, 77, 143, 221, 1001):
    for r in (5, 7, 9):
        if gcd(r, n) != 1:
            continue
        ok, detail = check_t7_aliasing(n, r)
        rows.append([n, r, str(detail["counts"]), detail["spread"], ok])
report.table(["n", "r", "support occupancy per class mod r", "spread", "even"], rows)
report.p("**Why it dies.** If `gcd(r, p) = 1` then `i -> i*p mod r` is a bijection "
         "of `Z/r`, so the multiples of `p` fall into every class equally. The "
         "occupancy is flat by *theorem*, not by luck. Positional information about "
         "`p` requires `gcd(r, p) > 1`, i.e. `r >= p` -- an aliasing barrier exactly "
         "like Nyquist: you cannot resolve a period-`p` signal by sampling it into "
         "fewer than `p` bins.")
report.p()

report.p("## The shared mechanism")
report.p()
report.p("All four die the same way. Modulo `p` the AKS object is "
         "`(x^(p^v) + a)^(n/p^v)`, whose only distinguishing feature is a period of "
         "`p^v`. Every cheap thing you can compute from it is either")
report.p()
report.p("- **symmetric** in the prime factors (A, and C for the unconstrained "
         "element) -- same value mod every `p`, so the gcd is `n`; or")
report.p("- **aliased** below the period (B, D) -- the period-`p` structure spreads "
         "evenly over all `r < p` buckets, leaving only the accidental vanishing of "
         "a bucket sum, a `1/p` event.")
report.p()
report.p("Extracting `p` is *period finding*. Classically that needs `Omega(p)` "
         "samples; it is the same problem Shor's algorithm solves in polylog time "
         "quantumly, by taking a Fourier transform of size `n` rather than size "
         "`r << p`. That is a sharp statement of what this framework is missing, "
         "and it is not something a cleverer choice of `r` or `a` can supply.")
report.write()
