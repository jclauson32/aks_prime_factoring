"""Correcting round 2: the domination argument has a gap, and a real method lives in it.

Round 2 concluded that generalising Pollard `p-1` into the AKS ring is strictly
dominated, because `p-1` divides `p^k - 1` for every `k`, so `p^k - 1` is
`B`-smooth only if `p-1` already was.

That is sound **about the full unit group** and too strong in general.  The
norm-one subgroup of `F_{p^d}*` has order `(p^d - 1)/(p - 1)`, which `p - 1` does
not divide -- and you can land inside it without knowing `p`, by choosing a
monic polynomial whose roots multiply to `1`.  For `d = 2` that is
`x^2 - a x + 1`, tracked by the Lucas sequence `V_k(a,1)`: Williams' `p+1`
method, recovered as a special case.

So the framework does contain an asymmetric statistic after all.
"""

import random
import time

from _common import Report

from aksfactor.arith import is_prime
from aksfactor.cyclo import (
    lucas_v,
    lucas_v_binomial,
    norm_one_search,
    pollard_pminus1,
    williams_pplus1,
)

report = Report(
    "exp10_norm_one",
    "The norm-one subgroup: where the domination argument fails",
    "A correction to round 2, and a working method in the gap.",
)


def largest_pf(m):
    d, worst = 2, 1
    while d * d <= m:
        while m % d == 0:
            m //= d
            worst = max(worst, d)
        d += 1
    return max(worst, m)


rng = random.Random(11)
SMALL = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59]

report.p("## The gap in the argument")
report.p()
report.p("| group | order | is it divisible by `p-1`? |")
report.p("|---|---|---|")
report.p("| full unit group of `F_{p^d}` | `p^d - 1` | yes -- hence dominated |")
report.p("| **norm-one subgroup** | `(p^d - 1)/(p - 1)` | **no** |")
report.p()
report.p("For `d = 2` the norm-one order is `p + 1`. So the question is whether "
         "primes with `p-1` rough and `p+1` smooth actually exist in useful "
         "density.")
report.p()

# --- build adversarial primes: p-1 rough, p+1 smooth
targets = []
attempts = 0
while len(targets) < 6 and attempts < 400000:
    attempts += 1
    v = 2
    while v < 10**9:
        v *= rng.choice(SMALL)
    p = v - 1
    if not is_prime(p):
        continue
    if largest_pf(p - 1) > 10**4 and largest_pf(p + 1) <= 60:
        targets.append(p)
report.p(f"Sampling smooth `p+1` candidates, **{len(targets)}** such primes turned "
         f"up in {attempts:,} attempts. They are not exotic.")
report.p()

# --- a q that is rough on BOTH sides, so only p is reachable
q = None
while q is None:
    c = rng.randrange(10**9, 3 * 10**9) | 1
    if is_prime(c) and largest_pf(c - 1) > 10**5 and largest_pf(c + 1) > 10**5:
        q = c
report.p(f"### Head to head at equal budget (bound = 500)")
report.p()
report.p(f"Cofactor `q = {q:,}` is rough on both sides "
         f"(`lpf(q-1) = {largest_pf(q - 1):,}`, `lpf(q+1) = {largest_pf(q + 1):,}`), "
         f"so only `p` is reachable at this budget.")
report.p()

rows = []
won = 0
for p in targets:
    n = p * q
    t0 = time.time()
    got_pm1 = pollard_pminus1(n, bound=500)
    t_pm1 = time.time() - t0
    t0 = time.time()
    found = norm_one_search(n, bound=500, degree=2, bases=range(3, 60))
    t_pp1 = time.time() - t0
    got_pp1 = found[0] if found else None
    if got_pp1 == p and got_pm1 is None:
        won += 1
    rows.append([f"{p:,}", f"{largest_pf(p - 1):,}", largest_pf(p + 1),
                 "fail" if got_pm1 is None else f"{got_pm1:,}",
                 "fail" if got_pp1 is None else f"{got_pp1:,}",
                 found[1] if found else "-",
                 f"{t_pm1:.2f}/{t_pp1:.2f}"])
report.table(["p", "largest prime factor of p-1", "of p+1", "Pollard p-1",
              "norm-one (Williams p+1)", "base used", "secs p-1/p+1"], rows)
report.p(f"**{won} of {len(targets)}** are factored by the norm-one method and by "
         "nothing on the `p-1` side. The round-2 claim of strict domination is "
         "refuted.")
report.p()
report.p("Base retries are not optional: for `d = 2` the root only lands in "
         "`F_{p^2}` when `base^2 - 4` is a non-residue mod `p`, otherwise it falls "
         "back into `F_p` and the method silently degenerates to `p-1`. That is a "
         "coin flip per base, so a single failing base means nothing. (This "
         "experiment's first draft tried seven bases that all happened to be "
         "residues and wrongly recorded a failure.)")
report.p()

report.p("### The statistic is still a binomial sum")
report.p()
report.p("Worth noting where this leaves the Pascal framing. The Lucas sequence "
         "driving the method has a closed form")
report.p()
report.p("```")
report.p("V_k(a,1) = sum_j (-1)^j * (k/(k-j)) * C(k-j, j) * a^(k-2j)")
report.p("```")
report.p()
bad = 0
checks = 0
for n_mod in (10**9 + 7, 1001 * 9973, 97):
    for a in (3, 5, 7, 11):
        for k in range(0, 45):
            checks += 1
            if lucas_v(k, a, n_mod) != lucas_v_binomial(k, a, n_mod):
                bad += 1
report.p(f"Verified on **{checks:,}** cases, **{bad}** mismatches. So the statistic "
         "that finally breaks the symmetry between `p` and `q` *is* a binomial "
         "sum -- just not a sum along row `n`. What changes is not leaving "
         "binomials behind; it is that the roots are constrained to multiply to "
         "`1`, which is what pins the order to `p+1` instead of `p-1`.")
report.p()

report.p("### Degree `d > 2`")
report.p()
rows = []
for degree in (2, 3, 4):
    ok = tried = 0
    for _ in range(12):
        a = rng.choice([x for x in range(2000, 20000) if is_prime(x)])
        b = rng.choice([x for x in range(2000, 20000) if is_prime(x)])
        if a == b:
            continue
        m = a * b
        tried += 1
        found = norm_one_search(m, bound=150, degree=degree, bases=range(3, 12))
        if found and 1 < found[0] < m and m % found[0] == 0:
            ok += 1
    rows.append([degree, f"(p^{degree} - 1)/(p - 1)", f"{ok}/{tried}"])
report.table(["degree d", "norm-one order", "split rate, random semiprimes, bound 150"],
             rows)
report.p("Higher `d` covers `Phi_e(p)` for every `e | d` at once, because the "
         "defining polynomial usually factors mod `p` and the norm vanishes if "
         "*any* component's order divides `M`. More coverage per run, at higher "
         "cost per run.")
report.p()

report.p("### What this does and does not change")
report.p()
report.p("It **does** answer the question round 2 left open -- \"find a statistic "
         "of the first digit that is asymmetric in the prime factors\" -- in the "
         "affirmative. The norm-one construction is exactly such a statistic, and "
         "it is reachable without knowing `p`.")
report.p()
report.p("It does **not** give polynomial time, and the reason is worth stating "
         "precisely, because it is the next real obstruction:")
report.p()
report.p("| method | group | order | varies with a parameter we control? |")
report.p("|---|---|---|---|")
report.p("| Pollard `p-1` | `F_p*` | `p - 1` | no |")
report.p("| Williams / norm-one, degree `d` | `ker N` in `F_{p^d}*` | `Phi_d(p)`-ish | no -- fixed once `p` is |")
report.p("| ECM | `E(F_p)` | `p + 1 - t` | **yes** -- a new order per curve |")
report.p()
report.p("Every method here has a **rigid** group order: once `p` is fixed, so is "
         "the number whose smoothness decides success. You get one lottery ticket "
         "per `p`, and if `Phi_d(p)` is rough for all small `d` you are stuck. "
         "ECM's advantage is not a better group -- it is a *family* of groups, one "
         "per curve, so you can keep drawing fresh orders for the same `p`. That "
         "is what takes it from `Theta(smoothness luck)` to `L[1/2]`.")
report.p()
report.p("So the open problem sharpens again: **is there a family of AKS-ring-like "
         "objects whose group order varies with a parameter, at fixed `p`?** The "
         "norm-one construction gives one order per `(p, d)`; elliptic curves give "
         "unboundedly many per `p`. Nothing found here bridges that.")
report.write()
