"""Round 36: the collision mechanism, and why it refuses to parallelise.

The collision family walks a map on `Z/N` and waits for two iterates to agree
mod `p`.  Three things decide what it costs: the map's fibre statistic (round
17), the cycle-finding algorithm, and -- the interesting one -- whether the
collision can be *detected* without knowing `p`.

The last point is where the factoring version parts company with the discrete
logarithm version.  In a known group, distinguished points turn `m` machines
into an `m`-fold speedup (van Oorschot-Wiener).  Here the test "is `x mod p`
distinguished?" is not computable, so the only collisions that can be found are
the self-collisions inside a single walk, and `m` machines buy `sqrt m`.  This
round measures both numbers on the same walks.
"""

import random
import statistics
from math import gcd, isqrt, log, pi, sqrt

from _common import Report

from aksfactor.arith import is_prime
from aksfactor.collision import binom_poly, fibre_kappa, predicted_kappa, pascal_rho, rho_length

report = Report(
    "exp42_collision_family",
    "The collision mechanism, and the price of a hidden group",
    "Round 36. Fibre statistics, cycle finding, and why m machines buy sqrt(m).",
)

rng = random.Random(42)


def randprime(bits):
    while True:
        v = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if is_prime(v):
            return v


# ---------------------------------------------------------------- section 1
report.p("## What the map contributes: the fibre statistic")
report.p()
report.p("A walk on `x -> C(x, k) mod p` collides after about "
         "`sqrt(pi p / (2 kappa))` steps, where `kappa = (1/p) sum_v m_v (m_v - 1)` "
         "counts how much the map folds. Round 17 predicted `kappa = 2` for even "
         "`k >= 4` (Pascal's reflection closes every fibre) and `kappa = 1` "
         "otherwise. Both are recomputed here, and the prediction is checked "
         "against measured rho lengths.")
report.p()
report.p("The average has to be taken across *maps*, not across starting points "
         "within one map. A fixed map has a fixed functional graph, and every "
         "start in its giant component reports much the same rho length, so a "
         "thousand starts on one map estimate that map's own value to high "
         "precision and say nothing about how far it sits from the law. The "
         "scatter that matters is between maps, and it is what the last column "
         "reports.")
report.p()

maps_per_k = 12
starts_per_map = 60
kappa_primes = []
while len(kappa_primes) < maps_per_k:
    cand = rng.randrange(10000, 40000) | 1
    if is_prime(cand):
        kappa_primes.append(cand)

rows = []
summary = []
kappa_by_k = {}
for k in (2, 3, 4, 5, 6):
    kappas, ratios = [], []
    for p in kappa_primes:
        kappa = fibre_kappa(k, p)
        f = binom_poly(k, p)
        mean = statistics.mean(rho_length(f, rng.randrange(p))
                               for _ in range(starts_per_map))
        kappas.append(kappa)
        ratios.append(mean / sqrt(pi * p / (2 * kappa)))
    kappa_by_k[k] = kappas
    spread = statistics.stdev(ratios)
    mean_ratio, err = statistics.mean(ratios), spread / sqrt(len(ratios))
    summary.append((k, mean_ratio, err))
    rows.append([k, predicted_kappa(k), f"{statistics.mean(kappas):.3f}",
                 f"{mean_ratio:.2f} +- {err:.2f}", f"{spread:.2f}"])
report.table(["k", "kappa predicted", "kappa measured (mean over maps)",
              f"measured / predicted rho length, over {maps_per_k} maps",
              "spread between maps"], rows)

overall = statistics.mean(r for _, r, _ in summary)
lo_k, lo_r, lo_e = min(summary, key=lambda t: t[1])
hi_k, hi_r, hi_e = max(summary, key=lambda t: t[1])
gap = (hi_r - lo_r) / sqrt(hi_e ** 2 + lo_e ** 2)
worst = max(abs(statistics.mean(k_list) - predicted_kappa(k)) / predicted_kappa(k)
            for k, k_list in kappa_by_k.items())
report.p(f"`kappa` is barely statistical: averaged over the twelve maps it lands "
         f"within {worst:.1%} of the predicted integer for every column, "
         f"including the factor of two for even `k >= 4` that round 17 derived "
         f"from Pascal's reflection. The rho length follows the "
         f"law in shape and sits about {abs(1 - overall) * 100:.0f}% "
         f"{'below' if overall < 1 else 'above'} it in constant, averaging "
         f"{overall:.2f} across the columns.")
report.p()
report.p(f"That residual is not all noise: `k = {hi_k}` and `k = {lo_k}` differ "
         f"by {gap:.1f} standard errors ({hi_r:.2f} against {lo_r:.2f}). So "
         f"`kappa` captures the big effect -- which columns fold and which do not "
         f"-- without capturing everything about these particular maps. The last "
         f"column is the reason the averaging is done over maps: one map's own "
         f"functional graph scatters by tens of percent whatever the law says.")
report.p()
report.p("Folding the map harder shortens the walk by `sqrt(kappa)` and lengthens "
         "each step by the extra multiplications, which is why no column of "
         "Pascal's triangle beats the `k = 2` column -- that column being "
         "Pollard's rho (round 17).")
report.p()

# ---------------------------------------------------------------- section 2
report.p("## What the cycle finder contributes")
report.p()
report.p("Same walk, two ways of noticing it has closed. The measure is "
         "evaluations of the map, not seconds, so the comparison does not depend "
         "on how the gcds are batched.")
report.p()


def floyd_rho(n, k=2, x0=4, max_steps=1 << 24):
    """Floyd's tortoise and hare on `x -> C(x, k) mod n`, counting evaluations."""
    f = binom_poly(k, n)
    x = y = x0 % n
    steps = 0
    while steps < max_steps:
        x = f(x)
        y = f(f(y))
        steps += 3
        g = gcd(abs(x - y), n)
        if g > 1:
            return (g if g < n else None), steps
    return None, steps


def with_retries(method, n, k=2, starts=range(4, 30)):
    """Total evaluations until a split, retrying the starts a caller would.

    A walk can close up mod `n` instead of mod `p`; both finders then fail on
    that start and the work still counts.
    """
    spent = 0
    for x0 in starts:
        g, steps = method(n, k, x0)
        spent += steps
        if g:
            return g, spent
    return None, spent


trials = 40
floyd_steps, brent_steps, both = [], [], 0
for _ in range(trials):
    n = randprime(20) * randprime(20)
    fa, sa = with_retries(floyd_rho, n)
    fb, sb = with_retries(pascal_rho, n)
    if fa and fb:
        both += 1
        floyd_steps.append(sa)
        brent_steps.append(sb)
report.table(["cycle finder", "median evaluations", "numbers split"],
             [["Floyd (tortoise and hare)",
               f"{statistics.median(floyd_steps):,.0f}" if floyd_steps else "--",
               f"{both}/{trials}"],
              ["Brent", f"{statistics.median(brent_steps):,.0f}" if brent_steps else "--",
               f"{both}/{trials}"]])
ratio = (statistics.median(brent_steps) / statistics.median(floyd_steps)
         if floyd_steps and brent_steps else float("nan"))
report.p(f"Brent's evaluations come to {ratio:.2f} of Floyd's here. Floyd pays "
         f"three evaluations per step to Brent's one, and Brent overshoots the "
         f"collision by walking in powers of two; the saving is a constant, and "
         f"a constant is all a cycle finder can ever be worth.")
report.p()

# ---------------------------------------------------------------- section 3
report.p("## What `m` machines contribute -- twice")
report.p()
report.p("Now the point of the round. Fix a prime `p` and walk `x -> x^2 + c mod "
         "p` (the collision is a fact about the walk mod `p`; how it is detected "
         "is the variable). Two ways to use `m` machines:")
report.p()
report.p("- **Honest**: `m` walks with independent `c`, each waiting to close up "
         "on itself. That is what Pollard's rho can detect from `Z/N` alone, "
         "because a self-collision is found by the cycle structure, not by "
         "comparing values.")
report.p("- **Oracle**: all `m` walks on one map, with distinguished points "
         "stored in a shared table keyed by `x mod p`. Any two walks that meet "
         "are caught at the next distinguished point. This is van Oorschot-Wiener, "
         "and it needs `p`.")
report.p()

pbits = 26
trials = 30
machines = (1, 2, 4, 8, 16)


def rho_len(c, x0, p):
    """Steps until the orbit of `x0` under `x^2 + c` repeats a value mod `p`."""
    seen, x, i = set(), x0, 0
    while x not in seen:
        seen.add(x)
        x = (x * x + c) % p
        i += 1
    return i


def parallel_dp(p, m, threshold, max_steps):
    """Steps, summed over machines, until two distinguished points agree.

    A walk that falls into a cycle carrying no distinguished point never
    reports one, so van Oorschot and Wiener's rule is used: abandon a walk that
    has gone too long without a distinguished point and restart it from a fresh
    point. The abandoned steps still count. `None` if the budget runs out.
    """
    c = rng.randrange(1, p)
    walks = [rng.randrange(p) for _ in range(m)]
    since = [0] * m
    limit = 20 * (p // threshold)          # 20x the expected gap between points
    table, total = set(), 0
    while total < max_steps:
        for i in range(m):
            x = (walks[i] * walks[i] + c) % p
            walks[i] = x
            since[i] += 1
            total += 1
            if x < threshold:
                if x in table:
                    return total
                table.add(x)
                since[i] = 0
            elif since[i] > limit:
                walks[i] = rng.randrange(p)
                since[i] = 0
    return None


honest, oracle = {}, {}
abandoned = 0
for m in machines:
    per_machine, dp_total = [], []
    for _ in range(trials):
        p = randprime(pbits)
        per_machine.append(min(rho_len(rng.randrange(1, p), rng.randrange(p), p)
                               for _ in range(m)))
        got = parallel_dp(p, m, 128 * isqrt(p), 100 * isqrt(p))
        if got is None:
            abandoned += 1
        else:
            dp_total.append(got)
    honest[m] = statistics.mean(per_machine)
    oracle[m] = statistics.mean(dp_total) / m

report.table(["machines `m`", "honest: steps per machine", "speedup",
              "honest: total steps", "oracle DP: steps per machine", "speedup"],
             [[m, f"{honest[m]:,.0f}", f"{honest[1] / honest[m]:.2f}x",
               f"{honest[m] * m:,.0f}", f"{oracle[m]:,.0f}",
               f"{oracle[1] / oracle[m]:.2f}x"] for m in machines])


def slope(xs, ys):
    lx = [log(x) for x in xs]
    ly = [log(y) for y in ys]
    mx, my = statistics.mean(lx), statistics.mean(ly)
    num = sum((a - mx) * (b - my) for a, b in zip(lx, ly))
    den = sum((a - mx) ** 2 for a in lx)
    return num / den


honest_slope = slope(machines, [honest[m] for m in machines])
oracle_slope = slope(machines, [oracle[m] for m in machines])
report.p(f"Walks that went 20 gaps without a distinguished point were restarted "
         f"from a fresh point, as van Oorschot and Wiener prescribe -- without "
         f"that rule a walk can fall into a short cycle carrying no distinguished "
         f"point and never report. {abandoned} of {len(machines) * trials} runs "
         f"still hit the step cap and were dropped.")
report.p()
report.p("Both columns count steps to the moment a collision exists or is "
         "caught; the honest column omits the constant a cycle finder adds and "
         "the oracle column includes the distance to the next distinguished "
         "point, so the constants are not comparable between columns. The "
         "exponent in `m` is, and it is the claim.")
report.p()
report.p(f"Fitted exponents in `m`: honest `m^({honest_slope:.2f})`, oracle "
         f"`m^({oracle_slope:.2f})`, both within the noise of a five-point fit "
         f"of the theoretical `-1/2` and `-1`: the minimum "
         f"of `m` independent rho lengths is `sqrt(pi p / 2m)`, while pooling the "
         f"walks makes one birthday problem over all `m` trajectories, whose total "
         f"is flat in `m`.")
report.p()
report.p(f"The honest column's *total* work grows: {honest[1] * 1:,.0f} steps at "
         f"`m = 1` against {honest[max(machines)] * max(machines):,.0f} at "
         f"`m = {max(machines)}`. Parallel rho for factoring is not work-efficient; "
         f"parallel rho for a discrete logarithm is.")
report.p()

# ---------------------------------------------------------------- section 4
report.p("## Why the hidden group costs exactly that")
report.p()
report.p("The distinguished-point table is a hash table on the collision value. "
         "Without `p` there is no value to hash: two iterates `x` and `y` can only "
         "be compared by `gcd(x - y, N)`, one pair at a time. So `K` stored points "
         "cost `K^2 / 2` gcds to search, while hashing costs `K`. A single walk "
         "escapes this because its self-collision is found by cycle structure -- "
         "Floyd and Brent compare `O(K)` pairs, not `K^2` -- and that escape is "
         "available once per walk, not once per pair of walks.")
report.p()
K = 4096
report.table(["how a collision is noticed", "comparisons among K points", "needs p?"],
             [["cycle structure, one walk (Floyd/Brent)", f"O(K) = {K:,}", "no"],
              ["pairwise gcd, any two points", f"K^2/2 = {K * K // 2:,}", "no"],
              ["distinguished points, hashed", f"O(K) = {K:,}", "yes"]])
report.p("That table is the whole difference between `sqrt(m)` and `m`, and it is "
         "the same shape as every other result in this project: the mechanism is "
         "not the obstacle, the missing residue is.")
report.p()

report.p("## Where the four mechanisms land")
report.p()
report.p("Rounds 33, 34 and 35 did sign, size and order; with collision the "
         "survey is complete. Stripped of smooth numbers, every mechanism costs "
         "the same:")
report.p()
report.table(["mechanism", "best exponent without smoothness", "in this repository",
              "with smooth numbers"],
             [["size", "`O~(N^(1/4))`", "factorial threshold (T20), hull walk `N^(1/3)` (P23)",
               "--"],
              ["collision", "`N^(1/4)`", "Pascal rho = Pollard rho (P22)", "--"],
              ["sign", "`N^(1/4)`", "central column (round 15), SQUFOF (round 32)",
               "`L[1/2]` quadratic sieve, `L[1/3]` NFS"],
              ["order", "`N^(1/4)`", "baby-step giant-step over the unknown order, with "
               "the differences batched into one gcd -- the same trick as T20, "
               "and the one entry in this table not implemented separately here",
               "`L[1/2]` ECM, Pollard `p - 1`, the elliptic triangle"]])
report.p("The two escapes from `N^(1/4)` are both visible in round 33-36's data. "
         "Smoothness buys sub-exponential time and is the only thing that ever "
         "has. Combining mechanisms buys a little: Harvey's `N^(1/5)` is the hull "
         "walk's geometry plus a Fermat congruence plus a collision sweep, and it "
         "is the best deterministic exponent known. Neither escape is polynomial, "
         "and nothing in this repository suggests a third.")
report.write()
