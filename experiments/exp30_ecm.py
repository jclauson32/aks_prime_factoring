"""Round 22: the order mechanism, redrawn -- Lenstra's ECM against one fixed group.

Rounds 1-5 and 9-10 put Pascal's triangle on the order mechanism: the row
(1+x)^N mod (N, f) lives in F_{p^d}*, whose order p^d - 1 is fixed by p.  That
is Pollard's p - 1 in disguise, and Proposition 14 said why it stalls: one
group, one ticket.  This experiment shows the same budget spent on groups that
can be redrawn.
"""

import math
import random
import statistics
import time

from _common import Report

from aksfactor.arith import is_prime
from aksfactor.collision import pascal_rho
from aksfactor.ecm import ecm, ecm_curve, is_b1b2_smooth, pminus1

report = Report(
    "exp30_ecm",
    "The order mechanism, redrawn: ECM against one fixed group",
    "Round 22. Same budget, one group per prime versus a fresh group per curve.",
)

rng = random.Random(30)
B1, B2 = 2000, 100000


def randprime(bits):
    while True:
        x = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if is_prime(x):
            return x


report.p("## 1. One ticket per prime, or one per curve")
report.p()
report.p(f"200 semiprimes `N = p q` with a 28-bit `p` and a 60-bit `q`. Pollard's "
         f"`p - 1` gets `B1 = {B1}`, `B2 = {B2}`: it succeeds when the order of its base "
         f"mod `p` -- a divisor of `p - 1` -- is smooth to those bounds, and a new base "
         f"only moves among those divisors. Each ECM curve gets the same bounds and "
         f"a fresh group order `p + 1 - t`.")
report.p()
trials = 200
pm1 = rescued = failures = 0
by_class = {"p - 1 is (B1, B2)-smooth": [0, 0, 0], "p - 1 is not": [0, 0, 0]}
curves_needed = []
for _ in range(trials):
    p, q = randprime(28), randprime(60)
    n = p * q
    got = pminus1(n, B1, B2)
    pm1 += got is not None
    if got is None:
        failures += 1
        rescued += any(pminus1(n, B1, B2, base=b) for b in (3, 5, 7, 11))
    key = "p - 1 is (B1, B2)-smooth" if is_b1b2_smooth(p - 1, B1, B2) else "p - 1 is not"
    c = 0
    while True:
        c += 1
        if ecm_curve(n, rng.randrange(6, n - 1), B1, B2):
            break
    curves_needed.append(c)
    by_class[key][0] += 1
    by_class[key][1] += got is not None
    by_class[key][2] += c == 1
rows = [["Pollard p - 1 (one fixed group)", f"{pm1} of {trials}",
         f"{rescued} of {failures} failures, with four new bases"],
        ["ECM, first curve", f"{sum(c == 1 for c in curves_needed)} of {trials}",
         f"all, within {max(curves_needed)} curves (median {statistics.median(curves_needed):.0f})"]]
report.table(["method, same B1 and B2", "succeeded", "retries that succeed"], rows)
rows = [[k, v[0], f"{v[1]} of {v[0]}", f"{v[2]} of {v[0]}"] for k, v in by_class.items()]
report.table(["primes where", "count", "p - 1 succeeded", "ECM's first curve succeeded"], rows)
rate_s = by_class["p - 1 is (B1, B2)-smooth"][2] / max(1, by_class["p - 1 is (B1, B2)-smooth"][0])
rate_n = by_class["p - 1 is not"][2] / max(1, by_class["p - 1 is not"][0])
report.p(f"ECM's first curve succeeds on {rate_s:.0%} of the primes where `p - 1` "
         f"is smooth and {rate_n:.0%} of those where `p - 1` is hopeless: the "
         f"curve's order has nothing to do with `p - 1`. Every group in the Pascal family -- the row mod "
         "`(N, f)`, the AKS fold, the norm-one torus, the q-deformation -- has "
         "order `p^d - 1` or a divisor, fixed by `p`; that is the whole reason they "
         "stalled.")
report.p()

report.p("## 2. The cost follows `p`, not `N`")
report.p()
rows = []
ecm_series = []
for pbits in (32, 40, 48, 56, 64):
    ecm_t, rho_t, curves = [], [], []
    for _ in range(3):
        p, q = randprime(pbits), randprime(160 - pbits)
        n = p * q
        lp = math.log(2 ** pbits)
        b1 = max(500, int(math.exp(0.7 * math.sqrt(lp * math.log(lp)))))
        st = {}
        t0 = time.perf_counter()
        g = ecm(n, b1=b1, b2=50 * b1, curves=5000, rng=rng, stats=st)
        ecm_t.append(time.perf_counter() - t0)
        curves.append(st["curves"])
        assert g in (p, q)
        if pbits <= 40:
            t0 = time.perf_counter()
            for x0 in range(3, 100):
                g, _ = pascal_rho(n, 2, x0)
                if g:
                    break
            rho_t.append(time.perf_counter() - t0)
            assert g in (p, q)
    ecm_series.append((pbits, statistics.median(ecm_t)))
    rows.append([pbits, 160, b1, f"{statistics.mean(curves):.0f}",
                 f"{statistics.median(ecm_t):.2f} s",
                 f"{statistics.median(rho_t):.2f} s" if rho_t else "--"])
report.table(["bits of p", "bits of N", "B1", "ECM curves (mean of 3)",
              "ECM time (median)", "Pollard rho (median)"], rows)
(b0, t0_), (b1_, t1_) = ecm_series[1], ecm_series[-1]
per8 = (t1_ / t0_) ** (8 / (b1_ - b0)) if t0_ > 0 else float("nan")
report.p(f"With `N` fixed at 160 bits -- beyond the pure-Python quadratic sieve "
         f"here -- both costs are set by `p`. Rho's grows as `sqrt p`, `2^4 = 16` "
         f"times per 8 bits; ECM's grew {per8:.1f} times per 8 bits between "
         f"{b0}- and {b1_}-bit `p`, the sub-exponential `L_p[1/2]` the order "
         f"mechanism reaches once its groups can be redrawn -- Proposition 17's "
         f"ceiling, attained.")
report.write()
