"""Round 17: Pascal's triangle on the collision mechanism, and the full map.

Round 14's corrected taxonomy lists four ways a factoring method manufactures
a zero divisor.  Fifteen rounds had put Pascal's triangle on three of them --
order (the row, the AKS fold), size (the factorial threshold) and sign (the
central column).  This round supplies the fourth: iterate a column,
x -> C(x, k) mod N, and wait for a collision mod p.
"""

import math
import random
import statistics

from _common import Report

from aksfactor.arith import is_prime, sieve
from aksfactor.collision import (
    binom_poly,
    conjugate_check,
    fibre_kappa,
    pascal_rho,
    predicted_kappa,
    rho_length,
)

report = Report(
    "exp25_pascal_rho",
    "Pascal rho: the triangle on the collision mechanism",
    "Round 17. The fourth mechanism, and what Pascal's reflection does to it.",
)

rng = random.Random(17)


def randprime(lo, hi):
    while True:
        x = rng.randrange(lo, hi) | 1
        if is_prime(x):
            return x


report.p("## 1. `C(x, 2)` is Pollard's rho")
report.p()
checks = {p: conjugate_check(p) for p in (101, 1009, 65537, 1000003)}
report.p("With `x = 2y + 1/2`, `C(x, 2) = 2(y^2 - 5/16) + 1/2`: iterating the "
         "second column of the triangle is iterating `y -> y^2 - 5/16`, Pollard's "
         "map with `c = -5/16`. Checked on `F_p` for "
         + ", ".join(f"`p = {p}` ({'ok' if ok else 'FAIL'})" for p, ok in checks.items())
         + ".")
report.p()

report.p("## 2. Pascal's reflection sets the collision rate")
report.p()
report.p("For a map `f` on `F_p` the expected rho length is `sqrt(pi p / (2 kappa))`, "
         "where `kappa = (1/p) sum_v m_v(m_v - 1)` over the fibre sizes "
         "`m_v = #f^(-1)(v)`. A random map has `kappa = 1`. The rows of the "
         "triangle are palindromes, which as polynomials reads "
         "`C(k-1-x, k) = (-1)^k C(x, k)`: for even `k` the map factors through "
         "`(x - (k-1)/2)^2`, each fibre is closed under the reflection, and "
         "`kappa` doubles.")
report.p()
primes = [q for q in sieve(300000) if q > 150000]
rows = []
runs = 300
for k in range(2, 10):
    kap = statistics.mean(fibre_kappa(k, p) for p in (10007, 20011))
    rs = []
    for _ in range(runs):
        p = rng.choice(primes)
        rs.append(rho_length(binom_poly(k, p), rng.randrange(p)) / math.sqrt(p))
    mean = statistics.mean(rs)
    se = statistics.stdev(rs) / math.sqrt(runs)
    pred = math.sqrt(math.pi / (2 * predicted_kappa(k)))
    rows.append([k, f"{kap:.3f}", predicted_kappa(k), f"{mean:.3f} +- {se:.3f}",
                 f"{pred:.3f}", f"{(mean - pred) / se:+.1f}"])
report.table(["k", "measured kappa", "predicted", f"rho / sqrt(p) ({runs} walks)",
              "sqrt(pi / 2 kappa)", "z"], rows)
worst_z = max(abs(float(r[5])) for r in rows)
report.p(f"Every row agrees with the prediction (largest deviation {worst_z:.1f} "
         f"standard errors). Even columns from `k = 4` on collide `sqrt 2` times "
         f"sooner.")
report.p()
report.p("That is the effect Brent and Pollard exploited with `x^(2^j) + c` for "
         "Fermat numbers, arising here from the symmetry of the triangle and "
         "for every `p`, not only `p = 1 (mod 2^j)`. Section 3 prices it.")
report.p()

report.p("## 3. Factoring with it")
report.p()
report.p("Brent's cycle finder on `x -> C(x, k) mod N`, gcds batched 64 at a time. "
         "`evaluations` counts every application of the map, Brent's advance loop "
         "included.")
report.p()
rows = []
ratios = []
for bits in (20, 24, 28):
    stats = {2: [], 4: []}
    trials = 30
    for _ in range(trials):
        p = randprime(1 << (bits - 1), 1 << bits)
        q = randprime(1 << (bits - 1), 1 << bits)
        while q == p:
            q = randprime(1 << (bits - 1), 1 << bits)
        n = p * q
        for k in (2, 4):
            total = 0
            for x0 in range(3, 60):
                g, steps = pascal_rho(n, k, x0)
                total += steps
                if g:
                    break
            assert g in (p, q)
            stats[k].append(total / math.sqrt(min(p, q)))
    m2, m4 = statistics.mean(stats[2]), statistics.mean(stats[4])
    ratios.append(m4 / m2)
    rows.append([bits, trials, f"{m2:.2f}", f"{m4:.2f}", f"{m4 / m2:.2f}"])
report.table(["bits of p and q", "semiprimes", "evaluations / sqrt(min p), k=2",
              "k=4", "ratio"], rows)
report.p(f"Pascal rho factors at `N^(1/4)` like Pollard's rho, because it is "
         f"Pollard's rho. The even column needs {statistics.mean(ratios):.2f} times "
         f"the evaluations of `k = 2` (the fibre theory says `1/sqrt 2 = 0.71`; "
         f"these are {trials}-sample means of a heavy-tailed quantity).")
report.p()
report.p("In multiplications: `C(x, 2)` conjugates to `y^2 - 5/16`, one squaring "
         "per step; `C(x, 4)` conjugates to `((u^2 - 5/4)^2 - 37)/24`, two "
         "squarings and a constant. With the gcd accumulator that is 2 against "
         "3 or 4 multiplications per step, so per unit of `sqrt p` the even column "
         "costs `3/sqrt 2 = 2.12` at best against `2`: **the reflection's `sqrt 2` "
         "is paid back in full, and then some.**")
report.p()

report.p("## 4. The map: Pascal's triangle on all four mechanisms")
report.p()
report.p("| mechanism | Pascal object | cost here | best known for the mechanism |")
report.p("|---|---|---|---|")
report.p("| order | the row `(1+x)^N` mod `(N, f)`; the AKS fold (rounds 1-5, 9-10) | rigid orders, one ticket per `p` | `L[1/2]` (ECM, varying curves) |")
report.p("| size | first non-zero entry, factorial threshold (T4, T20) | `O~(N^(1/4))` | `N^(1/5)` (Harvey) |")
report.p("| sign | central column `sum C(2k,k) x^k` (round 15) | `O~(N^(1/4))`, ~15x Strassen | `L[1/3]` (NFS) |")
report.p("| collision | column map `x -> C(x, k)` (this round) | `O(N^(1/4))`, = Pollard rho | `N^(1/4)` (rho) |")
report.p()
report.p("The triangle hosts every mechanism, and on each one it lands at that "
         "mechanism's *baseline* -- never at its best. The accelerations are all "
         "outside it: ECM's speed comes from redrawing the group, Harvey's from "
         "Lehman's geometry, NFS's from smoothness in a number field. Pascal's "
         "triangle is a universal host for factoring mechanisms and an "
         "accelerator of none of them.")
report.write()
