"""The sharpest form of the barrier: what computation *would* suffice, and its cost.

Theorem 2 says the first n-adic digit of `C(n,k)` vanishes whenever
`gcd(k,n) = 1`.  So look at the *second* one:  `C(n,k) = n * m`.

Modulo a prime `p | n`, Lucas makes `m` vanish exactly when some base-`p` digit
of `k-1` exceeds the corresponding digit of `n-1`.  For a random `k` that is a
**constant-probability** event -- not the `1/p` lottery every other route in this
project ran into.

Which means: an algorithm computing `C(n,k) mod n^2` for random `k` in
`poly(log n)` would factor `n` in `O(1)` expected tries.  That is a hardness
result, not an algorithm -- and it pins down exactly which primitive the whole
approach is missing.
"""

import random
from math import comb, gcd

from _common import Report

from aksfactor.arith import is_prime
from aksfactor.binom import binom_mod_prime

report = Report(
    "exp09_second_digit",
    "The second n-adic digit, and what its hardness implies",
    "`C(n,k) = n*m` for `gcd(k,n) = 1`.  How often does `gcd(m, n)` split `n`?",
)

rng = random.Random(17)


def randprime(lo, hi):
    while True:
        x = rng.randrange(lo, hi) | 1
        if is_prime(x):
            return x


report.p("### The reduction")
report.p()
report.p("From `k*C(n,k) = n*C(n-1,k-1)` and `C(n,k) = n*m` we get `m = C(n-1,k-1)/k`. "
         "Since `gcd(k,n) = 1`, `k` is invertible, so for a prime `p | n`")
report.p()
report.p("```")
report.p("m = 0 (mod p)   <=>   C(n-1, k-1) = 0 (mod p)")
report.p("                <=>   some base-p digit of k-1 exceeds that of n-1   [Lucas]")
report.p("```")
report.p()
report.p("`gcd(m mod n, n)` is a proper factor exactly when that digit condition holds "
         "for some prime factors of `n` but not all.")
report.p()

report.p("### The Lucas shortcut is exact")
report.p()
mismatch = 0
checked = 0
for _ in range(3):
    p = randprime(11, 45)
    q = randprime(p + 2, 3 * p)
    n = p * q
    for k in range(2, min(n, 220)):
        if gcd(k, n) != 1:
            continue
        m = (comb(n, k) % (n * n)) // n
        zp = binom_mod_prime(n - 1, k - 1, p) == 0
        zq = binom_mod_prime(n - 1, k - 1, q) == 0
        checked += 1
        if (1 < gcd(m, n) < n) != (zp != zq):
            mismatch += 1
report.p(f"Checked **{checked:,}** positions against exact `C(n,k) mod n^2`: "
         f"**{mismatch}** mismatches.")
report.p()

report.p("### Hit rate versus the size of `p`")
report.p()
rows = []
budget = {0: 3000, 1: 3000, 2: 3000, 3: 2000, 4: 800, 5: 300}
for idx, (lo, hi) in enumerate([(11, 40), (120, 400), (1500, 4000),
                                (20000, 50000), (120000, 300000),
                                (600000, 1200000)]):
    p = randprime(lo, hi)
    q = randprime(p + 2, 3 * p)
    n = p * q
    digits, t = 0, n
    while t:
        digits += 1
        t //= p
    hits = trials = 0
    while trials < budget[idx]:
        k = rng.randrange(2, n)
        if gcd(k, n) != 1:
            continue
        zp = binom_mod_prime(n - 1, k - 1, p) == 0
        zq = binom_mod_prime(n - 1, k - 1, q) == 0
        trials += 1
        hits += zp != zq
    rows.append([f"{p:,}", f"{q:,}", digits, f"{hits / trials:.3f}",
                 f"{1 / p:.2e}", f"{(hits / trials) / (1 / p):,.0f}x"])
report.table(["p", "q", "digits of n base p", "hit rate", "1/p",
              "ratio to the 1/p lottery"], rows)
report.p("Sample sizes shrink as `p` grows for a reason worth noting: measuring "
         "the rate needs Lucas at *random* base-`p` digits, and one digit binomial "
         "`C(a,b) mod p` with `b ~ p/2` costs `O(p)` modular operations. Even "
         "*observing* this event is priced at the barrier.")
report.p()
report.p("The hit rate is **constant in `p`** -- it depends on the base-`p` digit "
         "profile of `n-1`, not on the size of `p`. Every other route in this "
         "project produced a `1/p` rate; this one is orders of magnitude above it "
         "and does not decay.")
report.p()

report.p("### What this means")
report.p()
report.p("**Proposition.** If `C(n,k) mod n^2` can be computed in `poly(log n)` for "
         "uniform random `k`, then `n = pq` can be factored in `O(1)` expected "
         "iterations. So that computation is factoring-hard.")
report.p()
report.p("This is where the project's barrier finally becomes *sharp*. Theorem 5 "
         "gives genuine random access to `C(n,k) mod n` -- for `n` of any size, in "
         "microseconds -- because the closed form only ever needs `C(n-1,k-1)` "
         "modulo the **small** number `k`, where Lucas and Granville apply. "
         "Extending that access by a single n-adic digit, to modulus `n^2`, would "
         "break factoring.")
report.p()
report.p("Two consequences worth stating:")
report.p()
report.p("1. The requirement, in Granville-style algorithms, that the modulus be "
         "*factored* before binomials can be reduced modulo it is essential -- not "
         "an artifact of how those algorithms are written.")
report.p("2. The reason this project stalls is not that the Pascal row hides the "
         "factors. It does not: the factors are sitting in the second n-adic digit "
         "at constant density. The obstruction is purely that reading that digit is "
         "the same problem as factoring.")
report.write()
