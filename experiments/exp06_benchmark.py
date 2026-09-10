"""Timings: every Pascal route against conventional factoring.

The point of this table is *not* that one column wins. It is that all the
Pascal columns are `Theta(spf(n))`, and that the fastest of them is fast only
because Theorem 2 lets it replace a binomial residue with a gcd -- at which
point it has become trial division wearing a certificate.
"""

import random
import time

from _common import Report

from aksfactor.arith import factorize, is_prime, spf_trial
from aksfactor.factor import certificate, pascal_spf
from aksfactor.pascal import row_series

report = Report(
    "exp06_benchmark",
    "Benchmarks",
    "All Pascal routes are `Theta(spf(n))`; only the constants differ.",
)

rng = random.Random(7)


def randprime(lo, hi):
    while True:
        c = rng.randrange(lo, hi) | 1
        if is_prime(c):
            return c


report.p("### Finding `spf(n)` for balanced semiprimes")
report.p()
report.p("Four ways to reach the same answer:")
report.p()
report.p("1. **naive row scan** -- build the row prefix with `(1+x)^n mod (n, x^(k+1))` "
         "and look for the first non-zero. Assumes nothing, reads real residues.")
report.p("2. **residue scan** -- call `row_entry` at `k = 2, 3, 4, ...`. `row_entry` "
         "applies Theorem 2 and returns `0` immediately when `gcd(n,k) = 1`, so this "
         "loop *is* a gcd loop. That is the theorem doing the work, not a trick.")
report.p("3. **certified** -- 2-3-5 wheel to locate `p`, then one residue evaluation "
         "to emit the Pascal certificate `C(n,p) mod n = n/p`.")
report.p("4. **reference** -- ordinary trial division, and Pollard rho.")
report.p()

rows = []
for lo, hi in [(300, 900), (3000, 9000), (30000, 90000), (300000, 900000)]:
    p = randprime(lo, hi)
    q = randprime(p + 2, p * 2)
    n = p * q

    if p <= 5000:  # the naive scan is only affordable for small spf
        t0 = time.time()
        prefix = row_series(n, p)
        assert next(k for k in range(1, p + 1) if prefix[k]) == p
        t_naive = f"{time.time() - t0:.3f}"
    else:
        t_naive = "(too slow)"

    t0 = time.time()
    got_scan = pascal_spf(n, mode="scan")
    t_scan = time.time() - t0

    t0 = time.time()
    got_cert = pascal_spf(n, mode="certified")
    t_cert = time.time() - t0

    t0 = time.time()
    spf_trial(n)
    t_trial = time.time() - t0

    t0 = time.time()
    factorize(n)
    t_rho = time.time() - t0

    assert got_scan == got_cert == p, (n, got_scan, got_cert, p)
    assert certificate(n, p)["valid"]
    rows.append([p, t_naive, f"{t_scan:.5f}", f"{t_cert:.5f}", f"{t_trial:.5f}",
                 f"{t_rho:.5f}"])

report.table(["p = spf(n)", "1. naive row scan (s)", "2. residue scan (s)",
              "3. certified (s)", "4. trial division (s)", "Pollard rho (s)"], rows)
report.p("Columns 2, 3 and 4 converge, because Theorem 2 proves they are the same "
         "search. Column 1 shows the price of refusing to use the theorem. Pollard "
         "rho is `O(n^(1/4))` and pulls away as `p` grows -- the honest ceiling on "
         "this whole approach.")
report.p()
report.p("### Cost of reading a row prefix")
report.p()
rows = []
n = randprime(10**14, 10**15) * randprime(10**3, 10**4)
for kmax in (200, 800, 3200, 12800):
    t0 = time.time()
    row_series(n, kmax)
    rows.append([kmax, f"{time.time() - t0:.3f}"])
report.table(["kmax", "row_series time (s)"], rows)
report.p("`row_series` computes `(1+x)^n mod (n, x^(kmax+1))` in `O(log n)` truncated "
         "polynomial multiplications, Kronecker-packed into single big-integer "
         "products. It yields every residue up to `kmax` at once -- but the work "
         "still grows linearly in `kmax`, and by Theorem 2 the whole prefix is zero "
         "until `kmax` reaches `spf(n)`.")
report.write()
