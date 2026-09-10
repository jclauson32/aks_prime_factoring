"""The headline identity: first non-zero entry of row n sits at spf(n), value n/spf(n)."""

import random
import time

from _common import Report

from aksfactor.arith import is_prime, spf_trial
from aksfactor.pascal import row_entry
from aksfactor.theorems import check_t3_exact_value

UPTO = 3000

report = Report(
    "exp02_first_nonzero",
    "First non-zero residue: position = smallest prime factor, value = cofactor",
    "`min { k > 0 : C(n,k) != 0 mod n } = spf(n)` and the residue there is exactly "
    "`n / spf(n)`.",
)

fails = 0
checked = primes_checked = 0
for n in range(4, UPTO):
    if is_prime(n):
        # Theorem 4(1): the whole interior must vanish.
        primes_checked += 1
        if any(row_entry(n, k) for k in range(1, n)):
            fails += 1
        continue
    checked += 1
    p = spf_trial(n)
    # every position strictly before spf(n) must vanish ...
    if any(row_entry(n, k) for k in range(1, p)):
        fails += 1
    # ... and the entry at spf(n) must be exactly the cofactor.
    if row_entry(n, p) != n // p:
        fails += 1
    ok, _ = check_t3_exact_value(n)
    if not ok:
        fails += 1

report.p(f"- composite `n` checked: **{checked:,}**, prime `n` checked: "
         f"**{primes_checked:,}** (`n < {UPTO}`)")
report.p(f"- failures: **{fails}**")
report.p()
report.p("### `C(n, p) mod n = n/p` for every prime factor `p`, at scale")
report.p()
report.p("Random access (Lucas/Granville) evaluates a single residue without "
         "materialising the row, so the identity can be checked on numbers far "
         "beyond any tabulation.")
report.p()

rng = random.Random(20260909)
rows = []
for bits in (64, 128, 256, 512, 1024):
    q = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
    while not is_prime(q):
        q += 2
    p = 1000003
    n = p * q
    t0 = time.time()
    r = row_entry(n, p)
    dt = time.time() - t0
    rows.append([len(str(n)), p, f"{str(q)[:18]}...", "yes" if r == q else "NO",
                 f"{dt * 1e6:.0f} us"])
report.table(["digits of n", "p", "cofactor q", "C(n,p) mod n == q", "time"], rows)
report.p("A 300-digit cofactor drops out of one binomial residue in microseconds -- "
         "*provided you already know where to look*. Finding `p` is the hard part, "
         "and that is what `exp04` measures.")
report.write()
