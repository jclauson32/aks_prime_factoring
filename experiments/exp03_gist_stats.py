"""How often is the `x` in `residue = x * y` actually prime, and how often is y = 1?"""

from math import gcd

from _common import Report

from aksfactor.arith import is_prime, spf_trial
from aksfactor.pascal import row_prefix

UPTO = 600

report = Report(
    "exp03_gist_stats",
    "Statistics on the `x * y` split",
    "The conjecture said the remainder breaks into `x * y` with `x` a prime factor. "
    "`x` is *always* a proper divisor; here is how often it is also prime.",
)

nonzero = x_prime = y_one = gcd_nontrivial = 0
x_proper = x_divides_r = 0
first_x_prime = first_total = 0
for n in range(4, UPTO):
    if is_prime(n):
        continue
    row = row_prefix(n, n - 1)
    seen_first = False
    for k in range(1, n):
        r = row[k]
        if not r:
            continue
        x = n // gcd(n, k)
        nonzero += 1
        if 1 < x < n and n % x == 0:
            x_proper += 1
        if r % x == 0:
            x_divides_r += 1
        if is_prime(x):
            x_prime += 1
        if r == x:
            y_one += 1
        if 1 < gcd(r, n) < n:
            gcd_nontrivial += 1
        if not seen_first:
            seen_first = True
            first_total += 1
            if is_prime(x):
                first_x_prime += 1

report.table(
    ["property of a non-zero residue", "count", "share"],
    [
        ["`x = n/gcd(n,k)` is a proper divisor of `n`", f"{x_proper:,}",
         f"{x_proper / nonzero:.3%}"],
        ["`x` divides the residue", f"{x_divides_r:,}",
         f"{x_divides_r / nonzero:.3%}"],
        ["`gcd(residue, n)` is a proper divisor", f"{gcd_nontrivial:,}",
         f"{gcd_nontrivial / nonzero:.3%}"],
        ["`x` is prime", f"{x_prime:,}", f"{x_prime / nonzero:.3%}"],
        ["`y = 1` (residue is exactly `x`)", f"{y_one:,}", f"{y_one / nonzero:.3%}"],
    ],
)
report.p(f"At the **first** non-zero position specifically, `x = n/spf(n)` is prime for "
         f"**{first_x_prime}/{first_total} = {first_x_prime / first_total:.1%}** of "
         f"composites below {UPTO} -- exactly the semiprimes and prime squares, since "
         "`n/spf(n)` is prime iff `n` is a product of two primes.")
report.p()
report.p("So the conjecture is right in the form that matters: for a semiprime "
         "`n = p*q` the first non-zero remainder is **exactly `q`**, a prime factor, "
         "with `y = 1`.")
report.write()
