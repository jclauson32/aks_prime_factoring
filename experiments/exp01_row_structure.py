"""Exhaustively verify the structural theorems on every row up to a bound.

Checks, for every composite n in range and every interior k:
  T1  n/gcd(n,k) divides C(n,k) mod n            -- the "x * y" split
  T2  gcd(k,n) = 1  =>  C(n,k) = 0 mod n         -- support lives on shared factors
  T5  the closed form for the residue
  gist  gcd(residue, n) is always a proper divisor of n
"""

from _common import Report

from aksfactor.arith import is_prime
from aksfactor.theorems import (
    check_gist,
    check_t1_shape,
    check_t2_coprime_vanishing,
    check_t5_closed_form,
)

UPTO = 700

report = Report(
    "exp01_row_structure",
    "Row structure: every residue is `x * y` with `x` a proper divisor of `n`",
    "Exhaustive check over all interior positions of every row `n < %d`." % UPTO,
)

fails = []
rows = composites = interior = 0
for n in range(4, UPTO):
    rows += 1
    if is_prime(n):
        ok, detail = check_t2_coprime_vanishing(n)
        if not ok:
            fails.append(("prime interior nonzero", detail))
        continue
    composites += 1
    interior += n - 1
    for name, fn in (
        ("T1", check_t1_shape),
        ("T2", check_t2_coprime_vanishing),
        ("T5", check_t5_closed_form),
        ("gist", check_gist),
    ):
        ok, detail = fn(n)
        if not ok:
            fails.append((name, detail))

report.p(f"- rows examined: **{rows}** (`n = 4..{UPTO - 1}`)")
report.p(f"- composite rows: **{composites}**, interior positions checked: **{interior:,}**")
report.p(f"- counterexamples found: **{len(fails)}**")
report.p()
if fails:
    report.p("```")
    for name, detail in fails[:20]:
        report.p(f"{name}: {detail}")
    report.p("```")
else:
    report.p("T1, T2, T5 and the `x * y` conjecture hold at **every** position checked, "
             "and every prime row has an identically zero interior.")
report.write()
