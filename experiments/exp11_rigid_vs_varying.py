"""Counting the group orders reachable at a *fixed* prime.

Round 3 ended on a structural claim: every method reachable from the AKS ring
has a **rigid** group order -- fix `p` and the number whose smoothness decides
success is fixed too -- whereas elliptic curves supply a *family* of orders at
the same `p`, which is what buys ECM its subexponential running time.

That claim is countable, so this counts it.
"""

import random
from math import isqrt

from _common import Report

from aksfactor.arith import is_prime, sieve
from aksfactor.grouporder import elliptic_order, partitions, ring_unit_order

report = Report(
    "exp11_rigid_vs_varying",
    "Rigid versus varying group orders, counted",
    "How many distinct group orders can each construction reach at one fixed `p`?",
)

rng = random.Random(31)

report.p("### The setup")
report.p()
report.p("Every method here succeeds on a prime `p` exactly when some integer "
         "attached to `p` is smooth. The question is how many *different* such "
         "integers a construction can offer at one fixed `p`.")
report.p()
report.p("- **AKS ring**: `(F_p[x]/f)*` has order `prod_i (p^deg(f_i) - 1)` over the "
         "irreducible factors of `f`. Varying `f` over monic degree `d` can only "
         "permute the multiset of factor degrees, so the count of distinct orders "
         "is capped by `partitions(d)` -- a constant, independent of `p`.")
report.p("- **Elliptic curves**: `#E(F_p) = p + 1 - t` with `|t| <= 2*sqrt(p)`, so "
         "the orders range over an interval of width `4*sqrt(p)`.")
report.p()

report.p("### Measured")
report.p()
rows = []
for p in (211, 503, 1009, 2003):
    ring_orders = {}
    for d in (2, 3, 4):
        seen = set()
        for _ in range(400):
            f = [rng.randrange(p) for _ in range(d)] + [1]
            order = ring_unit_order(f, p)
            if order is not None:
                seen.add(order)
        ring_orders[d] = seen
    ec = set()
    trials = 0
    while trials < 400:
        a, b = rng.randrange(p), rng.randrange(p)
        order = elliptic_order(a, b, p)
        trials += 1
        if order is not None:
            ec.add(order)
    rows.append([
        p,
        f"{len(ring_orders[2])} / {partitions(2)}",
        f"{len(ring_orders[3])} / {partitions(3)}",
        f"{len(ring_orders[4])} / {partitions(4)}",
        len(ec),
        f"{4 * isqrt(p)}",
    ])
report.table(["p", "ring d=2 (found/cap)", "ring d=3", "ring d=4",
              "distinct #E(F_p)", "Hasse width 4*sqrt(p)"], rows)
report.p("The ring columns saturate at `partitions(d)` and **stop**: 2 orders at "
         "`d = 2`, 3 at `d = 3`, 5 at `d = 4`, no matter how large `p` grows or "
         "how many `f` are tried. The elliptic column grows with `p`, tracking the "
         "Hasse interval.")
report.p()

report.p("### Why the cap is structural, not an artifact of sampling")
report.p()
report.p("`F_p[x]/f` for squarefree `f` is a product of finite fields "
         "`F_{p^{d_1}} x ... x F_{p^{d_m}}` with `sum d_i = d`. Its unit group "
         "order is therefore `prod_i (p^{d_i} - 1)`, a function of the *partition* "
         "`(d_1, ..., d_m)` alone. Two polynomials with the same degree pattern "
         "give literally the same order. There is nothing left to vary.")
report.p()
report.p("The same holds for any subgroup singled out by a norm condition: "
         "Theorem 13's norm-one subgroup has order `(p^d - 1)/(p - 1)`, again a "
         "function of `p` and `d` only. Restricting to a subgroup changes *which* "
         "rigid number you get, never the fact that it is rigid.")
report.p()

report.p("### What it costs")
report.p()
report.p("If `Phi_e(p)` is rough for every `e <= d`, then every construction in "
         "this repository fails at that `p`, and no additional work at that `p` "
         "helps -- you have exhausted the finitely many orders available. The only "
         "recourse is to raise `d`, and the cost of working in degree `d` grows "
         "while the supply of new orders grows only like `partitions(d)`.")
report.p()
report.p("ECM never faces that wall: a fresh curve is a fresh order, drawn from an "
         "interval of `~4*sqrt(p)` candidates, so one keeps sampling until a "
         "smooth order appears. Balancing sampling cost against smoothness density "
         "is exactly what produces `L[1/2]`.")
report.p()
report.p("### The gap, stated as a requirement")
report.p()
report.p("A polynomial-time method in this family would need a construction over "
         "`Z/n` whose group order at a fixed prime `p` varies with a parameter "
         "under our control. Products of the multiplicative groups of finite "
         "fields -- which is all that quotients of `(Z/n)[x]` supply -- cannot do "
         "it: their orders are pinned to `prod (p^{d_i} - 1)` by the degree "
         "partition. Escaping requires a genuinely different algebraic group, "
         "which is what an elliptic curve is and what no rearrangement of Pascal's "
         "triangle will produce.")
report.write()
