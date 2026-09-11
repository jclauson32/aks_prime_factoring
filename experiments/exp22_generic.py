"""Round 14: what generic algebra buys, which is nothing, and what every method adds.

Thirteen rounds produced barriers, and the same sentence kept appearing in the
write-ups: *every computable construction turned out symmetric in p and q*.
This round stops proposing algorithms and measures that sentence instead.

Factoring is the manufacture of a **zero divisor**.  Under CRT
`Z/N = F_p x F_q`, and a non-trivial `gcd(x, N)` is an `x` vanishing in exactly
one coordinate.  Ring operations act *diagonally* -- the same `+`, `-`, `*`
applied to both coordinates -- so no arrangement of them distinguishes the two.

The measurement below puts a number on that, and on what the working methods add.
"""

import random
from math import gcd

from _common import Report

from aksfactor.arith import is_prime
from aksfactor.generic import alignment_gain, generic_rate, zero_divisor_density

report = Report(
    "exp22_generic",
    "Generic algebra is worth exactly nothing; one alignment is worth everything",
    "Round 14. A measurement of the algorithm space, not another algorithm.",
)

rng = random.Random(5)


def randprime(lo, hi):
    while True:
        x = rng.randrange(lo, hi) | 1
        if is_prime(x):
            return x


report.p("## The generic baseline")
report.p()
report.p("A straight-line program over `Z/N`: start from a random base, apply "
         "`+`, `-`, `*` in any arrangement. Because `*` can square a register, "
         "depth `d` reaches exponents up to `2^d` -- **this model contains "
         "Pollard's `p-1`**. The question is what a *random* such program achieves.")
report.p()
rows = []
for lo, hi in [(113, 300), (353, 900), (3559, 9000)]:
    p = randprime(lo, hi)
    q = randprime(p + 2, 3 * p)
    n = p * q
    base = zero_divisor_density(n, p, q)
    for depth in (6, 12, 24):
        rate = generic_rate(n, depth, 120000, random.Random(depth * 7 + 1))
        rows.append([p, q, depth, f"{rate:.6f}", f"{base:.6f}", f"{rate / base:.3f}"])
report.table(["p", "q", "program depth", "measured split rate",
              "1 - phi(N)/N", "ratio"], rows)
report.p("The rate converges to `1 - phi(N)/N = 1/p + 1/q`: **exactly the chance "
         "that an element drawn uniformly at random from `Z/N` happens to be a "
         "zero divisor.** Depth buys enormous exponents and buys nothing else. "
         "Generic algebra is worth precisely as much as reaching into a hat.")
report.p()

report.p("## What one alignment is worth")
report.p()
p = randprime(3559, 9000)
q = randprime(p + 2, 3 * p)
n = p * q
base = zero_divisor_density(n, p, q)
rows = []
for bound in (50, 100, 200):
    got = alignment_gain(n, bound, 1500, random.Random(1))
    rows.append([bound, f"{got['rate']:.4f}", f"{base:.6f}",
                 f"{got['rate'] / base:,.0f}x"])
report.table(["lcm bound B", "Pollard p-1 split rate", "generic baseline",
              "gain"], rows)
report.p("Same ring. Same operations. Same depth budget -- `lcm(1..B)` is reached "
         "by repeated squaring like anything else. The **only** difference is that "
         "the exponent is chosen to be divisible by everything small, i.e. to "
         "align with the *order* of the group mod `p`.")
report.p()
report.p("That single choice is worth three to four orders of magnitude, and it is "
         "the entire content of the method. Generic algebra cannot make it: the "
         "order of the group mod `p` is not visible to `+`, `-`, `*` on `Z/N`.")
report.p()

report.p("## The taxonomy")
report.p()
report.p("Every working factoring method manufactures its zero divisor through one "
         "of a small number of mechanisms. **An earlier version of this experiment "
         "claimed there were only two, order and size. That was wrong** -- it "
         "omitted Pollard's rho and the whole quadratic/number field sieve family, "
         "and the latter beats both caps it listed. The corrected table:")
report.p()
report.p("| mechanism | how the zero divisor arises | methods | best known |")
report.p("|---|---|---|---|")
report.p("| **order** | exponent divisible by `\\|G\\|` for a group `G` attached to `p` | Pollard `p-1`, Williams `p+1`, ECM, class groups | `L[1/2]` (ECM) -- capped by smooth-number density |")
report.p("| **size** | a window or interval that contains `p` | trial division, Fermat, Lehman, Strassen, Coppersmith, Harvey | `N^(1/5)` deterministic (Harvey) |")
report.p("| **collision** | two iterates agreeing mod `p` but not mod `q` | Pollard rho | `N^(1/4)` -- the birthday bound on a set of size `p` |")
report.p("| **sign** | a square root whose CRT signs differ between `p` and `q` | Dixon, CFRAC, QS, NFS; the central column of round 15 | `L[1/3]` (NFS, heuristic) |")
report.p()
report.p("The sign mechanism is the strongest known, and it is worth being precise "
         "about *where* its cost lives. Given a congruence `x^2 = y^2 (mod N)` with "
         "`x != +-y`, the split is a single gcd. Finding the congruence is the hard "
         "part, and the sieves do it through smooth relations -- so their running "
         "time is set by smoothness density, just in a smaller-number setting than "
         "the order methods. By Rabin's reduction, producing square roots mod `N` "
         "on demand is exactly as hard as factoring.")
report.p()
report.p("Mapped onto the fourteen rounds: the Pascal row, the AKS fold, the "
         "q-deformation, the norm-one subgroup and the class group are **order**; "
         "the gasket threshold, the factorial search, Lehman's fan and "
         "Coppersmith's window are **size**. None of them was collision or sign "
         "-- which is part of why none of them approached `L[1/3]`.")
report.p()

report.p("## What a polynomial-time algorithm has to do")
report.p()
report.p("Every mechanism above has a super-polynomial cost, and in each one the "
         "cost sits in the same place: **finding** the structure that makes the "
         "zero divisor, not using it. A smooth group order, a window containing "
         "`p`, a collision, a mixed-sign square root -- once any of these is in "
         "hand, the factor falls out in polynomial time. A polynomial-time "
         "algorithm needs one of them to be *findable* in polynomial time, or a "
         "mechanism not on the list.")
report.p()
report.p("Every barrier in this repository is a statement about one of those "
         "finding steps:")
report.p()
report.p("- symmetric constructions (Theorem 8, the norm; round 13, `N mod ell`) "
         "find nothing, because they see `N` and not `p`;")
report.p("- aliased constructions (Theorem 7) cannot locate a window below the "
         "scale `p`;")
report.p("- rigid group orders (Proposition 14) offer one order per `p` and no "
         "way to redraw;")
report.p("- dense families (round 5) redraw freely and run into smoothness "
         "density instead.")
report.p()
report.p("This is not a proof that no polynomial-time mechanism exists. It is a "
         "measurement that fourteen rounds of algebraic, geometric, "
         "group-theoretic, fractal, lattice and analytic attempts produced none, "
         "and a question precise enough to judge the next attempt in a sentence: "
         "*which mechanism does it use, and why would its finding step be cheaper "
         "than the best known one?*")
report.write()
