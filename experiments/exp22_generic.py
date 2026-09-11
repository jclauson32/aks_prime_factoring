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
report.p("Every working factoring method imports exactly one alignment with a "
         "structure attached to `p`, and only two such structures are known:")
report.p()
report.p("| alignment | what is aligned | methods | proven cap |")
report.p("|---|---|---|---|")
report.p("| **order** | exponent divisible by `\\|G\\|` for a group `G` attached to `p` | Pollard `p-1`, Williams `p+1`, ECM, class groups | `L[1/2]` -- the density of smooth numbers (Dickman) |")
report.p("| **size** | a window straddling the magnitude of `p` | Fermat, Lehman, Strassen, Coppersmith, Harvey | `N^(1/4)` by counting; `N^(1/5)` with Lehman + BSGS |")
report.p()
report.p("Each cap was measured or proved elsewhere in this repository: the "
         "smoothness ceiling in round 5, the counting bound in round 11, the "
         "`N^(1/5)` sweep in round 8 and its two closed escape routes in rounds 9 "
         "and 10.")
report.p()
report.p("And the fourteen rounds map onto the taxonomy exactly. The Pascal row, "
         "the AKS fold, the q-deformation, the norm-one subgroup and the class "
         "group are all **order** alignments. The gasket threshold, the factorial "
         "binary search, Lehman's fan and Coppersmith's window are all **size** "
         "alignments. Nothing examined here was anything else.")
report.p()

report.p("## What a polynomial-time algorithm has to do")
report.p()
report.p("It has to import a *third* alignment -- some structure attached to `p`, "
         "visible from `N` in polynomial time, that is neither the order of a "
         "group nor the magnitude of `p`. Every barrier in this repository is a "
         "consequence of there being only two:")
report.p()
report.p("- symmetric constructions (Theorem 8, the norm; round 13, `N mod ell`) "
         "see neither alignment and return information about `N` alone;")
report.p("- aliased constructions (Theorem 7) destroy the size alignment by "
         "sampling below the scale `p`;")
report.p("- rigid group orders (Proposition 14) offer one order alignment per `p` "
         "and no way to redraw;")
report.p("- dense families (round 5) offer unboundedly many order alignments, and "
         "run into the smoothness density instead.")
report.p()
report.p("This is not a proof that no third alignment exists. It is a measurement "
         "of the fact that fourteen rounds of algebraic, geometric, "
         "group-theoretic, fractal, lattice and analytic attempts produced no "
         "candidate for one -- and a statement precise enough that a fifteenth "
         "attempt can be judged in a sentence: *which structure attached to `p` "
         "does it align with, and is that structure order, size, or new?*")
report.write()
