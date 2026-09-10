"""Round 5: the density requirement is satisfiable -- and satisfying it caps out at L[1/2].

Round 4 ended with a precise target: a construction over `Z/n` whose governing
quantity at fixed `p` ranges over a *dense* set -- an interval, not a divisor
lattice -- while staying computable without knowing `p`.

Class groups of imaginary quadratic orders meet it.  For a discriminant
`D = -k*n`, the class number `h(D)` sits near `sqrt(|D|)` and moves essentially
arbitrarily with `k`.  Ambiguous forms (order dividing 2) *are* factorizations of
`D`, so a smooth `h(D)` hands back a factor of `n`.  That is Schnorr-Lenstra, and
it is implemented here so the contrast is measurable rather than asserted.
"""

import random
import time
from math import isqrt

from _common import Report

from aksfactor.arith import is_prime
from aksfactor.classgroup import (
    class_number_bruteforce,
    compose,
    identity_form,
    pow_form,
    reduce_form,
    schnorr_lenstra_split,
)
from aksfactor.grouporder import elliptic_order, partitions, ring_unit_order

report = Report(
    "exp13_classgroup",
    "Class groups: a dense family, and the ceiling it runs into",
    "Round 5. The density requirement is met -- and meeting it buys L[1/2], not P.",
)

rng = random.Random(2)

# --------------------------------------------------------------- correctness
report.p("## The group law, verified")
report.p()
bad = tested = 0
discs = [d for d in range(-3, -800, -1) if d % 4 in (0, 1)]
for d in discs[:120]:
    h = class_number_bruteforce(d)
    idf = identity_form(d)
    forms = []
    for a in range(1, isqrt(-d // 3) + 2):
        for b in range(-a + 1, a + 1):
            num = b * b - d
            if num % (4 * a):
                continue
            c = num // (4 * a)
            from math import gcd

            if c < a or gcd(gcd(a, b), c) != 1 or (a == c and b < 0):
                continue
            forms.append((a, b, c))
    assert len(forms) == h, (d, len(forms), h)
    for f in forms:
        tested += 1
        bad += compose(f, idf, d) != f
        bad += pow_form(f, h, d) != idf
        bad += compose(f, reduce_form(f[0], -f[1], f[2]), d) != idf
    for _ in range(15):
        x, y, z = (rng.choice(forms) for _ in range(3))
        bad += compose(compose(x, y, d), z, d) != compose(x, compose(y, z, d), d)
report.p(f"Identity, inverse, `f**h == 1`, and associativity over **120** "
         f"discriminants and **{tested}** forms: **{bad}** failures.")
report.p()
report.p("One bug is worth recording because it is invisible to casual testing: an "
         "extended-gcd returning a *negative* gcd makes the change-of-basis matrix "
         "have determinant `-1`. That is an **improper** equivalence -- it lands in "
         "a different class, so composition silently computes in the wrong group. "
         "Every structural check above failed until the sign was normalised.")
report.p()

# --------------------------------------------------------------- it factors
report.p("## It factors")
report.p()
rows = []
for n in (91, 187, 391, 667, 1147, 2021, 5183, 10403, 11021, 1009 * 1013):
    t0 = time.time()
    got = schnorr_lenstra_split(n, multipliers=range(1, 30), bound=200,
                                forms_per_disc=3)
    dt = time.time() - t0
    assert got and 1 < got[0] < n and n % got[0] == 0, n
    rows.append([f"{n:,}", f"{got[0]:,}", f"{n // got[0]:,}", got[1], f"{dt:.2f}"])
report.table(["n", "factor", "cofactor", "multiplier k that worked", "seconds"], rows)
report.p()

# --------------------------------------------------------------- density
report.p("## The density measurement")
report.p()
report.p("How many distinct values can the governing quantity take, at fixed input?")
report.p()
rows = []
for n in (143, 323, 667, 1147, 2021, 5183):
    hs = []
    for k in range(1, 81):
        d = -k * n
        if d % 4 in (0, 1):
            hs.append(class_number_bruteforce(d))
    rows.append([f"{n:,}", len(hs), len(set(hs)), f"{len(set(hs)) / len(hs):.0%}",
                 f"{min(hs)}..{max(hs)}"])
report.table(["n", "multipliers tried", "distinct h(-kn)", "distinct share",
              "range of h"], rows)
report.p()
p = 1009
rows = []
for degree in (2, 3, 4):
    seen = set()
    for _ in range(300):
        f = [rng.randrange(p) for _ in range(degree)] + [1]
        o = ring_unit_order(f, p)
        if o is not None:
            seen.add(o)
    rows.append([f"ring unit group, degree {degree}", len(seen),
                 f"partitions({degree}) = {partitions(degree)}"])
ec = set()
for _ in range(300):
    o = elliptic_order(rng.randrange(p), rng.randrange(p), p)
    if o is not None:
        ec.add(o)
rows.append(["elliptic curve", len(ec), f"Hasse width 4*sqrt(p) = {4 * isqrt(p)}"])
report.table(["construction (at fixed p = 1009)", "distinct orders seen",
              "theoretical reach"], rows)
report.p("Class numbers are dense: nearly every multiplier gives a fresh number, "
         "spread over a wide range. The ring is capped at `partitions(d)` -- two, "
         "three, five -- no matter how many polynomials are tried. **Round 4's "
         "requirement is met.**")
report.p()

# --------------------------------------------------------------- the ceiling
report.p("## And the ceiling that meeting it runs into")
report.p()
report.p("Meeting the density requirement does not give polynomial time. It gives "
         "`L[1/2]` -- the same class as ECM. The reason is structural, and it is "
         "the real conclusion of this project.")
report.p()
report.p("Every method in the ladder below is a **smoothness lottery**: build a "
         "group whose order is some integer attached to `p`, then hope that "
         "integer is `B`-smooth.")
report.p()
report.p("| the set the parameter ranges over | example | what you get |")
report.p("|---|---|---|")
report.p("| a single value | classical Pascal row: period `p` | `Theta(p)` -- trial division |")
report.p("| divisor lattice of `p-1` | q-Pascal row; Pollard `p-1` | fast only when `p-1` is smooth |")
report.p("| one value per `(p, d)` | norm-one subgroup: `Phi_d(p)` | Williams `p+1` and relatives |")
report.p("| `partitions(d)` values | ring unit groups | a constant, independent of `p` |")
report.p("| **dense interval near `p`** | ECM: `p + 1 - t` | `L[1/2]` |")
report.p("| **dense, near `sqrt(kn)`** | class groups: `h(-kn)` | `L[1/2]` |")
report.p()
report.p("Density is what separates the top two rows from the rest. But it is not "
         "enough, because the probability that a number of size `p` is `B`-smooth "
         "is itself governed by the Dickman function, and optimising `B` against "
         "the cost of sampling gives `L[1/2]` no matter how good the sampling is. "
         "**Denser sampling cannot beat the smoothness density itself.**")
report.p()
report.p("This is why the number field sieve is faster: `L[1/3]` comes not from "
         "sampling better but from making the numbers *smaller* -- testing "
         "smoothness of algebraic norms of size `L`-ish rather than of size `p`. "
         "It changes what is being tested, not how the test is sampled.")
report.p()
report.p("## Where this leaves the project")
report.p()
report.p("Five rounds have walked a ladder of requirements, each one answered and "
         "each answer revealing the next constraint:")
report.p()
report.p("1. *An asymmetric statistic* -- supplied by the norm-one subgroup (round 3).")
report.p("2. *A parameter that varies* -- supplied by the q-deformation (round 4).")
report.p("3. *Variation over a dense set* -- supplied by class groups (round 5).")
report.p()
report.p("Each was necessary. None was sufficient, and the reason is now visible: "
         "all three are refinements of the *same* mechanism, sampling for "
         "smoothness, whose ceiling is `L[1/2]` regardless of how well it is "
         "refined. The next requirement is therefore not a better group but a "
         "different mechanism -- and that is where this line of attack, which "
         "began with a conjecture about remainders in Pascal's triangle, honestly "
         "runs out.")
report.p()
report.p("What the Pascal framing did produce is exact and worth keeping: the "
         "first non-zero entry of row `n` mod `n` sits at `spf(n)` and equals "
         "`n/spf(n)`; the factors live in the second n-adic digit at constant "
         "density; and reading that digit is provably as hard as factoring.")
report.write()
