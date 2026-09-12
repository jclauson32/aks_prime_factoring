"""Round 38: the binary search for `p`, run to the end.

The obvious idea, and the one that keeps coming back: `p` is a number, its bits
can be searched one at a time, and `N = p q` should constrain each bit.  Lift
`x y == N (mod 2^k)` from `k` to `k + 1` and you have a search tree.

Round 24 already showed the branches never die: every odd residue of `p` has a
partner `q = N/p mod 2^k`, so all `2^(k-1)` of them survive.  This round builds
the tree, measures every prune available from `N` alone, and then finishes the
job with Coppersmith: the tree does terminate in a working
factoring algorithm, and the algorithm costs `N^(1/4)`, which is where the size
mechanism already was (round 34).  The tree is not a shortcut, it is the same
wall reached from the 2-adic side.
"""

import random
import time
from math import log2

from _common import Report

from aksfactor.arith import is_prime
from aksfactor.lattice import congruence_modulus_needed, factor_with_congruence

report = Report(
    "exp44_padic_tree",
    "The 2-adic search tree, and what it costs to finish it",
    "Round 38. Bitwise lifting has no pruning; with Coppersmith it is an "
    "`N^(1/4)` algorithm.",
)

rng = random.Random(44)


def randprime(bits):
    while True:
        v = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if is_prime(v):
            return v


# ---------------------------------------------------------------- section 1
report.p("## The tree")
report.p()
report.p("A node at depth `k` is a pair `(x, y)` of odd residues mod `2^k` with "
         "`x y == N`. Its children are the four ways to add a bit to each; the "
         "congruence keeps some of them. Round 24 counted the survivors from the "
         "residue side -- every odd `r` is a legal `p mod 2^k` -- and the counts "
         "below enumerate the tree itself.")
report.p()

n_demo = randprime(14) * randprime(14)
level = [(1, 1)]
rows = []
for k in range(1, 16):
    mod = 1 << (k + 1)
    children = []
    for x, y in level:
        for a in (0, 1):
            for b in (0, 1):
                nx, ny = x + (a << k), y + (b << k)
                if nx * ny % mod == n_demo % mod:
                    children.append((nx, ny))
    branch = len(children) / len(level)
    odd_x = len({x for x, _ in children})
    rows.append([k + 1, len(children), f"{branch:.2f}", odd_x, 1 << k])
    level = children
report.table(["depth k", "nodes", "children per node", "distinct `x` values",
              "odd residues mod 2^k"], rows)
report.p("Two of the four children survive at every node, so the tree doubles "
         "with depth, and the distinct `x` values are *all* the odd residues mod "
         "`2^k`. That is the whole difficulty in one line: for every odd `r` "
         "there is a `y` with `r y == N (mod 2^k)`, so the congruence says "
         "nothing whatever about `p` on its own. There is no branch to cut.")
report.p()

# ---------------------------------------------------------------- section 2
report.p("## Every prune available from `N` alone")
report.p()
report.p("The prune a bitwise search would actually use is the product bound: a "
         "node whose literal `x y` already exceeds `N` cannot be completed "
         "downward. Below, the share of nodes it removes, and the share removed "
         "by testing `x | N` (which is trial division wearing a hat).")
report.p()

bits_n = n_demo.bit_length()
level = [(1, 1)]
rows = []
for k in range(1, 19):
    mod = 1 << (k + 1)
    children = []
    for x, y in level:
        for a in (0, 1):
            for b in (0, 1):
                nx, ny = x + (a << k), y + (b << k)
                if nx * ny % mod == n_demo % mod:
                    children.append((nx, ny))
    kept_bound = [(x, y) for x, y in children if x * y <= n_demo]
    divides = [(x, y) for x, y in children if x > 1 and n_demo % x == 0]
    rows.append([k + 1, len(children),
                 f"{1 - len(kept_bound) / len(children):.0%}",
                 len(divides)])
    level = children
report.table(["depth k", "nodes", "removed by the product bound", "nodes with `x | N`"],
             rows)
report.p(f"The product bound removes nothing until the two halves together "
         f"exceed `N`, which happens at depth `log2(N)/2` = {bits_n // 2} here -- "
         f"and the next section needs the tree only to depth `log2(N)/4` = "
         f"{bits_n // 4}. In the range that matters the prune is inactive. The "
         f"`x | N` column is the other way to end the search: it finds `p` once "
         f"the tree is deep enough to contain it, having enumerated about `p` "
         f"nodes to get there, which is trial division with extra steps.")
report.p()

# ---------------------------------------------------------------- section 3
report.p("## Finishing the tree with Coppersmith")
report.p()
report.p("A node at depth `t` is exactly a claim about `p mod 2^t`, and round 12 "
         "showed Coppersmith finishes from such a claim once `2^t` is about "
         "`N^(1/4)`. So the tree does not have to be walked to the bottom: it has "
         "to be walked to depth `log2(N)/4`, and every leaf handed to a lattice. "
         "First, how much depth the lattice really needs, measured on instances "
         "where `p` is known.")
report.p()

cases = {}
rows = []
for bits in (28, 34, 40):
    pp, qq = randprime(bits // 2), randprime(bits - bits // 2)
    n = pp * qq
    base = congruence_modulus_needed(n).bit_length()
    found = None
    for extra in range(0, 8):
        if factor_with_congruence(n, pp % (1 << (base + extra)), 1 << (base + extra)):
            found = extra
            break
    cases[bits] = (n, pp, base, found)
    rows.append([bits, base, found if found is not None else "> 7",
                 base + found if found is not None else "--"])
report.table(["bits of N", "`log2(N^(1/4))`", "extra bits this instance needed",
              "depth used"], rows)

extras = [c[3] for c in cases.values() if c[3] is not None]
report.p(f"Coppersmith needs a little more than the bare quarter, and how much "
         f"more is a property of the instance rather than of its size -- "
         f"{', '.join(str(v) for v in extras)} bits for these three, because the "
         f"bound it has to beat is asymptotic and 28 bits is not asymptotic. Each "
         f"row below walks the same instance whose depth was just measured. "
         f"Measuring it needed `p`, which the search does not have: a real run "
         f"would deepen the tree one bit at a time until a leaf worked, which "
         f"costs about twice the final depth's leaves and changes nothing below.")
report.p()

rows = []
for bits, (n, pp, base, extra) in cases.items():
    t = base + (extra if extra is not None else 7)
    mod = 1 << t
    leaves = mod // 2
    if leaves <= 20000:                 # small enough to run to the end
        t0 = time.perf_counter()
        tried, got = 0, None
        for r in range(1, mod, 2):
            tried += 1
            out = factor_with_congruence(n, r, mod)
            if out is not None:
                got = out
                break
        dt = time.perf_counter() - t0
        assert got is not None and got[0] * got[1] == n, f"walk failed at {bits}"
        rows.append([bits, t, f"{leaves:,}", f"{log2(leaves) / log2(n):.3f}",
                     f"{tried:,}", f"{dt:.2f} s", "run to a split"])
    else:                               # time the leaf, project the walk
        probe = [rng.randrange(1, mod, 2) for _ in range(200)]
        t0 = time.perf_counter()
        for r in probe:
            factor_with_congruence(n, r, mod)
        per_leaf = (time.perf_counter() - t0) / len(probe)
        assert factor_with_congruence(n, pp % mod, mod) is not None
        rows.append([bits, t, f"{leaves:,}", f"{log2(leaves) / log2(n):.3f}",
                     f"{leaves // 2:,} expected",
                     f"{per_leaf * leaves / 2:.0f} s", f"{per_leaf * 1000:.1f} ms/leaf"])
report.table(["bits of N", "depth t", "leaves in the tree",
              "`log(leaves)/log(N)`", "leaves tried", "time", "how"], rows)
def joined(items):
    words = [f"{b} bits" for b in items]
    if len(words) < 2:
        return "".join(words)
    return ", ".join(words[:-1]) + " and " + words[-1]


walked = [r[0] for r in rows if r[6] == "run to a split"]
projected = [r[0] for r in rows if r[6] != "run to a split"]
report.p(f"The {'row' if len(walked) == 1 else 'rows'} at {joined(walked)} "
         f"{'was' if len(walked) == 1 else 'were'} walked to the end: every odd "
         f"residue mod `2^t` handed to the lattice until one of them factored "
         f"the number. "
         + (f"{joined(projected)} "
            f"{'was' if len(projected) == 1 else 'were'} not, because a tree "
            f"that size takes hours; there the cost of one leaf was measured "
            f"over 200 of them and multiplied by half the tree, and the lattice "
            f"was checked separately on the true residue so the row is not "
            f"projecting a walk that would fail. " if projected else "")
         + "A full walk stops at whichever of `p` and `q` has the earlier "
           "residue, so it runs about four times short of the whole tree on "
           "average.")
report.p()
exponents = [float(r[3]) for r in rows]
report.p(f"The algorithm is real -- every row factored its number -- and its "
         f"exponent is the one the size mechanism already had. The "
         f"`log(leaves)/log(N)` column runs {exponents[0]:.3f} to "
         f"{exponents[-1]:.3f}: the tree is `N^(1/4)` leaves plus the lattice's "
         f"fixed margin, and the margin is what keeps the measured column above "
         f"0.25 at sizes this small. The lattice itself is polynomial, so the "
         f"whole thing is `N^(1/4)` times a polynomial, and round 34's table "
         f"lists four other ways to pay exactly that.")
report.p()

# ---------------------------------------------------------------- section 4
report.p("## What would have to change")
report.p()
report.p("The tree is prunable the moment anything else is known. One fact -- the "
         "true `p mod 2^t` -- collapses it from `N^(1/4)` leaves to one, which is "
         "the measurement in the first table of section 3. That is the shape of "
         "Heninger and Shacham's attack on RSA keys with corrupted bits: their "
         "tree is the same tree, and what prunes it is the redundancy between "
         "`p`, `q`, `d`, `d_p` and `d_q`, all of which have to agree. Given only "
         "`N` there is no second equation, so nothing disagrees with anything, "
         "and every branch survives.")
report.p()
report.p("So the honest answer to 'can the bits of `p` be searched one at a "
         "time': yes, and it terminates, and it costs `N^(1/4)`. A binary search "
         "needs a comparison that eliminates half the space; the 2-adic "
         "congruence eliminates nothing, because it is satisfiable for every odd "
         "residue. What makes a search fast is not that the space is searchable "
         "in order -- it is that the order tells you which half to drop.")
report.write()
